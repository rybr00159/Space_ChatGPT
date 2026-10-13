"""Per-sector BLS screening. Scores are ranking statistics, not false-alarm probabilities."""
from pathlib import Path
import argparse, json, os
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get('TESS_DATA_DIR','/workspace/tess-research-data'))
DURATIONS = np.array([.045,.065,.09,.13,.18,.23])

def phase(t,p,t0): return (t-t0+.5*p)%p-.5*p

def read_curve(path, column='PDCSAP_FLUX'):
    with fits.open(path) as hdul:
        d=hdul[1].data; h=dict(hdul[0].header)
        timing=hdul[1].header
        if timing.get('TIMESYS') != 'TDB' or abs(timing.get('BJDREFI',2457000)+timing.get('BJDREFF',0)-2457000)>1e-8:
            raise ValueError('Unexpected time system/reference; do not assume BTJD')
        t=np.array(d['TIME'],dtype=float); f=np.array(d[column],dtype=float); e=np.array(d[column+'_ERR'],dtype=float); q=np.array(d['QUALITY'])
    keep=np.isfinite(t)&np.isfinite(f)&np.isfinite(e)&(e>0)&(q==0)&(f>0)
    t,f,e=t[keep],f[keep],e[keep]
    order=np.argsort(t); t,f,e=t[order],f[order],e[order]
    return t,f,e,h

def flatten(t,f,e,width=1.0,bin_minutes=10):
    """Gap-aware robust rolling median, then inverse-variance bins; no low clipping."""
    fs=np.empty_like(f); es=np.empty_like(e)
    edges=np.r_[0,np.flatnonzero(np.diff(t)>.2)+1,len(t)]
    for lo,hi in zip(edges[:-1],edges[1:]):
        dt=np.median(np.diff(t[lo:hi])) if hi-lo>2 else .0014
        n=max(3,int(width/dt)//2*2+1)
        trend=median_filter(f[lo:hi],size=n,mode='reflect')
        fs[lo:hi]=f[lo:hi]/trend; es[lo:hi]=e[lo:hi]/trend
    # Drop only high flare/cosmic-ray points; leave deep negative events for vetting.
    scatter=1.4826*np.median(np.abs(fs-np.median(fs)))
    keep=fs<1+max(7*scatter,.001)
    t,fs,es=t[keep],fs[keep],es[keep]
    idx=np.floor((t-t.min())/(bin_minutes/1440)).astype(int)
    uniq,inv=np.unique(idx,return_inverse=True)
    w=1/es**2; sw=np.bincount(inv,weights=w)
    tb=np.bincount(inv,weights=w*t)/sw
    fb=np.bincount(inv,weights=w*fs)/sw
    eb=1/np.sqrt(sw)
    # Conservative diagonal errors >= robust observed point scatter.
    sigma=1.4826*np.median(np.abs(fb-np.median(fb)))
    eb=eb*max(1,sigma/np.median(eb))
    return tb,fb,eb

def known_for(tic):
    p=ROOT/'data'/'known_planets.csv'
    if not p.exists(): return []
    table=pd.read_csv(p)
    return table[table['tic'].astype(int)==int(tic)].to_dict('records')

def mask_known(t,known):
    mask=np.ones(len(t),bool)
    for k in known:
        p=k.get('period'); t0=k.get('epoch_bjd'); dur=k.get('duration_hours')
        if pd.isna(p) or pd.isna(t0): continue
        if pd.isna(dur) or dur<=0: dur=3.0
        mask &= np.abs(phase(t,p,t0-2457000))>max(dur/24*.85,.05)
    return mask

def scan(t,f,e,minimum=.7,maximum=13.0,n_peaks=3):
    maximum=min(maximum,(t.max()-t.min())/2)
    if maximum<=minimum: return [],None
    bls=BoxLeastSquares(t,f,e)
    r=bls.autopower(DURATIONS,minimum_period=minimum,maximum_period=maximum,
                    frequency_factor=.7,oversample=10,objective='snr')
    score=np.array(r.power); chosen=[]; available=np.ones(len(score),bool)
    for _ in range(n_peaks):
        j=int(np.argmax(np.where(available,score,-np.inf))); p=float(r.period[j]);t0=float(r.transit_time[j]);duration=float(r.duration[j])
        intrans=np.abs(phase(t,p,t0))<duration/2
        cycles=np.round((t[intrans]-t0)/p).astype(int)
        centers=[]
        for cycle in np.unique(cycles):
            ii=intrans&(np.round((t-t0)/p).astype(int)==cycle)
            oot=(np.abs(t-(t0+cycle*p))<duration*2)&~intrans
            if ii.sum()<2 or oot.sum()<3: continue
            centers.append({'cycle':int(cycle),'time':float(t0+cycle*p),'n':int(ii.sum()),
                            'depth':float(np.median(f[oot])-np.mean(f[ii]))})
        stats=bls.compute_stats(p,duration,t0)
        odd,even=stats['depth_odd'],stats['depth_even']
        chosen.append({'period':p,'t0_btjd':t0,'duration_days':duration,'depth':float(r.depth[j]),
                       'depth_err':float(r.depth_err[j]),'snr':float(r.depth_snr[j]),
                       'n_events':len(centers),'events':centers,
                       'odd_even_sigma':float(abs(odd[0]-even[0])/np.hypot(odd[1],even[1])) if np.hypot(odd[1],even[1])>0 and len({v['cycle']%2 for v in centers})==2 else None,
                       'sde_robust':float((score[j]-np.median(score))/(1.4826*np.median(np.abs(score-np.median(score))))),
                       'period_grid_count':len(score)})
        # Suppress very close aliases only, preserving harmonics for explicit review.
        available &= np.abs(1/r.period-1/p)>.02/(t.max()-t.min())
        available &= np.abs(r.period-p)/p>.01
    return chosen,r

def process(record,width,mask=True,column='PDCSAP_FLUX'):
    t,f,e,h=read_curve(record['local_path'],column)
    tb,fb,eb=flatten(t,f,e,width)
    known=known_for(record['tic']); keep=mask_known(tb,known) if mask else np.ones(len(tb),bool)
    # Keep raw and flattened intermediate values in persistent cache.
    label=f"tic{record['tic']}_s{record['sector']:04}_w{width:g}_{column}_{'residual' if mask else 'control'}"
    inter=CACHE/'intermediate'/f'{label}.npz'; inter.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(inter,t=tb,f=fb,e=eb,mask=keep)
    peaks,r=scan(tb[keep],fb[keep],eb[keep])
    result={'tic':record['tic'],'sector':record['sector'],'width':width,'column':column,'known_masked':mask,
            'known_planets':known,'data_uri':record['data_uri'],'n_raw':len(t),'n_binned':len(tb),'n_searched':int(keep.sum()),
            't_min':float(tb.min()),'t_max':float(tb.max()),'crowdsap':h.get('CROWDSAP'),
            'intermediate':str(inter),'peaks':peaks}
    # CROWDSAP typically in extension header; read explicitly.
    with fits.open(record['local_path']) as hdul: result['crowdsap']=hdul[1].header.get('CROWDSAP')
    out=ROOT/'results'/'scans';out.mkdir(parents=True,exist_ok=True)
    (out/f'{label}.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    if peaks:
        best=peaks[0];fig,axes=plt.subplots(3,1,figsize=(10,9),layout='constrained')
        axes[0].plot(tb,fb,'.',ms=1,color='.6');axes[0].plot(tb[keep],fb[keep],'.',ms=1)
        axes[0].set(xlabel='BTJD (BJD_TDB - 2457000)',ylabel='Normalized flux',title=f'TIC {record["tic"]}, sector {record["sector"]}, width {width:g} d')
        axes[1].plot(r.period,r.power,lw=.5);axes[1].set(xlabel='Period (days)',ylabel='BLS S/N (diagonal errors)')
        ph=phase(tb[keep],best['period'],best['t0_btjd']);axes[2].plot(ph,fb[keep],'.',ms=1)
        axes[2].set(xlabel='Phase time (days)',ylabel='Normalized flux',xlim=(-best['duration_days']*3,best['duration_days']*3),
                    title=f"P={best['period']:.5f} d; depth={best['depth']*1e6:.0f} ppm; S/N={best['snr']:.1f}; {best['n_events']} events")
        fig.savefig(out/f'{label}.png',dpi=120);plt.close(fig)
    print(label,[(round(p['period'],4),round(p['snr'],1),p['n_events']) for p in peaks],flush=True)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--tic',type=int,nargs='+');p.add_argument('--width',type=float,nargs='+',default=[1,2]);p.add_argument('--control',action='store_true');p.add_argument('--column',default='PDCSAP_FLUX');p.add_argument('--force',action='store_true')
    a=p.parse_args();records=json.loads((ROOT/'data'/'observations.json').read_text())
    for rec in records:
        if a.tic and rec['tic'] not in a.tic: continue
        for width in a.width:
            label=f"tic{rec['tic']}_s{rec['sector']:04}_w{width:g}_{a.column}_{'control' if a.control else 'residual'}"
            if (ROOT/'results'/'scans'/f'{label}.json').exists() and not a.force: continue
            process(rec,width,not a.control,a.column)
if __name__=='__main__':main()
