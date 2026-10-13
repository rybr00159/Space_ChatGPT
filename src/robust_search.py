"""Gap-edge-excluded joint search with two robust median widths.

Exclusion reduces completeness at boundaries. This alternate preprocessing tests
whether rankings depend on boundary excursions; it cannot replace pixel vetting.
"""
import argparse,json
import numpy as np
from search import ROOT,CACHE,flatten,read_curve,known_for,mask_known,scan,phase
from joint_search import grouped
import matplotlib.pyplot as plt

def edge_mask(t,cut=.25,gap=.2):
    keep=np.ones(len(t),dtype=bool)
    edges=np.r_[0,np.flatnonzero(np.diff(t)>gap)+1,len(t)]
    for a,b in zip(edges[:-1],edges[1:]):
        keep[a:b]=(t[a:b]>=t[a]+cut)&(t[a:b]<=t[b-1]-cut)
    return keep

def run(group,width,cut=.25,maxperiod=25):
    tic=group[0]['tic'];sectors=[r['sector'] for r in group];chunks=[]
    for rec in group:
        t,f,e,h=read_curve(rec['local_path']); tb,fb,eb=flatten(t,f,e,width)
        # Boundaries found on quality-filtered raw samples, rather than after flare clipping.
        rawkeep=edge_mask(t,cut);edges=np.r_[0,np.flatnonzero(np.diff(t)>.2)+1,len(t)]
        good=np.zeros(len(tb),bool)
        for a,b in zip(edges[:-1],edges[1:]):good|=(tb>=t[a]+cut)&(tb<=t[b-1]-cut)
        good &= mask_known(tb,known_for(tic))
        chunks.append((tb[good],fb[good],eb[good]))
    t,f,e=(np.concatenate([c[j] for c in chunks]) for j in range(3))
    label=f"tic{tic}_s{'-'.join(map(str,sectors))}_w{width:g}_edge{cut:g}"
    p=ROOT/'results'/'robust';p.mkdir(parents=True,exist_ok=True)
    if (p/f'{label}.json').exists():return
    peaks,r=scan(t,f,e,maximum=maxperiod,n_peaks=5)
    intermediate=CACHE/'intermediate'/f'{label}.npz'
    np.savez_compressed(intermediate,t=t,f=f,e=e)
    result={'tic':tic,'sectors':sectors,'width':width,'edge_cut_days':cut,'n_searched':len(t),'t_min':float(t.min()),'t_max':float(t.max()),'peaks':peaks,'intermediate':str(intermediate)}
    (p/f'{label}.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    if peaks:
        best=peaks[0];fig,ax=plt.subplots(3,1,figsize=(10,9),layout='constrained')
        ax[0].plot(t,f,'.',ms=1);ax[0].set(xlabel='BTJD',ylabel='Flux',title=label)
        ax[1].plot(r.period,r.power,lw=.5);ax[1].set(xlabel='Period (days)',ylabel='BLS diagonal S/N')
        ax[2].plot(phase(t,best['period'],best['t0_btjd']),f,'.',ms=1);ax[2].set(xlabel='Phase time (days)',ylabel='Flux',xlim=(-.4,.4),title=f"P={best['period']:.5f}; SNR={best['snr']:.1f}; {best['n_events']} events")
        fig.savefig(p/f'{label}.png',dpi=120);plt.close(fig)
    print(label,[(round(p['period'],4),round(p['snr'],1),p['n_events']) for p in peaks],flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--tic',type=int,nargs='+');p.add_argument('--width',type=float,nargs='+',default=[1,2]);p.add_argument('--edge-cut',type=float,default=.25);p.add_argument('--max-period',type=float,default=25)
    a=p.parse_args();records=json.loads((ROOT/'data'/'observations.json').read_text())
    for tic in sorted({r['tic'] for r in records}):
        if a.tic and tic not in a.tic:continue
        rows=[r for r in records if r['tic']==tic];groups=grouped(rows)
        used={r['sector'] for g in groups for r in g}
        groups += [[r] for r in rows if r['sector'] not in used]
        for group in groups:
            for width in a.width:run(group,width,a.edge_cut,a.max_period)
if __name__=='__main__':main()
