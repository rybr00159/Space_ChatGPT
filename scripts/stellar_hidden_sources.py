"""Independent coherence and morphology audit of potential faint stellar blends.

Uses retained, hash-verified QLP fluxes. It does not assume that a large-radius
catalog host produces a box-search peak. Figures retain full orbital phase and
raw short-cadence structure, which distinguish giant troughs from eclipses.
"""
from pathlib import Path
import sys, os, json, hashlib
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.optimize import minimize_scalar
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from broad_vet import read_source, phase, local_quadratic, integrated_trapezoid
OUT=ROOT/'candidates/stellar_systems/hidden_sources'
CONFIG={260352037:(.8807,.05),406280158:(10.1258553,.16),271697430:(3.0593723,.16)}

def load(tic):
    qa=json.loads((ROOT/f'candidates/broad_s1/giant_survivors/tic{tic}_qa.json').read_text())
    curves=[]
    for rec in qa['independent_sector_qa']:
        d,checks=read_source(rec['source'])
        use=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.isfinite(d['SAP_FLUX'])&(d['SAP_FLUX']>0)
        t=d['TIME'][use];f=d['SAP_FLUX'][use];f=f/np.median(f)
        e=d['KSPSAP_FLUX_ERR'][use]/np.median(d['SAP_FLUX'][use])
        dt=float(np.median(np.diff(t)))
        n=max(3,int(2/dt)//2*2+1)
        trend=median_filter(f,size=n,mode='nearest')
        curves.append(dict(sector=rec['sector'],t=t,f=f,e=e,dt=dt,flat=f/trend,source=rec['source'],checks=checks,prior=rec['independent_local_search']))
    return curves,qa

def bins(x,y,n=100):
    edges=np.linspace(-.5,.5,n+1);which=np.digitize(x,edges)-1
    center=(edges[1:]+edges[:-1])/2
    med=np.array([np.median(y[which==i]) if np.sum(which==i)>0 else np.nan for i in range(n)])
    return center,med

def profile_score(t,y,p,epoch,duration,dt):
    xx=phase(t,p,epoch);shape=integrated_trapezoid(xx,duration,.015,dt)
    cols=np.column_stack([np.ones(len(t)),shape]);coef=np.linalg.lstsq(cols,y,rcond=None)[0]
    return float(np.sum((y-cols@coef)**2)),float(-coef[-1])

def analyze(tic):
    curves,qa=load(tic);p0,duration=CONFIG[tic];e0=curves[0]['prior']['t0_btjd']
    # First two consecutive sectors set the template period and epoch; later
    # sectors test that fixed prediction, without fitting their phases.
    train=[c for c in curves if c['sector']<=2]
    if len(train)<2: train=curves[:1]
    t=np.concatenate([c['t'] for c in train]);y=np.concatenate([c['flat']-1 for c in train]);dt=np.median([c['dt'] for c in train])
    # Grid resolution limits phase smearing to <1/20 template duration.
    span=t.max()-t.min();half=.002*p0
    step=duration*p0/(20*span)
    ps=np.linspace(p0-half,p0+half,int(2*half/step)+3)
    es=np.linspace(e0-duration,e0+duration,51)
    best=(np.inf,None,None,None)
    for p in ps:
        for ep in es:
            s,dep=profile_score(t,y,p,ep,duration,dt)
            if dep>0 and s<best[0]:best=(s,p,ep,dep)
    _,p,epoch,depth=best
    rows=[];spectra=[]
    for c in curves:
        fixed,data=local_quadratic(c['t'],c['f'],c['e'],p,epoch,duration)
        prior=c['prior'];floating,fdata=local_quadratic(c['t'],c['f'],c['e'],prior['period'],prior['t0_btjd'],prior['duration_days'])
        freq=np.linspace(.08,8.,16000)
        power=LombScargle(c['t'],c['flat']-1).power(freq)
        top=np.argsort(power)[-8:][::-1]
        # Harmonic amplitudes, normalized raw flux, sector-specific constant.
        xx=2*np.pi*(c['t']-epoch)/p
        cols=[np.ones(len(xx))]
        for k in range(1,9):cols.extend([np.cos(k*xx),np.sin(k*xx)])
        coef=np.linalg.lstsq(np.column_stack(cols),c['flat']-1,rcond=None)[0]
        amps=[float(np.hypot(coef[2*k-1],coef[2*k])*1e6) for k in range(1,9)]
        # Fixed ephemeris fold and fitted-width validation are distinct from
        # selecting the strongest phase independently in each sector.
        rows.append(dict(sector=c['sector'],n=len(c['t']),fixed_ephemeris=fixed,independently_optimized_ephemeris=floating,phase_drift_of_prior_days=float(phase(np.array([prior['t0_btjd']]),p,epoch)[0]),top_frequency_per_day=freq[top].tolist(),top_lomb_power=power[top].tolist(),harmonic_amplitudes_ppm=amps,source=c['source'],checksum_status=c['checks']))
        spectra.append((freq,power))
    selected=sorted(set([0,min(1,len(curves)-1),len(curves)//2,len(curves)-1]))
    fig,axes=plt.subplots(len(selected),3,figsize=(14,3*len(selected)),layout='constrained')
    for ax,c,spec in zip(axes, [curves[i] for i in selected],[spectra[i] for i in selected]):
        pha=phase(c['t'],p,epoch)/p;x,b=bins(pha,c['flat'])
        ax[0].plot(pha,c['flat'],'.',ms=1.5,alpha=.25);ax[0].plot(x,b,'r.-',ms=3);ax[0].set(xlabel='Phase: fixed early-sector ephemeris',ylabel=f'S{c["sector"]}: 2-day normalized SAP')
        mid=c['t'].min()+5;use=abs(c['t']-mid)<1.8
        ax[1].plot(c['t'][use],c['f'][use],'.-',ms=2);ax[1].set(xlabel='BTJD',ylabel='Raw SAP / sector median')
        ax[2].plot(spec[0],spec[1]);ax[2].axvline(1/p,color='r');ax[2].set(xlabel='Frequency / day',ylabel='Lomb–Scargle power',xlim=(.1,5))
    fig.suptitle(f'TIC {tic}: early-sector P={p:.9f} d; morphology and fixed-phase holdout')
    fig.savefig(OUT/f'tic{tic}_coherence.png',dpi=150);plt.close(fig)
    result=dict(tic=tic,checked_utc=datetime.now(timezone.utc).isoformat(),training_sectors=[c['sector'] for c in train],period_days=float(p),epoch_btjd=float(epoch),duration_days=duration,early_template_depth_ppm=depth*1e6,sectors=rows,limitations=['Selected peaks have no global false-alarm calibration.','A two-day running median can reshape broad stellar troughs; raw SAP and harmonic information must be inspected.','Training grid may not establish a unique multi-year cycle count. No source is assigned before pixel localization.'],code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/f'tic{tic}_coherence.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(tic,'P',p,'epoch',epoch,'sector fixed-depthSNR',[(r['sector'],round(r['fixed_ephemeris']['snr'],2),round(r['phase_drift_of_prior_days'],3)) for r in rows],flush=True)
    return result

def pixel_audit(tic=260352037, sector=28, bootstrap=200, event_phase=0.):
    from astropy.io import fits
    from astropy.wcs import WCS
    from archive import verify_fits_checksums
    from pixel_vet_49617264 import gaussfit, combine, safe
    import io, urllib.parse, urllib.request, zipfile
    ra,dec=94.787911339,-57.078220725
    cache=Path('/workspace/tess-research-data/stellar_systems/hidden_sources')/f'tic{tic}'
    cache.mkdir(parents=True,exist_ok=True)
    url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urllib.parse.urlencode(dict(ra=ra,dec=dec,x=11,y=11,units='px',sector=sector))
    files=list(cache.glob(f'tess-s{sector:04}-*.fits'))
    if not files:
        with urllib.request.urlopen(url,timeout=60) as response:blob=response.read(60000001)
        if len(blob)>60000000:raise RuntimeError('60 MB bound exceeded')
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            for name in z.namelist():
                if name.endswith('.fits'):(cache/Path(name).name).write_bytes(z.read(name))
        files=list(cache.glob(f'tess-s{sector:04}-*.fits'))
    if len(files)!=1:raise RuntimeError('Expected one cube')
    path=files[0];out=OUT/f'tic{tic}_pixels'/f'sector{sector}'
    if event_phase:out=out/'phase_half'
    out.mkdir(parents=True,exist_ok=True)
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    old=out/'raw_manifest.json'
    if old.exists() and json.loads(old.read_text())['sha256']!=digest:raise RuntimeError('Original byte hash changed')
    with fits.open(path) as h:
        checks=verify_fits_checksums(path,h);w=WCS(h[2].header);header=dict(h[0].header);lh=dict(h[1].header)
        d=h[1].data;t=np.array(d['TIME']);q=np.array(d['QUALITY']);f=np.array(d['FLUX'],float);err=np.array(d['FLUX_ERR'],float);ffi=np.array(d['FFI_FILE']).astype(str)
    dt=float(lh['TIMEDEL']);host=w.all_world2pix([[ra,dec]],0)[0]
    column=lh['1CRV5P']+(host[0]+1-lh['1CRP5P'])*lh['1CDL5P'];row=lh['2CRV5P']+(host[1]+1-lh['2CRP5P'])*lh['2CDL5P']
    if not (44<=column<2092 and 0<=row<2048):raise RuntimeError('Not science silicon')
    good=(q==0)&np.isfinite(t)&np.all(np.isfinite(f)&np.isfinite(err),axis=(1,2))
    t,f,err=t[good],f[good],err[good]
    reference=np.median(f,axis=0);yy,xx=np.indices(reference.shape);blank=reference<np.percentile(reference,20)
    bg=np.median(f[:,blank],axis=1);f=f-bg[:,None,None];reference=reference-np.median(reference[blank])
    # A reference gradient provides pointing nuisance columns; all stars inside
    # the stamp contribute, so pointing is a sensitivity test, not a calibrated
    # registration model or physically exact PRF correction.
    gy,gx=np.gradient(reference);reg=np.column_stack([reference.ravel(),gx.ravel(),gy.ravel(),np.ones(reference.size)])
    cc=np.linalg.lstsq(reg,f.reshape(len(t),-1).T,rcond=None)[0]
    shifts=(-cc[1]/cc[0],-cc[2]/cc[0])
    eph=json.loads((OUT/f'tic{tic}_coherence.json').read_text());p=eph['period_days'];epoch=eph['epoch_btjd']+event_phase*p;duration=.065
    modes={};sets={}
    for mode in ['quadratic','background','pointing']:
        events=[]
        for cycle in range(int(np.ceil((t.min()-epoch)/p)),int(np.floor((t.max()-epoch)/p))+1):
            center=epoch+cycle*p;use=abs(t-center)<.24;tx=t[use]-center
            shape=integrated_trapezoid(tx,duration,.022,dt)
            min_in=max(2,int(.55*duration/dt))
            min_out=max(5,int(.07/dt))
            if np.sum(shape>.05)<min_in or np.sum((shape==0)&(tx<0))<min_out or np.sum((shape==0)&(tx>0))<min_out:continue
            if np.max(np.diff(tx))>.035:continue
            cols=[np.ones(len(tx)),tx,tx**2]
            if mode in ['background','pointing']:
                nuisances=[bg[use]]
                if mode=='pointing':nuisances.extend([shift[use] for shift in shifts])
                for cov in nuisances:
                    cov=cov-np.mean(cov)
                    if np.std(cov)>1e-8:cols.append(cov/np.std(cov))
            design=np.column_stack(cols+[-shape]);inv=np.linalg.pinv(design);ff=f[use].reshape(len(tx),-1);coef=inv@ff;res=ff-design@coef
            noise=np.maximum(np.std(res,axis=0,ddof=design.shape[1]),np.sqrt(np.mean(err[use].reshape(len(tx),-1)**2,axis=0)))
            de=noise*np.sqrt(np.sum(inv[-1]**2))
            events.append(dict(cycle=cycle,time=center,difference=coef[-1].reshape(reference.shape),error=de.reshape(reference.shape),residual=res,projection=inv[-1],n_in=int(np.sum(shape>.05)),n_out=int(np.sum(shape==0))))
        if not events:raise RuntimeError(f'No covered events in sector {sector}; cannot localize source')
        diff,de=combine(events);roi=(xx-host[0])**2+(yy-host[1])**2<4.7**2
        jy,jx=np.unravel_index(np.nanargmax(diff/de),diff.shape)
        gf,_=gaussfit(diff,de,[jx,jy],roi)
        gf['sky_ra_dec']=w.all_pix2world([[gf['x'],gf['y']]],0)[0];gf['target_offset_pixels']=np.hypot(gf['x']-host[0],gf['y']-host[1])
        modes[mode]=dict(n_events=len(events),source_fit=gf,events=[{k:e[k] for k in ['cycle','time','n_in','n_out']} for e in events]);sets[mode]=events
        np.savez_compressed(out/f'difference_{mode}.npz',difference=diff,error=de,per_event_difference=np.array([e['difference'] for e in events]),per_event_error=np.array([e['error'] for e in events]))
    rng=np.random.default_rng(tic);events=sets['quadratic'];gf=modes['quadratic']['source_fit'];roi=(xx-host[0])**2+(yy-host[1])**2<4.7**2;boots=[]
    for _ in range(bootstrap):
        draws=[]
        for j in rng.integers(0,len(events),len(events)):
            e=events[j];res=e['residual'];n=len(res);starts=rng.integers(0,n,(n+4)//5);ix=np.concatenate([(s+np.arange(5))%n for s in starts])[:n]
            draws.append(dict(difference=e['difference']+(e['projection']@res[ix]).reshape(reference.shape),error=e['error']))
        dd,ee=combine(draws);bf,_=gaussfit(dd,ee,[gf['x'],gf['y']],roi);boots.append([bf['x'],bf['y']])
    boots=np.array(boots);stars=pd.read_csv(ROOT/f'data/broad_catalogs/tic{tic}_gaia_nearby.csv',dtype={'source_id':str})
    stars[['x','y']]=w.all_world2pix(stars[['ra','dec']].to_numpy(),0);stars['distance_difference_pixels']=np.hypot(stars.x-gf['x'],stars.y-gf['y']);stars['distance_target_pixels']=np.hypot(stars.x-host[0],stars.y-host[1]);stars.to_csv(out/'gaia_pixel_positions.csv',index=False)
    modes['quadratic']['bootstrap']=dict(n=bootstrap,median_xy=np.median(boots,axis=0),percentile16_84_xy=np.percentile(boots,[16,84],axis=0),target_offset_16_50_84=np.percentile(np.linalg.norm(boots-host,axis=1),[16,50,84]))
    dd=np.load(out/'difference_quadratic.npz');fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for ax,img,title in zip(axes,[np.log10(np.maximum(reference,1)),dd['difference'],dd['difference']/dd['error']],['Reference log10 e-/s','Eclipse deficit e-/s','Deficit / diagonal error']):
        im=ax.imshow(img,origin='lower');ax.plot(host[0],host[1],'rx',ms=10);ax.plot(gf['x'],gf['y'],'c+',ms=12)
        for _,s in stars.iterrows():
            if 0<=s.x<11 and 0<=s.y<11 and s.phot_g_mean_mag<18:ax.plot(s.x,s.y,'wo',mfc='none',ms=6)
        ax.set_title(title);fig.colorbar(im,ax=ax,fraction=.05)
    fig.savefig(out/'source_localization.png',dpi=170);plt.close(fig)
    np.savez_compressed(out/'reference_bootstrap.npz',reference=reference,target_xy=host,centroids=boots,time_btjd=t)
    (out/'ffi_product_ids.txt').write_text('\n'.join(sorted(set(ffi)))+'\n')
    manifest=dict(path=str(path),sha256=digest,bytes=path.stat().st_size,url=url,checked_utc=datetime.now(timezone.utc).isoformat(),checksums=checks,sector=sector,camera=header.get('CAMERA'),ccd=header.get('CCD'),processing_header=header,table_header=lh)
    old.write_text(json.dumps(safe(manifest),indent=2)+'\n')
    result=dict(tic=tic,sector=sector,event_phase=event_phase,target_xy=host,science_column=column,science_row=row,period_days=p,epoch_btjd=epoch,duration_days=duration,modes=modes,nearest_gaia_to_difference=stars.sort_values('distance_difference_pixels')[['source_id','phot_g_mean_mag','distance_difference_pixels','distance_target_pixels','x','y']].head(8).to_dict('records'),raw_manifest='raw_manifest.json',limitations=['Gaussian centroid and event/block bootstrap exclude PRF, WCS, aperture, and proper-motion systematics.','G17 neighbor only 0.2 pixel from giant cannot be conclusively separated without calibrated scene modeling.','Period fixed using early two sectors; no independent-period search used for source image.'],code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'pixel_vetting.json').write_text(json.dumps(safe(result),indent=2)+'\n');print(json.dumps(safe(dict(tic=tic,sector=sector,event_phase=event_phase,modes={k:dict(n_events=v['n_events'],source_fit=v['source_fit']) for k,v in modes.items()})),indent=2))


def detailed_photometry(tic=260352037):
    from broad_vet import fit_shape, fit_event_timings
    from pixel_vet_49617264 import safe
    curves,_=load(tic)
    eph=json.loads((OUT/f'tic{tic}_coherence.json').read_text())
    p,epoch=eph['period_days'],eph['epoch_btjd']
    local=[]
    for c in curves:
        stats,data=local_quadratic(c['t'],c['f'],c['e'],p,epoch,.065)
        data['exposure']=np.full(len(data['t']),c['dt']);local.append(data)
    joined={k:np.concatenate([v[k] for v in local]) for k in local[0]}
    shape=fit_shape(joined,p,epoch,.065,joined['exposure'])
    timings=fit_event_timings(local,[c['sector'] for c in curves],p,epoch,shape)
    fit=timings['fitted_ephemeris'];p=fit['period'];epoch=fit['epoch_btjd']
    rows=[]
    for c in curves:
        primary,data=local_quadratic(c['t'],c['f'],c['e'],p,epoch,shape['duration_days'])
        secondary,_=local_quadratic(c['t'],c['f'],c['e'],p,epoch+p/2,shape['duration_days'])
        rows.append(dict(sector=c['sector'],primary=primary,secondary_half_phase=secondary))
    pools={}
    for label in ['primary','secondary_half_phase']:
        depths=np.array([r[label]['depth'] for r in rows]);errors=np.array([r[label]['depth_err'] for r in rows]);weight=1/errors**2
        mean=np.average(depths,weights=weight);chi=np.sum(((depths-mean)/errors)**2);formal=1/np.sqrt(weight.sum())
        err=max(formal*np.sqrt(max(1,chi/(len(depths)-1))),np.std(depths,ddof=1)/np.sqrt(len(depths)))
        pools[label]=dict(depth_ppm=mean*1e6,error_ppm=err*1e6,conditional_snr=mean/err,sector_chi2=chi,n_sectors=len(depths),n_events=sum(r[label]['n_events'] for r in rows))
    # Sector-binned fold for the refined full-phase shape, with both event
    # families measured at pre-specified phases rather than selected peaks.
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for c in curves:
        x,b=bins(phase(c['t'],p,epoch)/p,c['flat'],160)
        axes[0].plot(x,(b-1)*1e6,color='.6',lw=.5,alpha=.3)
    t=np.concatenate([c['t'] for c in curves]);f=np.concatenate([c['flat'] for c in curves]);x,b=bins(phase(t,p,epoch)/p,f,200)
    axes[0].plot(x,(b-1)*1e6,'k.-',ms=3);axes[0].set(xlabel='Phase at refined 0.88-day period',ylabel='Normalized raw SAP − 1 (ppm)')
    sectors=[r['sector'] for r in rows]
    for label,color in [('primary','tab:blue'),('secondary_half_phase','tab:orange')]:
        axes[1].errorbar(sectors,[r[label]['depth']*1e6 for r in rows],yerr=[r[label]['depth_err']*1e6 for r in rows],fmt='.',color=color,label=label)
    axes[1].legend();axes[1].set(xlabel='TESS sector',ylabel='Local quadratic window mean depth (ppm)')
    fig.savefig(OUT/f'tic{tic}_eclipse_photometry.png',dpi=160);plt.close(fig)
    result=dict(tic=tic,shape=shape,timings=timings,refined_ephemeris=fit,sectors=rows,sector_pooled_depths=pools,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitations=['Conditional statistical errors; correlated giant variability and template errors remain.','Local mean depths differ from fitted central depths because eclipse ingress contributes.','Host giant parameters cannot be assigned to the unresolved eclipse source.'])
    (OUT/f'tic{tic}_photometry.json').write_text(json.dumps(safe(result),indent=2)+'\n')
    print(json.dumps(safe(dict(ephemeris=fit,shape=shape,pools=pools)),indent=2),flush=True)


def baseline_controls(tic=260352037):
    """Compare raw-polynomial baselines and non-eclipse orbital phases.

    Orbital control phases share data and correlated giant noise. Their range
    is a useful falsification check, not an independent null distribution or
    a survey-wide false-alarm probability.
    """
    from pixel_vet_49617264 import safe
    curves,_=load(tic)
    phot=json.loads((OUT/f'tic{tic}_photometry.json').read_text())
    p=phot['refined_ephemeris']['period'];epoch=phot['refined_ephemeris']['epoch_btjd'];duration=phot['shape']['duration_days']
    def pool(rows):
        values=[r for r in rows if r['depth'] is not None and r['n_events']>=2]
        d=np.array([r['depth'] for r in values]);e=np.array([r['depth_err'] for r in values]);w=1/e**2
        mean=np.average(d,weights=w);chi=np.sum(((d-mean)/e)**2)
        err=max(np.sqrt(1/w.sum())*np.sqrt(max(1,chi/(len(d)-1))),np.std(d,ddof=1)/np.sqrt(len(d)))
        return dict(depth_ppm=mean*1e6,error_ppm=err*1e6,n_sectors=len(d),n_events=sum(r['n_events'] for r in values))
    alternatives=[]
    for degree in [1,2,3]:
        row=dict(polynomial_degree=degree)
        for ph,label in [(0.,'primary'),(.5,'secondary')]:
            measures=[local_quadratic(c['t'],c['f'],c['e'],p,epoch+ph*p,duration,degree=degree)[0] for c in curves]
            row[label]=pool(measures)
        alternatives.append(row)
    phases=np.r_[np.linspace(.12,.38,12),np.linspace(.62,.88,12)];controls=[]
    for ph in phases:
        values=[local_quadratic(c['t'],c['f'],c['e'],p,epoch+ph*p,duration)[0] for c in curves]
        controls.append(dict(phase=float(ph),**pool(values)))
    result=dict(tic=tic,alternative_polynomial_baselines=alternatives,non_eclipse_phase_controls=controls,limitations=['Phase controls overlap in time and share giant variability; they are not independent samples.','Some control-phase out-of-event baselines include the real primary or secondary; this can bias their apparent depths negative.','These targeted controls have no survey-wide false-alarm interpretation.'],code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/f'tic{tic}_baseline_controls.json').write_text(json.dumps(safe(result),indent=2)+'\n')
    print(json.dumps(safe(alternatives),indent=2));print('Control depth rangeppm',min(r['depth_ppm'] for r in controls),max(r['depth_ppm'] for r in controls))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['coherence','photometry','pixels','controls'],default='coherence')
    parser.add_argument('--sectors',nargs='+',type=int,default=[28])
    parser.add_argument('--bootstrap',type=int,default=200)
    parser.add_argument('--event-phase',type=float,choices=[0.,.5],default=0.)
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.mode=='coherence':
        for tic in CONFIG:analyze(tic)
    elif args.mode=='photometry':detailed_photometry()
    elif args.mode=='controls':baseline_controls()
    else:
        for sector in args.sectors:pixel_audit(sector=sector,bootstrap=args.bootstrap,event_phase=args.event_phase)
