"""Independent falsification audit of the 0.8808-day blend near TIC 260352037.

Retains all observations and uses the published old ephemeris as a holdout
prediction before any refit. Raw local polynomial fits and empirical image
registration supplement, rather than replace, the historical results.
"""
from pathlib import Path
import argparse, hashlib, json, os, sys
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from scipy.optimize import minimize_scalar, least_squares
from scipy.ndimage import shift, median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import verify_fits_checksums
from broad_vet import integrated_trapezoid, phase
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'candidates/stellar_systems/tic260352037_validation'
DATA=ROOT/'data/stellar_systems/tic260352037_validation'
CACHE=Path('/workspace/tess-research-data/stellar_systems/tic260352037_validation')
P0=.880788466
E0=1326.593218
D0=.07718209037083432
ING0=D0*.4335096557700692
TIC=260352037

def safe(x):
    if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [safe(v) for v in x]
    if isinstance(x,np.ndarray):return safe(x.tolist())
    if isinstance(x,(np.floating,float)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.integer):return int(x)
    return x

def save(name,obj,folder=OUT):
    folder.mkdir(parents=True,exist_ok=True)
    (folder/name).write_text(json.dumps(safe(obj),indent=2,allow_nan=False)+'\n')

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def retrieve():
    DATA.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
    url='https://mast.stsci.edu/tesscut/api/v0.1/sector?ra=94.787911339&dec=-57.078220725'
    r=requests.get(url,timeout=60);r.raise_for_status();(DATA/'tesscut_coverage.json').write_bytes(r.content)
    sectors=[int(x['sector']) for x in r.json()['results'] if int(x['sector'])>39]
    # Do not use known affected QLP v01 timing. The present sky position has
    # no coverage in 74-79 or 99-104, but explicit version/skip guards persist.
    def one(s):
        version='v02' if 74<=s<=79 else 'v01'
        if 99<=s<=104:return dict(sector=s,status='skipped_known_v01_time_warning')
        tic=str(TIC).zfill(16);sub='/'.join(tic[i:i+4] for i in range(0,16,4))
        name=f'hlsp_qlp_tess_ffi_s{s:04}-{tic}_tess_{version}_llc.fits'
        url=f'https://archive.stsci.edu/hlsps/qlp/s{s:04}/{sub}/{name}'
        path=CACHE/name
        if not path.exists():
            rr=requests.get(url,timeout=60)
            if rr.status_code==404:return dict(sector=s,status='not_found',url=url)
            rr.raise_for_status();temp=path.with_suffix('.part');temp.write_bytes(rr.content);temp.replace(path)
        with fits.open(path) as h:
            h.verify('exception');checks=verify_fits_checksums(path,h)
            meta={k:h[1].header.get(k,h[0].header.get(k)) for k in ['TIMESYS','BJDREFI','BJDREFF','TIMEDEL','TESSMAG','RADIUS','TEFF','CAMERA','CCD','PROCVER']}
            if meta['TIMESYS']!='TDB' or (meta['BJDREFI'] or 0)+(meta['BJDREFF'] or 0)!=2457000:raise ValueError('Time system not BTJD TDB')
        rec=dict(tic=TIC,sector=s,status='retained',version=version,url=url,path=str(path),sha256=digest(path),bytes=path.stat().st_size,checks=checks,metadata=meta,retrieved_utc=datetime.now(timezone.utc).isoformat())
        print('retained later sector',s,flush=True);return rec
    with ThreadPoolExecutor(4) as pool:records=list(pool.map(one,sectors))
    old=DATA/'later_observations.json'
    if old.exists():
        prior={r['sector']:r for r in json.loads(old.read_text())}
        for r in records:
            if r['status']=='retained' and prior.get(r['sector'],{}).get('sha256',r['sha256'])!=r['sha256']:raise ValueError('Retained raw bytes changed')
    save('later_observations.json',records,DATA)
    return records

def load():
    qa=json.loads((ROOT/'candidates/broad_s1/giant_survivors/tic260352037_qa.json').read_text())
    records=[dict(r['source'],sector=r['sector'],era='old') for r in qa['independent_sector_qa']]
    if (DATA/'later_observations.json').exists():records += [dict(r,era='later') for r in json.loads((DATA/'later_observations.json').read_text()) if r['status']=='retained']
    curves=[]
    for r in records:
        path=Path(r['path'])
        if digest(path)!=r['sha256']:raise ValueError('Immutable hash mismatch')
        with fits.open(path) as h:
            verify_fits_checksums(path,h);d=h[1].data
            errcol='KSPSAP_FLUX_ERR' if 'KSPSAP_FLUX_ERR' in d.names else 'DET_FLUX_ERR'
            cols={key:np.asarray(d[key],float) for key in ['TIME','QUALITY','SAP_FLUX']}
            cols['KSPSAP_FLUX_ERR']=np.asarray(d[errcol],float)
            dt=float(h[1].header.get('TIMEDEL',np.nan))
        good=(cols['QUALITY']==0)&np.isfinite(cols['TIME'])&np.isfinite(cols['SAP_FLUX'])&(cols['SAP_FLUX']>0)&np.isfinite(cols['KSPSAP_FLUX_ERR'])&(cols['KSPSAP_FLUX_ERR']>0)
        t=cols['TIME'][good];f=cols['SAP_FLUX'][good];e=cols['KSPSAP_FLUX_ERR'][good];norm=np.median(f)
        if not np.isfinite(dt):dt=float(np.median(np.diff(t)))
        curves.append(dict(t=t,f=f/norm,e=e/norm,dt=dt,sector=r['sector'],era=r['era'],source=r))
    return curves

def eventfits(c,p=P0,epoch=E0,event_phase=0.,degree=2,duration=D0,ingress=ING0,window=.24,mask_actual=False):
    t,f,e,dt=c['t'],c['f'],c['e'],c['dt'];epoch+=event_phase*p
    rows=[]
    for cy in range(int(np.ceil((t.min()-epoch)/p)),int(np.floor((t.max()-epoch)/p))+1):
        center=epoch+cy*p;use=abs(t-center)<window
        if mask_actual:use &= (abs(phase(t,P0,E0))>D0*.7)&(abs(phase(t,P0,E0+P0/2))>D0*.7)
        x=t[use]-center;y=f[use];er=e[use]
        shape=integrated_trapezoid(x,duration,ingress,dt)
        if len(x)<degree+7 or np.sum(shape>.1)<2 or np.sum((shape==0)&(x<0))<3 or np.sum((shape==0)&(x>0))<3:continue
        if not mask_actual and np.max(np.diff(x))>max(.04,3.01*dt):continue
        design=np.column_stack([x**k for k in range(degree+1)]+[-shape])
        weight=1/np.maximum(er,np.median(er))**2
        inv=np.linalg.pinv(design*np.sqrt(weight[:,None]))*np.sqrt(weight[None,:])
        beta=inv@y;res=y-design@beta
        floor=max(np.std(res,ddof=design.shape[1]),np.median(er))
        diagonal=floor*np.sqrt(np.sum(inv[-1]**2))
        # Analytic Bartlett/HAC-style residual covariance, estimating a
        # five-cadence correlation length. Data are never permuted in this fit.
        projection=inv[-1];var=diagonal**2
        for lag in range(1,min(6,len(res)//3)):
            cov=np.mean(res[:-lag]*res[lag:]);var+=2*(1-lag/6)*cov*np.sum(projection[:-lag]*projection[lag:])
        hac=np.sqrt(max(diagonal**2,var))
        rows.append(dict(sector=c['sector'],era=c['era'],cycle=cy,time=center,depth=beta[-1],error=max(diagonal,hac),diagonal_error=diagonal,n_in=int(np.sum(shape>.1)),n_points=len(x)))
    return rows

def pooled(rows):
    if not rows:return None
    table=pd.DataFrame(rows);a=[]
    for sec,d in table.groupby('sector'):
        w=1/d.error.to_numpy()**2;dep=np.average(d.depth,weights=w)
        err=max(1/np.sqrt(w.sum()),np.std(d.depth,ddof=1)/np.sqrt(len(d)))
        a.append(dict(sector=int(sec),depth=dep,error=err,n=len(d)))
    v=pd.DataFrame(a);w=1/v.error.to_numpy()**2;mean=np.average(v.depth,weights=w);chi=np.sum(((v.depth-mean)/v.error)**2)
    err=max(1/np.sqrt(w.sum())*np.sqrt(max(1,chi/max(1,len(v)-1))),np.std(v.depth,ddof=1)/np.sqrt(len(v)) if len(v)>1 else 0)
    return dict(depth_ppm=mean*1e6,error_ppm=err*1e6,conditional_ratio=mean/err,n_events=len(rows),n_sectors=len(v),sector_chi2=chi,sectors=a)

def profiles(curves,p=P0,epoch=E0):
    """Sector-constant raw-local eclipse and harmonic morphology comparisons."""
    tt=[];ff=[];ee=[];ss=[];dt=[]
    for c in curves:
        # Raw quadratic baseline fit independently at *all* phases is avoided:
        # a two-day median is used only for broad display/model comparison.
        n=max(5,int(2/c['dt'])//2*2+1);trend=median_filter(c['f'],size=n,mode='nearest')
        tt.extend(c['t']);ff.extend(c['f']/trend-1);ee.extend(c['e']);ss.extend([c['sector']]*len(c['t']));dt.extend([c['dt']]*len(c['t']))
    t,y,e,sectors,dt=map(np.asarray,[tt,ff,ee,ss,dt]);pha=phase(t,p,epoch)
    sec=[(sectors==s).astype(float) for s in np.unique(sectors)]
    s1=integrated_trapezoid(pha,D0,ING0,dt);s2=integrated_trapezoid(phase(t,p,epoch+p/2),D0,ING0,dt)
    models={};indices={}
    # Sector offsets and sector-specific eclipse depths absorb different apertures.
    x0=np.column_stack(sec);shapes=[np.column_stack([-s1*v,-s2*v]) for v in sec]
    x=np.column_stack([x0]+shapes);b=np.linalg.lstsq(x,y,rcond=None)[0];res=y-x@b
    for label,cols in [('constant',[]),('sinusoid_one',[np.cos(2*np.pi*pha/p),np.sin(2*np.pi*pha/p)]),('harmonic_two',[v for k in range(1,3) for v in [np.cos(2*np.pi*k*pha/p),np.sin(2*np.pi*k*pha/p)]]),('harmonic_eight',[v for k in range(1,9) for v in [np.cos(2*np.pi*k*pha/p),np.sin(2*np.pi*k*pha/p)]])]:
        design=np.column_stack([x0]+cols);coef=np.linalg.lstsq(design,y,rcond=None)[0];rr=y-design@coef;models[label]=dict(rss=np.sum(rr**2),n_parameters=design.shape[1])
    models['two_eclipse_templates']=dict(rss=np.sum(res**2),n_parameters=x.shape[1]);models['limitations']='RSS is descriptive: median normalization, correlated noise, unequal model parameter counts; no formal model odds.'
    return (t,y,sectors,dt),models

def photometry():
    OUT.mkdir(parents=True,exist_ok=True);curves=load();later=[c for c in curves if c['era']=='later'];old=[c for c in curves if c['era']=='old']
    rows=[];summary=[]
    for degree in [1,2,3]:
        primary=sum([eventfits(c,degree=degree) for c in curves],[]);secondary=sum([eventfits(c,event_phase=.5,degree=degree) for c in curves],[])
        for group,eras in [('old',['old']),('later',['later']),('all',['old','later'])]:
            summary.append(dict(baseline_degree=degree,group=group,primary=pooled([v for v in primary if v['era'] in eras]),secondary=pooled([v for v in secondary if v['era'] in eras])))
        if degree==2:rows=[dict(r,family='primary') for r in primary]+[dict(r,family='secondary') for r in secondary]
    table=pd.DataFrame(rows);table.to_csv(OUT/'event_depths.csv',index=False)
    # P/2 forces all eclipse families to identical shape/depth, which can be
    # falsified by persistent alternation. 2P adds another primary/secondary
    # parity; absent parity alone cannot prove the actual orbital period.
    parity=[]
    for family in ['primary','secondary']:
        ff=table[table.family==family]
        for era in ['old','later','all']:
            d=ff if era=='all' else ff[ff.era==era]
            # Per-sector contrasts retain arbitrary sector depths. Bootstrap
            # whole sectors provides a conservative correlated-noise error.
            contrasts=[]
            for sec,r in d.groupby('sector'):
                a=[r[r.cycle%2==par] for par in [0,1]]
                if min(map(len,a))<3:continue
                means=[np.average(v.depth,weights=1/v.error**2) for v in a];errs=[max(np.sqrt(1/np.sum(1/v.error**2)),np.std(v.depth,ddof=1)/np.sqrt(len(v))) for v in a]
                contrasts.append(dict(sector=int(sec),contrast=means[0]-means[1],error=np.hypot(*errs)))
            if not contrasts:continue
            z=pd.DataFrame(contrasts);w=1/z.error.to_numpy()**2;mean=np.average(z.contrast,weights=w);chi=np.sum(((z.contrast-mean)/z.error)**2);err=max(np.sqrt(1/w.sum())*np.sqrt(max(1,chi/max(1,len(z)-1))),np.std(z.contrast,ddof=1)/np.sqrt(len(z)))
            parity.append(dict(family=family,era=era,contrast_ppm=mean*1e6,error_ppm=err*1e6,conditional_ratio=mean/err,n_sectors=len(z),sector_contrasts=contrasts))
    # Observed phase controls are sensitivity checks and are explicitly not
    # an independent null distribution or global false-alarm probability.
    controls=[]
    for ph in np.arange(.1,.951,.05):
        if abs(ph-.5)<.09:continue
        rr=sum([eventfits(c,event_phase=float(ph)) for c in later],[])
        controls.append(dict(phase=ph,summary=pooled(rr)))
    masked_controls=[]
    for ph in [.15,.2,.25,.3,.35,.65,.7,.75,.8,.85]:
        rr=sum([eventfits(c,event_phase=float(ph),mask_actual=True) for c in later],[])
        masked_controls.append(dict(phase=ph,summary=pooled(rr)))
    raw,models=profiles(curves)
    # Derive per-sector eclipse timing shifts from *all* raw events, avoiding
    # selection on individual-event S/N. Sector profile scans use raw local
    # polynomial regression and actual exposure integration.
    sector_time=[]
    for c in curves:
        values=[]
        for fam,ph in [('primary',0.),('secondary',.5)]:
            # Data windows fixed to old centers, maintaining no phase search
            # outside ±.04 day; discovery ephemeris remains unmodified above.
            blocks=[]
            for cy in range(int(np.ceil((c['t'].min()-E0-ph*P0)/P0)),int(np.floor((c['t'].max()-E0-ph*P0)/P0))+1):
                center=E0+(cy+ph)*P0;use=abs(c['t']-center)<.24;x=c['t'][use]-center;y=c['f'][use]
                if len(x)<10 or np.sum(abs(x)<D0/2)<2 or np.max(np.diff(x))>max(.04,3.01*c['dt']):continue
                base=np.column_stack([np.ones(len(x)),x,x*x]);proj=np.eye(len(x))-base@np.linalg.pinv(base);blocks.append((x,proj@y,proj,cy))
            if len(blocks)<4:continue
            shifts=np.linspace(-.045,.045,181);rss=[]
            for ds in shifts:
                total=0
                for x,y,proj,cy in blocks:
                    shape=proj@integrated_trapezoid(x-ds,D0,ING0,c['dt']);dep=max(0,np.dot(shape,-y)/np.dot(shape,shape));total+=np.sum((y+dep*shape)**2)
                rss.append(total)
            best=int(np.argmin(rss));noise=np.median(np.concatenate([y for x,y,proj,cy in blocks])**2)*1.5
            # Correlated-noise floor from leave-one-event timing profiles would
            # be very costly; 3-minute floor prevents pixel-level precision.
            sigma=0.5*(shifts[np.array(rss)<=min(rss)+noise][-1]-shifts[np.array(rss)<=min(rss)+noise][0]);sigma=max(sigma,3/1440)
            values.append(dict(family=fam,sector=c['sector'],era=c['era'],mean_cycle=np.mean([cy for x,y,proj,cy in blocks])+ph,shift_days=shifts[best],error_days=sigma,n_events=len(blocks),bound=best in [0,len(shifts)-1]))
        sector_time+=values
    pd.DataFrame(sector_time).to_csv(OUT/'sector_timings.csv',index=False)
    primary_times=[v for v in sector_time if v['family']=='primary' and not v['bound']];dd=pd.DataFrame(primary_times);cy=dd.mean_cycle.to_numpy();tt=dd.shift_days.to_numpy();er=dd.error_days.to_numpy();ref=np.mean(cy);design=np.column_stack([np.ones(len(cy)),cy-ref]);a=design/er[:,None];b=tt/er;coef=np.linalg.lstsq(a,b,rcond=None)[0];res=tt-design@coef;chi=np.sum((res/er)**2);cov=np.linalg.pinv(a.T@a)*max(1,chi/(len(tt)-2));period=P0+coef[1];epoch=E0+coef[0]-coef[1]*ref
    # Cluster bootstrap entire sectors, with no exposure-level shuffle.
    rng=np.random.default_rng(260352037);draw=[]
    for _ in range(1000):
        ids=rng.integers(0,len(cy),len(cy));bb=np.linalg.lstsq(a[ids],b[ids],rcond=None)[0];draw.append(bb)
    draw=np.array(draw);perr=max(np.sqrt(cov[1,1]),np.std(draw[:,1],ddof=1));eperr=max(np.sqrt(np.array([1,-ref])@cov@np.array([1,-ref])),np.std(draw[:,0]-draw[:,1]*ref,ddof=1))
    eph=dict(period_days=period,period_error_days=perr,epoch_btjd=epoch,epoch_error_days=eperr,n_sector_timings=len(tt),conditional_chi2=chi,dof=len(tt)-2,old_fixed_period_days=P0,old_fixed_epoch_btjd=E0,sector_bootstrap_n=1000,limitations='Fixed trapezoid, giant variability and conditional correlation floor. This is not an orbital dynamics model.')
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained');t,y,secs,dt=raw
    for era,color in [('old','tab:blue'),('later','tab:orange')]:
        selected=[c['sector'] for c in curves if c['era']==era];use=np.isin(secs,selected);ph=phase(t[use],P0,E0)/P0;edges=np.linspace(-.5,.5,151);inds=np.digitize(ph,edges)-1;med=np.array([np.median(y[use][inds==i])*1e6 if np.sum(inds==i)>0 else np.nan for i in range(150)])
        axes[0,0].plot((edges[:-1]+edges[1:])/2,med,'.-',label=era,color=color,ms=3)
    axes[0,0].set(xlabel='Phase: old ephemeris held fixed',ylabel='Median-normalized raw SAP (ppm)');axes[0,0].legend()
    for family,color in [('primary','tab:blue'),('secondary','tab:orange')]:
        d=table[table.family==family];pool=pooled(d.to_dict('records'));se=pool['sectors'];axes[0,1].errorbar([r['sector'] for r in se],[r['depth']*1e6 for r in se],yerr=[r['error']*1e6 for r in se],fmt='.',label=family,color=color)
        tim=[v for v in sector_time if v['family']==family];axes[1,0].errorbar([v['mean_cycle'] for v in tim],[v['shift_days']*1440 for v in tim],yerr=[v['error_days']*1440 for v in tim],fmt='.',label=family,color=color)
    axes[0,1].set(xlabel='Sector',ylabel='Raw simultaneous template depth (ppm)');axes[0,1].legend()
    axes[1,0].set(xlabel='Eclipse cycle since old epoch',ylabel='O − C against old prediction (minutes)');axes[1,0].legend()
    con=[r for r in masked_controls if r['summary']];axes[1,1].errorbar([r['phase'] for r in con],[r['summary']['depth_ppm'] for r in con],yerr=[r['summary']['error_ppm'] for r in con],fmt='.');axes[1,1].axhline(0,color='k',lw=.5);axes[1,1].set(xlabel='Later control phase: actual eclipses masked',ylabel='Apparent template depth (ppm)')
    fig.savefig(OUT/'holdout_diagnostics.png',dpi=170);plt.close(fig)
    save('photometry.json',dict(tic=TIC,checked_utc=datetime.now(timezone.utc).isoformat(),n_sectors=len(curves),sectors=[c['sector'] for c in curves],later_sectors=[c['sector'] for c in later],frozen_old_ephemeris=dict(period_days=P0,epoch_btjd=E0),template=dict(duration_days=D0,ingress_days=ING0),depth_summaries=summary,twice_period_parity=parity,later_phase_controls=controls,masked_actual_eclipses_later_phase_controls=masked_controls,morphology_descriptive_models=models,sector_timing_refit=eph,code_sha256=digest(__file__),limitations=['All depth ratios are targeted and conditional, not calibrated global false-alarm probabilities.','QLP SAP errors use KSPSAP reported errors supplemented by event residual/HAC and event/sector scatter floors.','Native exposure smearing integrated; empirical five-cadence covariance excludes arbitrary long-time noise.','Unmasked quarter-phase controls include real eclipse tails in their baseline windows and manufacture template-like deficits; corrected controls mask both eclipse families.','2P with equal orbital cycles remains an unavoidable photometric alias in the absence of RV.']))
    print(json.dumps(safe(dict(summary=summary,eph=eph,parity=parity)),indent=2),flush=True)



def morphology():
    from broad_vet import local_quadratic,fit_shape
    curves=[c for c in load() if c['era']=='later'];families={}
    for family,ph in [('primary',0.),('secondary',.5)]:
        locals_=[];individual=[]
        for c in curves:
            stats,data=local_quadratic(c['t'],c['f'],c['e'],P0,E0+ph*P0,.10)
            if data is None:continue
            data['exposure']=np.full(len(data['t']),c['dt']);locals_.append(data)
            fit=fit_shape(data,P0,E0+ph*P0,D0,data['exposure']);individual.append(dict(sector=c['sector'],shape=fit))
        joined={key:np.concatenate([d[key] for d in locals_]) for key in locals_[0]}
        fit=fit_shape(joined,P0,E0+ph*P0,D0,joined['exposure'])
        fields=[('duration_days','duration_err'),('epoch_btjd','epoch_err'),('depth','depth_err')]
        scatter={}
        for field,errkey in fields:
            values=np.array([r['shape'][field] for r in individual]);errs=np.array([r['shape'][errkey] for r in individual]);weights=1/errs**2;mean=np.average(values,weights=weights);chi=np.sum(((values-mean)/errs)**2);err=max(np.sqrt(1/weights.sum())*np.sqrt(max(1,chi/(len(values)-1))),np.std(values,ddof=1)/np.sqrt(len(values)))
            scatter[field]=dict(weighted_mean=mean,sector_heterogeneity_error=err,sector_chi2=chi)
        families[family]=dict(joint_shape_fit=fit,individual_sector_shapes=individual,sector_parameter_summary=scatter)
    table=pd.read_csv(OUT/'event_depths.csv');table=table[table.era=='later'].copy();table['block']=table.sector.astype(str)+'_'+(table.cycle//4).astype(str)
    # Four-cycle synchronized sign reversals preserve magnitudes and the
    # primary/secondary red-noise association in each block; this is a
    # conditional symmetry diagnostic, not a discovery-wide randomization.
    keys=sorted(table.block.unique());rng=np.random.default_rng(TIC);draws=[]
    for _ in range(512):
        sg=dict(zip(keys,rng.choice([-1,1],len(keys))));r=table.copy();r['depth']=r.depth*r.block.map(sg)
        vals={}
        for fam in ['primary','secondary']:
            ff=r[r.family==fam];w=1/ff.error.to_numpy()**2;vals[fam]=np.average(ff.depth,weights=w)*1e6
        draws.append(vals)
    null={}
    for fam in ['primary','secondary']:
        ff=table[table.family==fam];w=1/ff.error.to_numpy()**2;observed=np.average(ff.depth,weights=w)*1e6;vals=np.array([d[fam] for d in draws]);null[fam]=dict(observed_inverse_variance_mean_ppm=observed,null_percentile1_50_99_ppm=np.percentile(vals,[1,50,99]),max_abs_null_ppm=np.max(abs(vals)),n_null_at_least_observed=int(np.sum(abs(vals)>=abs(observed))))
    save('morphology.json',dict(tic=TIC,families=families,sign_block_control=dict(n=512,cycles_per_block=4,n_blocks=len(keys),result=null),code_sha256=digest(__file__),limitations=['Shape covariances are conditional diagonal approximations. Sector heterogeneity summaries do not include all template/model uncertainties.','The later-sector trapezoids are phenomenological and do not determine component radii or temperatures.','Sign-reversal validity requires symmetric block errors; selecting on the known signal and few cadence aliases excludes global false-alarm interpretation.']))
    print(json.dumps(safe(dict(families={k:dict(joint=v['joint_shape_fit'],sector=v['sector_parameter_summary']) for k,v in families.items()},sign_block=null)),indent=2),flush=True)


def timing_stability():
    t=pd.read_csv(OUT/'sector_timings.csv');t=t[(t.family=='primary')&~t.bound];rows=[]
    for era in ['old','later','all']:
        d=t if era=='all' else t[t.era==era];x=d.mean_cycle.to_numpy();y=d.shift_days.to_numpy();e=d.error_days.to_numpy();ref=np.mean(x);A=np.column_stack([np.ones(len(x)),x-ref]);a=A/e[:,None];b=y/e;beta=np.linalg.lstsq(a,b,rcond=None)[0];chi=np.sum(((y-A@beta)/e)**2);cov=np.linalg.pinv(a.T@a)*max(1,chi/(len(y)-2))
        rows.append(dict(era=era,period_days=P0+beta[1],period_error_days=np.sqrt(cov[1,1]),chi2=chi,dof=len(y)-2))
    x=t.mean_cycle.to_numpy();y=t.shift_days.to_numpy();e=t.error_days.to_numpy();ref=np.mean(x);A=np.column_stack([np.ones(len(x)),x-ref,(x-ref)**2]);a=A/e[:,None];b=y/e;beta=np.linalg.lstsq(a,b,rcond=None)[0];chi=np.sum(((y-A@beta)/e)**2);cov=np.linalg.pinv(a.T@a)*max(1,chi/(len(y)-3));pdot=2*beta[2]/P0;perr=2*np.sqrt(cov[2,2])/P0
    save('timing_stability.json',dict(era_fits=rows,quadratic_fit=dict(chi2=chi,dof=len(y)-3,quadratic_day_per_cycle2=beta[2],quadratic_error=np.sqrt(cov[2,2]),period_derivative_sec_per_year=pdot*86400*365.25,error_sec_per_year=perr*86400*365.25),code_sha256=digest(__file__),limitations='Conditional sector profiles and 3-minute correlation floors, no significance correction or binary dynamics inference.'))


def wide_pair_dynamics():
    from astropy import constants as const,units as u
    source=ROOT/'data/stellar_systems/hidden_sources/tic260352037_astrometry.csv'
    stars=pd.read_csv(source,dtype={'source_id':str});a=stars[stars.source_id=='5496090260794096128'].iloc[0];b=stars[stars.source_id=='5496090260793489024'].iloc[0]
    sep=np.hypot((a.ra-b.ra)*np.cos(np.deg2rad(a.dec)),a.dec-b.dec)*3600;mu=np.array([b.pmra-a.pmra,b.pmdec-a.pmdec]);cov=np.zeros((2,2))
    for v in [a,b]:cov+=np.array([[v.pmra_error**2,v.pmra_pmdec_corr*v.pmra_error*v.pmdec_error],[v.pmra_pmdec_corr*v.pmra_error*v.pmdec_error,v.pmdec_error**2]])
    dist=1000/a.parallax;vel=4.74047*np.linalg.norm(mu)*dist/1000;proj=sep*dist;mass=(proj*u.au*(vel*u.km/u.s)**2/(2*const.G)).to_value(u.M_sun)
    rng=np.random.default_rng(TIC);md=rng.multivariate_normal(mu,cov,100000);pd_=rng.normal(a.parallax,a.parallax_error,100000);d=1000/pd_;v=4.74047*np.linalg.norm(md,axis=1)*d/1000;mm=(sep*d*u.au*(v*u.km/u.s)**2/(2*const.G)).to_value(u.M_sun)
    save('wide_pair_dynamics.json',dict(separation_arcsec=sep,projected_separation_au=proj,relative_pm_mas_year=mu,relative_pm_covariance=cov,nominal_transverse_speed_km_s=vel,transverse_speed_16_50_84_km_s=np.percentile(v,[16,50,84]),minimum_binding_mass_nominal_msolar=mass,minimum_binding_mass_16_50_84_msolar=np.percentile(mm,[16,50,84]),input_sha256=digest(source),code_sha256=digest(__file__),limitations=['Only transverse motion and projected separation; actual separation and radial motion increase the binding requirement.','Monte Carlo includes within-source proper-motion covariance but assumes zero cross-source angular covariance and neglects parallax/PM covariance.','No probability of association is inferred; minimum mass is conditional on a common distance, which is not established. Nonlinear norm statistics are not a binding posterior.']))

def empirical_fit(diff,error,reference,roi,initial=(0.,0.),fixed=None):
    """Fit translated giant-dominated reference image, avoiding Gaussian PRF.

    Translating a sampled PRF is approximate: color and intra-pixel response
    differ between sources. The fit is a sensitivity test, not calibrated PRF
    astrometry or independent source confirmation.
    """
    template=reference/np.sum(reference)
    err=np.maximum(error,np.percentile(error[roi],10))
    def linear(dx,dy):
        ps=shift(template,(dy,dx),order=3,mode='constant',cval=0.,prefilter=True)
        x=np.column_stack([ps[roi],np.ones(np.sum(roi))]);w=1/err[roi]
        beta=np.linalg.lstsq(x*w[:,None],diff[roi]*w,rcond=None)[0]
        return (diff[roi]-x@beta)*w,beta,ps
    if fixed is None:
        fit=least_squares(lambda z:linear(*z)[0],initial,bounds=([-.75,-.75],[.75,.75]),max_nfev=70)
        xy=fit.x;res,beta,ps=linear(*xy);cov=np.linalg.pinv(fit.jac.T@fit.jac)*max(1,np.sum(res**2)/max(1,len(res)-4));errxy=np.sqrt(np.diag(cov))
    else:
        xy=np.array(fixed);res,beta,ps=linear(*xy);errxy=None
    return dict(dx=xy[0],dy=xy[1],shift_bound_hit=bool(np.any(abs(xy)>.749)),diagonal_error_xy=errxy,total_difference_flux=beta[0],background=beta[1],chi2_diagonal=np.sum(res**2),n_pixels=int(np.sum(roi))),ps*beta[0]+beta[1]


def new_pixel_cube(sector=95):
    import io,urllib.parse,zipfile
    folder=CACHE/f'pixels_s{sector:04}';folder.mkdir(parents=True,exist_ok=True)
    files=list(folder.glob('*.fits'))
    if not files:
        files=list(Path('/workspace/tess-research-data/stellar_systems/hidden_sources/tic260352037').glob(f'tess-s{sector:04}-*.fits'))
    url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urllib.parse.urlencode(dict(ra=94.787911339,dec=-57.078220725,x=11,y=11,units='px',sector=sector))
    if not files:
        r=requests.get(url,timeout=90);r.raise_for_status()
        if len(r.content)>100000000:raise ValueError('100 MB cube bound exceeded')
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            for name in z.namelist():
                if name.endswith('.fits'):(folder/Path(name).name).write_bytes(z.read(name))
        files=list(folder.glob('*.fits'))
    if len(files)!=1:raise ValueError('Expected one cube')
    path=files[0]
    from astropy.wcs import WCS
    with fits.open(path) as h:
        h.verify('exception');checks=verify_fits_checksums(path,h);w=WCS(h[2].header);head=dict(h[0].header);tab=dict(h[1].header);d=h[1].data
        t=np.array(d['TIME']);q=np.array(d['QUALITY']);flux=np.array(d['FLUX'],float);err=np.array(d['FLUX_ERR'],float);ffi=np.array(d['FFI_FILE']).astype(str)
    target=w.all_world2pix([[94.787911339,-57.078220725]],0)[0];col=tab['1CRV5P']+(target[0]+1-tab['1CRP5P'])*tab['1CDL5P'];row=tab['2CRV5P']+(target[1]+1-tab['2CRP5P'])*tab['2CDL5P']
    if not (44<=col<2092 and 0<=row<2048):raise ValueError('Target not science silicon')
    good=(q==0)&np.isfinite(t)&np.all(np.isfinite(flux)&np.isfinite(err),axis=(1,2));t=t[good];flux=flux[good];err=err[good]
    reference=np.median(flux,axis=0);blank=reference<np.percentile(reference,20);bg=np.median(flux[:,blank],axis=1);flux=flux-bg[:,None,None];reference=np.median(flux,axis=0);dt=tab['TIMEDEL']
    out=OUT/f'sector{sector}';out.mkdir(parents=True,exist_ok=True)
    manifest=dict(sector=sector,path=str(path),sha256=digest(path),bytes=path.stat().st_size,url=url,checks=checks,camera=head.get('CAMERA'),ccd=head.get('CCD'),science_column=col,science_row=row,time_table_header=tab,primary_header=head,checked_utc=datetime.now(timezone.utc).isoformat())
    prior=out/'raw_manifest.json'
    if prior.exists() and json.loads(prior.read_text())['sha256']!=manifest['sha256']:raise ValueError('Raw cube changed')
    save('raw_manifest.json',manifest,out)
    (out/'ffi_product_ids.txt').write_text('\n'.join(sorted(set(ffi)))+'\n')
    gy,gx=np.gradient(reference);reg=np.column_stack([reference.ravel(),gx.ravel(),gy.ravel(),np.ones(reference.size)])
    cc=np.linalg.lstsq(reg,flux.reshape(len(t),-1).T,rcond=None)[0];shifts=[-cc[1]/cc[0],-cc[2]/cc[0]]
    results=[]
    for family,ph in [('primary',0.),('secondary',.5)]:
        for mode in ['quadratic','background','pointing']:
            events=[]
            for cy in range(int(np.ceil((t.min()-E0)/P0)),int(np.floor((t.max()-E0)/P0))+1):
                center=E0+(cy+ph)*P0;use=abs(t-center)<.24;x=t[use]-center;shape=integrated_trapezoid(x,D0,ING0,dt)
                if len(x)<10 or np.sum(shape>.1)<max(2,int(.55*D0/dt)) or np.sum((shape==0)&(x<0))<max(3,int(.07/dt)) or np.sum((shape==0)&(x>0))<max(3,int(.07/dt)):continue
                if np.max(np.diff(x))>max(.04,3.01*dt):continue
                cols=[np.ones(len(x)),x,x*x]
                if mode in ['background','pointing']:
                    nuisance=[bg[use]]
                    if mode=='pointing':nuisance += [z[use] for z in shifts]
                    for cov in nuisance:
                        cov=cov-np.mean(cov)
                        if np.std(cov)>1e-8:cols.append(cov/np.std(cov))
                design=np.column_stack(cols+[-shape]);inv=np.linalg.pinv(design);ff=flux[use].reshape(len(x),-1);beta=inv@ff;res=ff-design@beta
                noise=np.maximum(np.std(res,axis=0,ddof=design.shape[1]),np.sqrt(np.mean(err[use].reshape(len(x),-1)**2,axis=0)))
                de=noise*np.sqrt(np.sum(inv[-1]**2));events.append(dict(cycle=cy,time=center,difference=beta[-1].reshape(reference.shape),error=de.reshape(reference.shape)))
            if not events:raise ValueError('No covered pixel events')
            pdiff=np.array([e['difference'] for e in events]);perr=np.array([e['error'] for e in events]);weight=1/perr**2;diff=np.sum(pdiff*weight,axis=0)/weight.sum(axis=0);error=np.maximum(np.sqrt(1/weight.sum(axis=0)),np.std(pdiff,axis=0,ddof=1)/np.sqrt(len(events)))
            suffix='' if mode=='quadratic' else '_'+mode
            np.savez_compressed(out/f'{family}_difference{suffix}.npz',reference=reference,difference=diff,error=error,per_event_difference=pdiff,per_event_error=perr,target_xy=target)
            if mode=='quadratic':save(f'{family}_events.json',[dict(cycle=e['cycle'],time=e['time']) for e in events],out)
        results.append((sector,family,out/f'{family}_difference.npz',path))
    return results


def pixels():
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astropy import units as u
    stars=pd.read_csv(ROOT/'data/stellar_systems/hidden_sources/tic260352037_astrometry.csv',dtype={'source_id':str})
    giant=stars[stars.source_id=='5496090260794096128'].iloc[0]
    items=[]
    for sector in [1,28,39,95]:items += new_pixel_cube(sector)
    rows=[];rng=np.random.default_rng(260352037);fig,axes=plt.subplots(len(items),3,figsize=(13,3.5*len(items)),layout='constrained')
    for ii,(sector,family,product,path) in enumerate(items):
        with fits.open(path) as h:w=WCS(h[2].header);mid=np.nanmedian(h[1].data['TIME'])+2457000;dt=h[1].header['TIMEDEL']
        obs=Time(mid,format='jd',scale='tdb');shifted=[]
        for _,star in stars.iterrows():
            co=SkyCoord(ra=star.ra*u.deg,dec=star.dec*u.deg,pm_ra_cosdec=star.pmra*u.mas/u.yr,pm_dec=star.pmdec*u.mas/u.yr,obstime=Time(2016.,format='jyear')).apply_space_motion(new_obstime=obs)
            shifted.append([co.ra.deg,co.dec.deg])
        starxy=w.all_world2pix(shifted,0);giantindex=stars.index[stars.source_id=='5496090260794096128'][0];host=starxy[giantindex]
        offsets={str(star.source_id):(starxy[j]-host) for j,star in stars.iterrows()}
        if product.is_dir():
            ref=np.load(product/'reference_bootstrap.npz')['reference'];z=np.load(product/'difference_quadratic.npz');diff=z['difference'];err=z['error'];pdiff=z['per_event_difference'];perr=z['per_event_error']
            modes={m:product/f'difference_{m}.npz' for m in ['quadratic','background','pointing']}
        else:
            z=np.load(product);ref=z['reference'];diff=z['difference'];err=z['error'];pdiff=z['per_event_difference'];perr=z['per_event_error'];modes={m:product.with_name(product.stem+('' if m=='quadratic' else '_'+m)+product.suffix) for m in ['quadratic','background','pointing']}
        yy,xx=np.indices(ref.shape);tests=[]
        for mode,fn in modes.items():
            zmode=np.load(fn);dd=zmode['difference'];de=zmode['error']
            for radius in [2.5,3.5,4.5]:
                roi=(xx-host[0])**2+(yy-host[1])**2<radius**2
                free,model=empirical_fit(dd,de,ref,roi)
                hypotheses={sid:empirical_fit(dd,de,ref,roi,fixed=xy)[0] for sid,xy in offsets.items()}
                # Background/pointing variants and ROI choices are retained
                # even if their likelihoods disagree: that is diagnostic.
                tests.append(dict(mode=mode,roi_radius_pixels=radius,free=free,fixed_source_hypotheses=hypotheses))
        roi=(xx-host[0])**2+(yy-host[1])**2<3.5**2;best,model=empirical_fit(diff,err,ref,roi);draw=[]
        for _ in range(200):
            idx=rng.integers(0,len(pdiff),len(pdiff));weight=1/perr[idx]**2;dd=np.sum(pdiff[idx]*weight,axis=0)/weight.sum(axis=0);de=np.maximum(np.sqrt(1/weight.sum(axis=0)),np.std(pdiff[idx],axis=0,ddof=1)/np.sqrt(len(idx)))
            b,_=empirical_fit(dd,de,ref,roi,initial=(best['dx'],best['dy']));draw.append([b['dx'],b['dy']])
        draw=np.array(draw);np.savez_compressed(OUT/f'centroid_bootstrap_s{sector}_{family}.npz',shift_draws_xy=draw)
        row=dict(sector=sector,family=family,n_events=len(pdiff),native_exposure_seconds=float(dt)*86400,host_proper_motion_corrected_xy=host,source_offsets_from_giant_xy=offsets,empirical_reference_fit=best,event_bootstrap=dict(n=200,percentile16_50_84_xy=np.percentile(draw,[16,50,84],axis=0),covariance_xy=np.cov(draw.T)),tests=tests,raw_path=str(path),raw_sha256=digest(path));rows.append(row)
        for ax,img,title in zip(axes[ii],[np.log10(np.maximum(ref,1)),diff,(diff-model)/err],[f'S{sector} {family} reference',f'Deficit: {len(pdiff)} events','Translated-reference residual / diagonal error']):
            im=ax.imshow(img,origin='lower');ax.plot(host[0],host[1],'r+',ms=10);ax.plot(host[0]+best['dx'],host[1]+best['dy'],'cx',ms=10)
            for j,star in stars.iterrows():ax.plot(starxy[j,0],starxy[j,1],'wo',mfc='none',ms=7)
            ax.set_title(title);fig.colorbar(im,ax=ax,fraction=.04)
    fig.savefig(OUT/'empirical_source_diagnostics.png',dpi=160);plt.close(fig)
    save('pixels.json',dict(tic=TIC,checked_utc=datetime.now(timezone.utc).isoformat(),rows=rows,code_sha256=digest(__file__),limitations=['Empirical reference is giant-dominated and contains faint neighbors; it is not an isolated calibrated PRF.','Subpixel interpolation ignores color, focus and intra-pixel response. Hypothesis diagonal chi-square values are sensitivity diagnostics, not decisive odds.','Whole-event bootstrap excludes reference-template systematic errors; mode/ROI variation is retained as a systematic test.','Relative empirical source offsets cancel absolute WCS bias but not proper-motion uncertainty or temporal image registration errors.','A putative unresolved companion to the giant shares the giant PSF: rejecting physical giant eclipses does not reject an unresolved faint source at that position.']))
    print(json.dumps(safe([dict(sector=r['sector'],family=r['family'],empirical_reference_fit=r['empirical_reference_fit'],bootstrap=r['event_bootstrap']) for r in rows]),indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['retrieve','photometry','pixels','morphology'],required=True);args=ap.parse_args()
    if args.mode=='retrieve':retrieve()
    elif args.mode=='photometry':photometry()
    elif args.mode=='pixels':pixels()
    elif args.mode=='morphology':
        morphology();timing_stability();wide_pair_dynamics()
