"""Native SPOC TPF source tests for the already-public TIC49617264 TCE.

No discovery/planet interpretation is inferred from catalog absence or a dip.
Public products are retained, raw integrity is verified, and compact evidence
is written under candidates/tce107/tic49617264_pixels only.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,os,urllib.parse,urllib.request
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.special import erf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import ROOT,verify_fits_checksums

TIC=49617264
OUT=ROOT/'candidates/tce107/tic49617264_pixels'
CACHE=Path('/workspace/tess-research-data/tce107/pixels')
P,T0,DURATION=4.83194208145142,4265.644028635728,.09524246875
RA,DEC=37.8931844608584,-40.0058723004587


def safe(v):
    if isinstance(v,dict):return {str(k):safe(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):return [safe(x) for x in v]
    if isinstance(v,np.generic):return safe(v.item())
    if isinstance(v,float) and not np.isfinite(v):return None
    return v


def template(t,dt,kind='trapezoid'):
    off=np.linspace(-dt/2,dt/2,15)
    ph=(t[:,None]+off-T0+P/2)%P-P/2
    if kind=='box_broad':return np.mean(abs(ph)<.11544541666666668/2,axis=1)
    if kind=='box_core':return np.mean(abs(ph)<.065/2,axis=1)
    return np.mean(np.clip((DURATION/2-abs(ph))/(DURATION*.20),0,1),axis=1)


def fit_events(t,f,err,bg,px,py,dt,mode):
    events=[];shape=template(t,dt,mode if mode.startswith('box') else 'trapezoid')
    for cycle in range(int(np.ceil((t.min()-T0)/P)),int(np.floor((t.max()-T0)/P))+1):
        center=T0+P*cycle;use=abs(t-center)<.40;tx=t[use]-center;s=shape[use]
        if np.sum(s>.05)<15 or np.sum((s==0)&(tx<0))<30 or np.sum((s==0)&(tx>0))<30:continue
        cols=[np.ones(len(tx)),tx,tx**2]
        if mode not in ('quadratic',):
            for c in [bg[use],px[use],py[use]]:
                c=c-np.mean(c)
                if np.std(c)>1e-10:cols.append(c/np.std(c))
        a=np.column_stack(cols+[-s]);inv=np.linalg.pinv(a);ff=f[use].reshape(len(tx),-1);coef=inv@ff;res=ff-a@coef
        sigma=np.maximum(np.std(res,axis=0,ddof=a.shape[1]),np.sqrt(np.mean(err[use].reshape(len(tx),-1)**2,axis=0)))
        ee=sigma*np.sqrt(np.sum(inv[-1]**2))
        events.append({'cycle':cycle,'center_btjd':center,'n_in':int((s>.05).sum()),'n_oot':int((s==0).sum()),
                       'difference':coef[-1].reshape(f.shape[1:]),'error':ee.reshape(f.shape[1:]),
                       'residual':res,'projection':inv[-1],'time_btjd':t[use],'model_condition':np.linalg.cond(a),
                       'a':a,'coef':coef,'s':s,'flux':ff,'shape':f.shape[1:]})
    return events


def combine(events):
    d=np.array([e['difference'] for e in events]);ee=np.array([e['error'] for e in events]);w=1/ee**2
    return np.sum(w*d,axis=0)/np.sum(w,axis=0),np.sqrt(1/np.sum(w,axis=0))


def gaussfit(image,error,seed,roi=None,fixed_center=False):
    yy,xx=np.indices(image.shape)
    if roi is None:roi=np.ones(image.shape,bool)
    x,y=xx[roi],yy[roi];e=np.maximum(error[roi],.03)
    def model(a):
        flux,cx,cy,sx,sy,b,bx,by=a
        gx=.5*(erf((x+.5-cx)/(np.sqrt(2)*sx))-erf((x-.5-cx)/(np.sqrt(2)*sx)))
        gy=.5*(erf((y+.5-cy)/(np.sqrt(2)*sy))-erf((y-.5-cy)/(np.sqrt(2)*sy)))
        return flux*gx*gy+b+bx*(x-np.mean(x))+by*(y-np.mean(y))
    start=[max(1,float(np.nanmax(image[roi])))*5,*seed,.65,.65,0,0,0]
    lower=[0,0,0,.2,.2,-np.inf,-np.inf,-np.inf];upper=[np.inf,image.shape[1]-1,image.shape[0]-1,2.5,2.5,np.inf,np.inf,np.inf]
    if fixed_center:
        lower[1:3]=[seed[0]-1e-9,seed[1]-1e-9];upper[1:3]=[seed[0]+1e-9,seed[1]+1e-9]
    fit=least_squares(lambda a:(model(a)-image[roi])/e,start,bounds=(lower,upper),max_nfev=400)
    a=fit.x
    return {'source_flux_e_s':a[0],'x':a[1],'y':a[2],'sigma_x':a[3],'sigma_y':a[4],
            'chi2_diagonal':np.sum(fit.fun**2),'dof':int(roi.sum())-len(a),'width_bound_hit':bool(min(a[3:5])<.2001 or max(a[3:5])>2.4999)},a


def boot_centroids(events,seed,n):
    rng=np.random.default_rng(49617264);loc=[]
    for _ in range(n):
        ev=[]
        for j in rng.integers(0,len(events),len(events)):
            e=events[j];res=e['residual'];length=len(res)
            starts=rng.integers(0,length,(length+9)//10)
            inds=np.concatenate([(start+np.arange(10))%length for start in starts])[:length]
            noise=(e['projection']@res[inds]).reshape(e['difference'].shape)
            ev.append({'difference':e['difference']+noise,'error':e['error']})
        d,de=combine(ev);g,_=gaussfit(d,de,seed);loc.append([g['x'],g['y']])
    return np.array(loc)


def photometry(events,positions,reference,optimal):
    yy,xx=np.indices(reference.shape);results=[];arrays={}
    masks={'spoc_optimal':optimal}
    for label,position in positions.items():
        for radius in [1.,1.5,2.]:masks[f'{label}_r{radius}']=(xx-position[0])**2+(yy-position[1])**2<radius**2
    for name,mask in masks.items():
        measures=[];times=[];fluxes=[]
        for e in events:
            if not mask.any():continue
            selector=mask.ravel();model=e['a'];raw=e['flux'][:,selector].sum(axis=1);coef=np.linalg.lstsq(model,raw,rcond=None)[0]
            residual=raw-model@coef;n=len(raw);scale=np.std(residual,ddof=model.shape[1]);ee=scale*np.sqrt(np.sum(np.linalg.pinv(model)[-1]**2))
            baseline=float(np.median(model[:,:-1]@coef[:-1]));dep=coef[-1]/baseline;err=ee/abs(baseline)
            measures.append({'cycle':e['cycle'],'depth':dep,'error_residual_scaled':err,'snr':dep/err,'baseline_e_s':baseline})
            times.extend(e['time_btjd']);fluxes.extend((raw-model[:,:-1]@coef[:-1])/baseline+1)
        weights=1/np.array([m['error_residual_scaled']**2 for m in measures]);deps=np.array([m['depth'] for m in measures]);mean=np.average(deps,weights=weights)
        formal=np.sqrt(1/weights.sum());scatter=np.std(deps,ddof=1)/np.sqrt(len(deps));err=max(formal,scatter)
        results.append({'aperture':name,'pixels':int(mask.sum()),'depth':mean,'error_max_diagonal_event_scatter':err,'snr':mean/err,
                        'events':measures,'limitation':'Diluted observed-aperture depth; correlated residual covariance incomplete; no physical model.'})
        arrays[name]=[np.array(times),np.array(fluxes)]
    return results,arrays


def analyze(sector,ndraws,gaia):
    out=OUT/f'sector{sector}';out.mkdir(parents=True,exist_ok=True)
    obs=[r for r in json.loads((ROOT/'data/tce107_extra_observations.json').read_text()) if r['tic']==TIC and r['sector']==sector][0]
    name=obs['data_uri'].rsplit('/',1)[1].replace('_lc.fits','_tp.fits');raw=CACHE/name
    uri=obs['data_uri'].replace('_lc.fits','_tp.fits');url='https://mast.stsci.edu/api/v0.1/Download/file?'+urllib.parse.urlencode({'uri':uri})
    if not raw.exists():
        CACHE.mkdir(parents=True,exist_ok=True)
        with urllib.request.urlopen(url,timeout=90) as response,raw.with_suffix('.part').open('wb') as stream:
            while chunk:=response.read(1024**2):stream.write(chunk)
        raw.with_suffix('.part').replace(raw)
    byte_hash=hashlib.sha256(raw.read_bytes()).hexdigest()
    known_manifest=ROOT/'data/tce107_pixel_manifest.json'
    if known_manifest.exists():
        expected=[r for r in json.loads(known_manifest.read_text()) if r['data_uri']==uri]
        if expected and expected[0]['sha256']!=byte_hash:
            raise ValueError('Raw TPF differs from recorded SHA256; retain mismatch for diagnosis rather than rewriting evidence.')
    with fits.open(raw) as h:
        integrity=verify_fits_checksums(raw,h);header=dict(h[0].header);lh=dict(h[1].header);wcs=WCS(h[2].header)
        data=h[1].data;mask=(h[2].data&2)>0
        time=np.array(data['TIME']);f=np.array(data['FLUX'],float);err=np.array(data['FLUX_ERR'],float);q=np.array(data['QUALITY'])
        # All-NaN quality-flagged rows are discarded before evaluating median.
        braw=np.array(data['FLUX_BKG'],float);bg=np.full(len(time),np.nan)
        bgfinite=np.all(np.isfinite(braw),axis=(1,2));bg[bgfinite]=np.median(braw[bgfinite],axis=(1,2))
        px=np.array(data['POS_CORR1'],float);py=np.array(data['POS_CORR2'],float)
        keep=(q==0)&np.isfinite(time)&np.all(np.isfinite(f)&np.isfinite(err),axis=(1,2))&np.isfinite(bg)&np.isfinite(px)&np.isfinite(py)
        quality=dict(zip(*np.unique(q,return_counts=True)));original_n=len(time);time,f,err,bg,px,py=[a[keep] for a in [time,f,err,bg,px,py]]
    reference=np.median(f,axis=0);yy,xx=np.indices(reference.shape);target=wcs.all_world2pix([[RA,DEC]],0)[0]
    stars=gaia.copy() if gaia is not None else pd.DataFrame(columns=['ra','dec','phot_g_mean_mag'])
    if len(stars):stars[['x','y']]=wcs.all_world2pix(stars[['ra','dec']].to_numpy(),0);stars.to_csv(out/'gaia_pixel_positions.csv',index=False)
    bright_y,bright_x=np.unravel_index(np.argmax(reference),reference.shape)
    # Median-image brightest component is a comparison source, not assumed to
    # be the host. Gaia attribution remains separate from pixel-peak labeling.
    positions={'target':target,'bright_pixel_peak':np.array([bright_x,bright_y])}
    if len(stars):
        for _,s in stars[stars.phot_g_mean_mag<15].iterrows():
            if 0<=s.x<reference.shape[1] and 0<=s.y<reference.shape[0] and np.hypot(s.x-target[0],s.y-target[1])>.5:
                positions[f"gaia_{s.source_id}"]=np.array([s.x,s.y])
    modes={};evs={}
    for mode in ['quadratic','pointing','box_broad','box_core']:
        events=fit_events(time,f,err,bg,px,py,float(lh['TIMEDEL']),mode);evs[mode]=events;d,de=combine(events)
        maxy,maxx=np.unravel_index(np.argmax(d/de),d.shape)
        gf,_=gaussfit(d,de,[maxx,maxy]);gf['target_distance_pixels']=np.hypot(gf['x']-target[0],gf['y']-target[1]);gf['sky_ra_dec']=wcs.all_pix2world([[gf['x'],gf['y']]],0)[0]
        modes[mode]={'n_events':len(events),'events':[{k:e[k] for k in ['cycle','center_btjd','n_in','n_oot','model_condition']} for e in events], 'integrated_gaussian_source_fit':gf}
        np.savez_compressed(out/f'difference_{mode}.npz',difference=d,error=de,event_differences=np.array([e['difference'] for e in events]),event_errors=np.array([e['error'] for e in events]))
    gf=modes['pointing']['integrated_gaussian_source_fit'];events=evs['pointing'];bootstrap=boot_centroids(events,[gf['x'],gf['y']],ndraws)
    bsum={'n':len(bootstrap),'median_xy':np.median(bootstrap,axis=0),'percentile_16_84_xy':np.percentile(bootstrap,[16,84],axis=0),
          'covariance':np.cov(bootstrap.T),'target_offset_pixels_16_50_84':np.percentile(np.linalg.norm(bootstrap-target,axis=1),[16,50,84]),
          'limitation':'Whole event plus synchronized 10-cadence residual block bootstrap; Gaussian/PRF, WCS, proper-motion and unresolved-blend uncertainties excluded.'}
    if len(stars):
        stars['distance_difference_pixels']=np.hypot(stars.x-gf['x'],stars.y-gf['y']);stars['distance_target_pixels']=np.hypot(stars.x-target[0],stars.y-target[1]);stars.to_csv(out/'gaia_pixel_positions.csv',index=False)
    source_hypotheses={};dd=np.load(out/'difference_pointing.npz')
    for name,position in positions.items():
        fg,_=gaussfit(dd['difference'],dd['error'],position,fixed_center=True)
        fg['delta_chi2_vs_free_source']=fg['chi2_diagonal']-gf['chi2_diagonal']
        fg['bootstrap_fraction_closer_target_than_this_position']=None if name=='target' else np.mean(np.linalg.norm(bootstrap-target,axis=1)<np.linalg.norm(bootstrap-position,axis=1))
        source_hypotheses[name]=fg
    aperture,arrays=photometry(events,positions,reference,mask)
    pd.DataFrame([{'aperture':a['aperture'],**e} for a in aperture for e in a['events']]).to_csv(out/'per_event_aperture_depths.csv',index=False)
    np.savez_compressed(out/'photometry_arrays.npz',time_btjd=time,background=bg,pointing_dx=px,pointing_dy=py,reference=reference,optimal_aperture=mask,bootstrap_xy=bootstrap,
                        **{f'{key}_{label}':v for key,vals in arrays.items() for label,v in zip(['time_btjd','flux'],vals)})
    row={'tic':TIC,'sector':sector,'obsid':obs['obsid'],'obs_id':obs['obs_id'],'data_uri':uri,'retrieval_url':url,'local_path':str(raw),
         'sha256':byte_hash,'bytes':raw.stat().st_size,'processing_version':header['PROCVER'],'time_system':lh['TIMESYS'],
         'bjd_reference':lh['BJDREFI']+lh['BJDREFF'],'cadence_days':lh['TIMEDEL'],'camera':header['CAMERA'],'ccd':header['CCD'],'fits_checksums':integrity,
         'checked_utc':datetime.now(timezone.utc).isoformat()}
    summary={'tic':TIC,'sector':sector,'ephemeris':{'period_days':P,'epoch_btjd':T0,'duration_days':DURATION,'note':'Public TCE period, initial epoch shifted ~6 minutes following native sector107 fine BLS; fixed trapezoid with20%ingress.'},
             'input':row,'target_xy':target,'original_cadences':original_n,'usable_cadences':len(time),'quality_distribution':quality,
             'difference_modes':modes,'bootstrap':bsum,'apertures':aperture,'crowdsap':lh.get('CROWDSAP'),'flfrcsap':lh.get('FLFRCSAP'),
             'background_subtracted_in_original':lh.get('BACKAPP'),'fixed_center_gaussian_source_hypotheses':source_hypotheses,
             'assessment':'source_favors_TIC49617264_neighborhood_already_ground_followed_TOI2417.11_candidate_no_new_discovery',
             'prior_recognition':{'reported_aliases':['TOI2417.11','TIC49617264.11'],'ground_report_date':'2024-10-05',
                                 'notes_path':'data/broad_catalogs/tic49617264_followup20241005_notes.txt',
                                 'note':'Ground report explicitly identifies NPC to original bright TIC49617263/TOI2417.01 and on-target 13.1ppt in 6.6arcsec aperture; separately cites TRAPPIST 2024-09-15 12.8ppt in 4.5arcsec aperture.'},
             'limitations':['Public TCE already recognizes the signal. No independent discovery claim.',
                            'Gaussian is an empirical source-position diagnostic, not an official TESS PRF model.',
                            'Gaia positions not proper-motion propagated; unresolved companion/source incomplete.',
                            'Pixel extractions share photons with SPOC light curves; not independent significance.']}
    (out/'pixel_vetting.json').write_text(json.dumps(safe(summary),indent=2)+'\n')
    dd=np.load(out/'difference_pointing.npz');fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,im,title in zip(axes,[np.log10(np.maximum(reference,.1)),dd['difference'],dd['difference']/dd['error']],['Median image log10 e-/s','Local transit difference e-/s','Difference/diagonal error']):
        limit=None if title.startswith('Median') else np.percentile(abs(im),98)
        hh=ax.imshow(im,origin='lower',cmap='viridis' if limit is None else 'RdBu_r',vmin=None if limit is None else -limit,vmax=limit)
        ax.plot(*target,'gx',ms=12);ax.plot(gf['x'],gf['y'],'m+',ms=15)
        if len(stars):
            for _,s in stars[stars.phot_g_mean_mag<15].iterrows():
                if 0<=s.x<im.shape[1] and 0<=s.y<im.shape[0]:ax.plot(s.x,s.y,'o',mfc='none',mec='orange',ms=8)
        ax.set_title(title);fig.colorbar(hh,ax=ax,shrink=.65)
    fig.suptitle(f'TIC{TIC} sector{sector}: green target, magenta difference, orange Gaia G<15');fig.tight_layout();fig.savefig(out/'source_localization.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(int(np.ceil((len(events)+1)/3)),3,figsize=(12,4*int(np.ceil((len(events)+1)/3))),squeeze=False)
    lim=np.percentile(abs(np.array([e['difference'] for e in events])),98)
    for ax,e in zip(axes.ravel(),events):
        ax.imshow(e['difference'],origin='lower',cmap='RdBu_r',vmin=-lim,vmax=lim);ax.plot(*target,'gx',ms=10);ax.set_title(f"Cycle{e['cycle']},BTJD{e['center_btjd']:.4f}")
    for ax in axes.ravel()[len(events):-1]:ax.axis('off')
    ax=axes.ravel()[-1];ax.scatter(bootstrap[:,0],bootstrap[:,1],s=4,alpha=.2);ax.plot(*target,'gx',ms=10)
    for name,pos in positions.items():
        if name!='target':ax.plot(*pos,'o',mfc='none',mec='orange')
    ax.set_title('Bootstrap source locations');fig.tight_layout();fig.savefig(out/'event_source_images.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(len(arrays),1,figsize=(9,2.7*len(arrays)),squeeze=False)
    for ax,(name,(at,af)) in zip(axes.ravel(),arrays.items()):
        ph=(at-T0+P/2)%P-P/2;ax.plot(ph*24,af,'.',ms=2,alpha=.3);ax.set_xlim(-5,5)
        r=[r for r in aperture if r['aperture']==name][0];ax.set_title(f"{name}: {r['depth']*1e3:.2f} ppt, approx SNR{r['snr']:.1f}")
    axes[-1,0].set_xlabel('Hours from event');fig.tight_layout();fig.savefig(out/'aperture_photometry.png',dpi=140);plt.close(fig)
    return summary


def main():
    global P,T0,DURATION
    parser=argparse.ArgumentParser();parser.add_argument('--sectors',type=int,nargs='+',default=[106,107]);parser.add_argument('--bootstrap',type=int,default=500)
    parser.add_argument('--period',type=float,default=P);parser.add_argument('--epoch',type=float,default=T0);parser.add_argument('--duration',type=float,default=DURATION)
    parser.add_argument('--gaia',type=Path);args=parser.parse_args();P,T0,DURATION=args.period,args.epoch,args.duration
    OUT.mkdir(parents=True,exist_ok=True)
    gaia=pd.read_csv(args.gaia,dtype={'source_id':str}) if args.gaia else None
    manifest=ROOT/'data/tce107_pixel_manifest.json';records=json.loads(manifest.read_text()) if manifest.exists() else []
    for sector in args.sectors:
        s=analyze(sector,args.bootstrap,gaia);records=[r for r in records if r['data_uri']!=s['input']['data_uri']]+[s['input']]
        manifest.write_text(json.dumps(safe(records),indent=2)+'\n')
        print(json.dumps(safe({'sector':sector,'target_xy':s['target_xy'],'gaussian':s['difference_modes']['pointing']['integrated_gaussian_source_fit'],'bootstrap':s['bootstrap'],
                               'apertures':[{k:a[k] for k in ['aperture','depth','error_max_diagonal_event_scatter','snr']} for a in s['apertures']]}),indent=2),flush=True)


if __name__=='__main__':main()
