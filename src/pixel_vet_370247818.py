"""Independent FFI source-localization checks for TIC 370247818.

Retains per-event images, alternate detrending, independent-aperture photometry,
and residual/event bootstraps. Gaussian localization is a diagnostic, not an
official TESS PRF fit, and cannot discriminate unresolved companions.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import argparse
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.timeseries import BoxLeastSquares
from scipy.optimize import least_squares
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'candidates/broad_s1/tic370247818_pixels'
RAW = Path('/workspace/tess-research-data/broad/pixels/tess-s0001-3-3_55.455693_-79.113281_15x15_astrocut.fits')
P, T0, DUR, INGRESS = 3.5225004449826516, 1326.915994193442, .0939041044544696, .33111363574530256
CADENCE_DAYS=.0208333333333333
Y, X = np.indices((15,15))


def safe(x):
    if isinstance(x,dict): return {str(k):safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [safe(v) for v in x]
    if isinstance(x,np.ndarray): return safe(x.tolist())
    if isinstance(x,np.generic): return safe(x.item())
    if isinstance(x,float) and not np.isfinite(x): return None
    return x


def shape(t):
    # Finite exposure average over 30-min cadence span; deadtime is distributed
    # throughout the cadence so EXPOSURE does not describe a single short sample.
    offsets=np.linspace(-CADENCE_DAYS/2,CADENCE_DAYS/2,31)
    ph=(t[:,None]+offsets-T0+P/2)%P-P/2
    ingress=DUR*INGRESS
    return np.mean(np.clip((DUR/2-np.abs(ph))/ingress,0,1),axis=1)


def gaussian_fit(im, err, seed, roi=None):
    if roi is None: roi=np.ones(im.shape,bool)
    # Elliptical Gaussian + affine background: empirical localization only.
    xx,yy=X[roi],Y[roi]
    e=np.maximum(err[roi],.03)
    def model(a):
        amp,cx,cy,sx,sy,b,bx,by=a
        return amp*np.exp(-.5*((xx-cx)/sx)**2-.5*((yy-cy)/sy)**2)+b+bx*(xx-7)+by*(yy-7)
    a0=[max(.2,float(np.nanmax(im[roi]))),seed[0],seed[1],.85,.85,0,0,0]
    lo=[0,0,0,.4,.4,-np.inf,-np.inf,-np.inf]
    hi=[np.inf,14,14,2.8,2.8,np.inf,np.inf,np.inf]
    fit=least_squares(lambda a:(model(a)-im[roi])/e,a0,bounds=(lo,hi),max_nfev=300)
    a=fit.x
    return {'amplitude_e_s':a[0],'x':a[1],'y':a[2],'sigma_x':a[3],'sigma_y':a[4],
            'background':a[5],'chi2_diagonal':np.sum(fit.fun**2),'dof':len(xx)-len(a),'success':fit.success}, a


def pointing(flux, ref, exclude):
    # Amplitude freedom absorbs brightness changes of the bright registration
    # star. Excluding target pixels prevents target transits influencing shifts.
    gy,gx=np.gradient(ref)
    use=~exclude
    design=np.column_stack([ref[use],gx[use],gy[use],np.ones(use.sum())])
    coeff=np.linalg.lstsq(design,flux[:,use].T,rcond=None)[0]
    return -coeff[1]/coeff[0],-coeff[2]/coeff[0]


def event_images(t,flux,err,dx,dy,bg,mode,halfwindow=.30):
    first=int(np.ceil((t.min()-T0)/P));last=int(np.floor((t.max()-T0)/P))
    events=[]
    for cycle in range(first,last+1):
        center=T0+cycle*P
        use=abs(t-center)<halfwindow
        tt=t[use]-center;s=shape(t[use]);ff=flux[use].reshape(-1,225)
        if np.sum(s>.05)<2 or np.sum(s==0)<12: continue
        nuis=[np.ones(len(tt)),tt,tt**2]
        if mode in ('background','pointing'):
            cov=bg[use]-np.mean(bg[use]);nuis.append(cov/(np.std(cov)+1e-10))
        if mode=='pointing':
            for cov in [dx[use],dy[use]]:
                cov=cov-np.mean(cov);nuis.append(cov/(np.std(cov)+1e-10))
        design=np.column_stack(nuis+[-s])
        inv=np.linalg.pinv(design)
        coef=inv@ff
        residual=ff-design@coef
        sigma=np.maximum(np.std(residual,axis=0,ddof=design.shape[1]),np.sqrt(np.mean(err[use].reshape(-1,225)**2,axis=0)))
        differr=sigma*np.sqrt(np.sum(inv[-1]**2))
        events.append({'cycle':cycle,'center_btjd':center,'n_in':int(np.sum(s>.05)),
                       'n_out':int(np.sum(s==0)),'n':len(tt),'diff':coef[-1].reshape(15,15),
                       'err':differr.reshape(15,15),'residual':residual,'projection':inv[-1],
                       'tt':tt,'design_condition':np.linalg.cond(design)})
    return events


def combine(events):
    d=np.array([e['diff'] for e in events]);ee=np.array([e['err'] for e in events]);w=1/ee**2
    return np.sum(d*w,axis=0)/np.sum(w,axis=0),np.sqrt(1/np.sum(w,axis=0))


def aperture_photometry(t,flux,err,bg,dx,dy,positions):
    result={};arrays={}
    for label,pos in positions.items():
        for radius in [1.,1.5,2.]:
            mask=(X-pos[0])**2+(Y-pos[1])**2<radius**2
            values=np.sum(flux[:,mask],axis=1)-mask.sum()*bg
            errors=np.sqrt(np.sum(err[:,mask]**2,axis=1))
            measurements=[];loc_t=[];loc_f=[];loc_e=[]
            for cy in range(int(np.ceil((t.min()-T0)/P)),int(np.floor((t.max()-T0)/P))+1):
                center=T0+cy*P;use=abs(t-center)<.30;tt=t[use]-center;s=shape(t[use]);v=values[use];ee=errors[use]
                if np.sum(s>.05)<2 or np.sum(s==0)<12:continue
                med=np.median(v[s==0]);yy=v/med;ee=ee/med
                # Simultaneous nuisance terms preserve the transit response.
                cols=[np.ones(len(tt)),tt,tt**2]
                for cov in [dx[use],dy[use]]:
                    cov=cov-np.mean(cov);cols.append(cov/(np.std(cov)+1e-10))
                a=np.column_stack(cols+[-s]);w=1/ee**2;cov=np.linalg.pinv(a.T@(w[:,None]*a));coef=cov@(a.T@(w*yy));res=yy-a@coef
                scale=max(1,float(np.sum(w*res**2)/(len(tt)-a.shape[1])))
                de=np.sqrt(cov[-1,-1]*scale)
                measurements.append({'cycle':cy,'depth':float(coef[-1]),'error':float(de),'snr_diagonal':float(coef[-1]/de),'n_in':int((s>.05).sum())})
                loc_t.extend(t[use]);loc_f.extend(yy-a[:,:-1]@coef[:-1]+1);loc_e.extend(ee*np.sqrt(scale))
            w=np.array([1/m['error']**2 for m in measurements]);d=np.array([m['depth'] for m in measurements])
            dep=np.sum(w*d)/np.sum(w);de=max(np.sqrt(1/np.sum(w)),float(np.std(d,ddof=1)/np.sqrt(len(d))))
            key=f'{label}_r{radius}'
            result[key]={'center':pos,'radius_pixels':radius,'n_pixels':int(mask.sum()),'median_flux_bg_sub_e_s':float(np.median(values)),
                         'depth':dep,'error_max_diagonal_event_scatter':de,'snr':dep/de,'events':measurements,
                         'limitation':'Uncorrected aperture dilution; empirical pointing regression; no physical transit or PRF model.'}
            arrays[key]=[np.array(loc_t),np.array(loc_f),np.array(loc_e)]
    return result,arrays


def aperture_period_check(t,flux,err,bg,positions):
    """BLS on independently extracted all-sector light curves, with no event
    window selection. Regional search is a period consistency check only.
    """
    result={};curves={}
    for label,pos in positions.items():
        if label=='intermediate':continue
        mask=(X-pos[0])**2+(Y-pos[1])**2<2.0**2
        f=np.sum(flux[:,mask],axis=1)-mask.sum()*bg
        e=np.sqrt(np.sum(err[:,mask]**2,axis=1));flat=np.empty_like(f)
        boundaries=np.r_[0,np.flatnonzero(np.diff(t)>.10)+1,len(t)]
        for lo,hi in zip(boundaries[:-1],boundaries[1:]):
            width=min(49,(hi-lo)//2*2-1)
            trend=median_filter(f[lo:hi],size=max(3,width),mode='nearest')
            flat[lo:hi]=f[lo:hi]/trend
        ef=e/np.median(f)
        curves[label]=(flat,ef)
    neighbor_flat=curves['bright_neighbor'][0]
    neighbor_sigma=1.4826*np.median(abs(neighbor_flat-np.median(neighbor_flat)))
    # An ephemeris-independent bright comparison-star outlier criterion rejects
    # rare FFI image disruptions that QUALITY=0 misses. It must not be treated
    # as safe for a survey of eclipses originating on the comparison source.
    common_keep=abs(neighbor_flat-np.median(neighbor_flat))<6*neighbor_sigma
    for label,(flat,ef) in curves.items():
        noise=max(float(np.median(abs(np.diff(flat)-np.median(np.diff(flat))))*1.4826/np.sqrt(2)),float(np.median(ef)))
        ef=np.maximum(ef,noise)
        result[label]={}
        for mode,keep in [('all_quality_zero',np.ones(len(t),bool)),('comparison_star_image_outlier_mask',common_keep)]:
            model=BoxLeastSquares(t[keep],flat[keep],ef[keep])
            bls=model.power(np.linspace(2.8,4.3,4001),[.035,.055,.075,.095],oversample=10,objective='snr')
            j=int(np.argmax(bls.power))
            result[label][mode]={'period_days':float(bls.period[j]),'epoch_btjd':float(bls.transit_time[j]),'duration_days':float(bls.duration[j]),
                           'depth':float(bls.depth[j]),'depth_err_diagonal':float(bls.depth_err[j]),'snr_diagonal':float(bls.depth_snr[j]),
                           'noise_floor':noise,'search_range_days':[2.8,4.3],'n_periods':4001,'n_kept':int(keep.sum()),
                           'limitation':'Regional consistency search, not independent significance test; median detrending uncorrected pointing and event-correlated noise.'}
            np.savez_compressed(OUT/f'aperture_bls_{label}_{mode}.npz',time_btjd=t,flattened_flux=flat,relative_error=ef,keep=keep,period=bls.period,power=bls.power)
    result['comparison_star_mask']={'sigma_threshold':6,'mad_sigma':neighbor_sigma,'n_removed':int((~common_keep).sum()),
                                  'n_candidate_in_transit_removed':int(np.sum((~common_keep)&(shape(t)>.05))),
                                  'removed_times_btjd':t[~common_keep],
                                  'limitation':'Useful only for testing target signal, not for searching the chosen bright comparison source; all original data retained.'}
    return result


def main():
    global RAW,OUT,P,T0,CADENCE_DAYS
    parser=argparse.ArgumentParser();parser.add_argument('--bootstrap',type=int,default=500)
    parser.add_argument('--sector',type=int,default=1);parser.add_argument('--raw',type=Path)
    parser.add_argument('--period',type=float,default=P);parser.add_argument('--epoch',type=float,default=T0)
    args=parser.parse_args();P,T0=args.period,args.epoch
    if args.raw:RAW=args.raw
    elif args.sector!=1:
        matching=list(RAW.parent.glob(f'tess-s{args.sector:04d}-*_55.455693_-79.113281_15x15_astrocut.fits'))
        if len(matching)!=1:raise RuntimeError(f'Expected one sector{args.sector} cutout, found {matching}')
        RAW=matching[0]
    if args.sector!=1:OUT=OUT/f'sector{args.sector}'
    OUT.mkdir(parents=True,exist_ok=True)
    with fits.open(RAW,checksum=True) as hdus:
        header=dict(hdus[0].header);wcs=WCS(hdus[2].header)
        CADENCE_DAYS=float(hdus[1].header['TIMEDEL'])
        checksum=[{'hdu':i,'checksum':h.verify_checksum(),'datasum':h.verify_datasum()} for i,h in enumerate(hdus)]
        if any(c['checksum']==0 or c['datasum']==0 for c in checksum):
            raise RuntimeError('FITS integrity verification failed; preserve raw file and stop this analysis.')
        data=hdus[1].data
        ffi_files=np.array(data['FFI_FILE']).astype(str)
        science_digest=hashlib.sha256()
        for column,dtype in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4'),('FLUX_ERR','<f4')]:
            science_digest.update(column.encode('ascii'))
            science_digest.update(np.asarray(data[column],dtype=dtype).tobytes())
        science_digest.update(np.asarray(ffi_files,dtype='S38').tobytes())
        t=np.array(data['TIME']);q=np.array(data['QUALITY']);f=np.array(data['FLUX'],float);err=np.array(data['FLUX_ERR'],float)
        keep=(q==0)&np.isfinite(t)&np.all(np.isfinite(f),axis=(1,2))&np.all(np.isfinite(err),axis=(1,2))
        counts={'total':len(t),'quality_zero_and_finite':int(keep.sum()),'quality_distribution':dict(zip(*np.unique(q,return_counts=True)))}
        t,f,err=t[keep],f[keep],err[keep]
    (OUT/'ffi_product_ids.txt').write_text('\n'.join(ffi_files)+'\n')
    stars=pd.read_csv(ROOT/'data/broad_catalogs/tic370247818_gaia_nearby.csv',dtype={'source_id':str})
    stars[['x','y']]=wcs.all_world2pix(stars[['ra','dec']].to_numpy(),0)
    stars.to_csv(OUT/'gaia_pixel_positions.csv',index=False)
    target=stars.loc[stars.source_id=='4625697305555206912',['x','y']].to_numpy()[0]
    neighbor=stars.loc[stars.source_id=='4625697168116255360',['x','y']].to_numpy()[0]
    third=stars.loc[stars.source_id=='4625697717872066432',['x','y']].to_numpy()[0]
    pos={'target':target,'bright_neighbor':neighbor,'third_star':third,'intermediate':(target+neighbor)/2}
    ref=np.median(f,axis=0)
    # Blank sky excludes catalog stars brighter than G=19 and all edge rows.
    blank=np.ones((15,15),bool)
    for _,star in stars[stars.phot_g_mean_mag<19].iterrows(): blank&=(X-star.x)**2+(Y-star.y)**2>3.2**2
    bg=np.median(f[:,blank],axis=1)
    exclude=(X-target[0])**2+(Y-target[1])**2<2.5**2
    dx,dy=pointing(f,ref,exclude)
    modes={};event_sets={}
    for mode in ['quadratic','background','pointing']:
        ev=event_images(t,f,err,dx,dy,bg,mode);event_sets[mode]=ev;diff,de=combine(ev)
        modes[mode]={'n_events':len(ev),'events':[{k:v for k,v in e.items() if k in ['cycle','center_btjd','n_in','n_out','n','design_condition']} for e in ev]}
        # Source-centered ROI restricts bright-neighbor dipole artifacts.
        roi=(X-target[0])**2+(Y-target[1])**2<3.0**2
        gf,_=gaussian_fit(diff,de,target,roi)
        gf['offset_target_pixels']=np.hypot(gf['x']-target[0],gf['y']-target[1])
        gf['offset_neighbor_pixels']=np.hypot(gf['x']-neighbor[0],gf['y']-neighbor[1])
        gf['sky_ra_dec']=wcs.all_pix2world([[gf['x'],gf['y']]],0)[0]
        modes[mode]['target_roi_gaussian_fit']=gf
        gf2,_=gaussian_fit(diff,de,neighbor,(X-neighbor[0])**2+(Y-neighbor[1])**2<2.5**2)
        modes[mode]['neighbor_roi_gaussian_fit']=gf2
        np.savez_compressed(OUT/f'difference_{mode}.npz',difference=diff,error=de,
                            per_event_difference=np.array([e['diff'] for e in ev]),per_event_error=np.array([e['err'] for e in ev]))
    results,aparrays=aperture_photometry(t,f,err,bg,dx,dy,pos)
    period_checks=aperture_period_check(t,f,err,bg,pos)
    sensitivity=[]
    dd=np.load(OUT/'difference_pointing.npz')
    for radius in [2.,2.2,2.5,3.]:
        rroi=(X-target[0])**2+(Y-target[1])**2<radius**2
        rg,_=gaussian_fit(ref,np.ones((15,15)),target,rroi)
        dg,_=gaussian_fit(dd['difference'],dd['error'],target,rroi)
        sensitivity.append({'radius_pixels':radius,'n_pixels':int(rroi.sum()),'reference_x':rg['x'],'reference_y':rg['y'],
                            'difference_x':dg['x'],'difference_y':dg['y'],'reference_sigma_x':rg['sigma_x'],'reference_sigma_y':rg['sigma_y'],
                            'difference_minus_reference_offset_pixels':np.hypot(dg['x']-rg['x'],dg['y']-rg['y']),
                            'reference_width_bound_hit':bool(min(rg['sigma_x'],rg['sigma_y'])<.4001 or max(rg['sigma_x'],rg['sigma_y'])>2.7999),
                            'difference_width_bound_hit':bool(min(dg['sigma_x'],dg['sigma_y'])<.4001 or max(dg['sigma_x'],dg['sigma_y'])>2.7999)})
    pd.DataFrame(safe(sensitivity)).to_csv(OUT/'localization_roi_sensitivity.csv',index=False)
    boot=[];rng=np.random.default_rng(370247818);ev=event_sets['pointing'];roi=(X-target[0])**2+(Y-target[1])**2<3.0**2
    # Resample events and residuals in 3-cadence blocks, synchronously across all
    # pixels, retaining pixel covariance and short-time correlations.
    for _ in range(args.bootstrap):
        replicate=[]
        for j in rng.integers(0,len(ev),len(ev)):
            e=ev[j];r=e['residual'];n=len(r);resampled=[]
            while len(resampled)<n:
                start=int(rng.integers(0,n));resampled.extend(r[(start+np.arange(3))%n])
            noise=(e['projection']@np.array(resampled[:n])).reshape(15,15)
            replicate.append({'diff':e['diff']+noise,'err':e['err']})
        bd,be=combine(replicate);gf,_=gaussian_fit(bd,be,target,roi);boot.append([gf['x'],gf['y']])
    boot=np.array(boot);center=np.array([modes['pointing']['target_roi_gaussian_fit'][k] for k in ['x','y']])
    bootstrap={'n':len(boot),'median_xy':np.median(boot,axis=0),'interval_16_84_xy':np.percentile(boot,[16,84],axis=0),
               'covariance_xy':np.cov(boot.T),'offset_target_percentile_16_50_84':np.percentile(np.linalg.norm(boot-target,axis=1),[16,50,84]),
               'fraction_closer_to_target_than_bright_neighbor':np.mean(np.linalg.norm(boot-target,axis=1)<np.linalg.norm(boot-neighbor,axis=1)),
               'limitation':'Event plus 3-cadence synchronized residual block bootstrap. Does not cover Gaussian/PRF systematics or unresolved blending.'}
    np.savez_compressed(OUT/'pixel_photometry.npz',time_btjd=t,background=bg,pointing_dx=dx,pointing_dy=dy,reference_image=ref,
                        bootstrap_xy=boot,**{f'{k}_{label}':a for k,vals in aparrays.items() for label,a in zip(['time','flux','error'],vals)})
    row={'tic':370247818,'sector':header['SECTOR'],'camera':header['CAMERA'],'ccd':header['CCD'],'path':str(RAW),
         'archive_uri':f"https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=55.455693&dec=-79.113281&x=15&y=15&units=px&sector={header['SECTOR']}",
         'sha256':hashlib.sha256(RAW.read_bytes()).hexdigest(),'bytes':RAW.stat().st_size,'retrieved_utc':datetime.fromtimestamp(RAW.stat().st_mtime,timezone.utc).isoformat(),
         'science_array_sha256_time_quality_flux_fluxerr_ffifiles':science_digest.hexdigest(),
         'science_hash_recipe':'SHA256 of each ASCII column name followed by C-order little-endian bytes: TIME float64; QUALITY int32; FLUX float32; FLUX_ERR float32; then C-order FFI_FILE strings S38 without a column-name prefix.',
         'ffi_product_identifiers':str((OUT/'ffi_product_ids.txt').relative_to(ROOT)),
         'processing_version':header['PROCVER'],'data_release':header['DATA_REL'],'time_system':header['TIMESYS'],'bjd_reference':header['BJDREFI']+header['BJDREFF'],
         'ffi_type':header['FFI_TYPE'],'fits_checksum_verification':checksum,
         'retrieval_note':'TESSCut public FFI cutout; API regenerates FITS CHECKSUM timestamps so future retrieval may differ byte-for-byte; verify science arrays and recorded FFI identifiers.'}
    manifest=ROOT/'data/broad_pixel_manifest.json'
    previous=json.loads(manifest.read_text()) if manifest.exists() else {'products':[]}
    previous['products']=[r for r in previous['products'] if r.get('path')!=str(RAW)]+[row]
    manifest.write_text(json.dumps(safe(previous),indent=2)+'\n')
    summary={'checked_utc':datetime.now(timezone.utc).isoformat(),'ephemeris':{'period_days':P,'epoch_btjd':T0,'duration_days':DUR,'ingress_fraction':INGRESS,'cadence_days':CADENCE_DAYS},
             'input':row,'cadences':counts,'target_xy':target,'bright_neighbor_xy':neighbor,'third_star_xy':third,
             'background_pixels':int(blank.sum()),'pointing_range_pixels':{'dx':np.percentile(dx,[1,50,99]),'dy':np.percentile(dy,[1,50,99])},
             'method':'Per-event joint integrated-trapezoid + quadratic time + measured sky + reference-gradient pointing; simultaneous nuisance and transit regression. All QUALITY != 0 removed.',
             'difference_modes':modes,'localization_bootstrap':bootstrap,'apertures':results,'independent_aperture_period_checks':period_checks,
             'reference_vs_difference_roi_sensitivity':sensitivity,
             'localization_caution':'Absolute Gaussian vs WCS positions depend on PSF asymmetry and ROI/background; reference-source comparison is essential. Target-neighborhood source attribution cannot rule out unresolved companions.',
             'limitations':['Gaussian is not calibrated TESS PRF; WCS/catalog epochs not proper-motion propagated.',
                            'Gaia cone may miss close unresolved companions; source localization does not validate planet.',
                            'Short per-event windows and regressors permit correlated residual uncertainty; no survey-wide false-alarm statistic.',
                            'Diluted aperture depths are not stellar transit depths.']}
    (OUT/'pixel_vetting.json').write_text(json.dumps(safe(summary),indent=2)+'\n')
    pd.DataFrame([{'aperture':name,**event} for name,value in results.items() for event in value['events']]).to_csv(OUT/'per_event_aperture_depths.csv',index=False)
    diff=np.load(OUT/'difference_pointing.npz')['difference'];de=np.load(OUT/'difference_pointing.npz')['error']
    fig,axes=plt.subplots(1,3,figsize=(16,5))
    images=[np.log10(np.maximum(ref-np.median(ref[blank]),1)),diff,diff/de]
    titles=['Median image (log10 sky-subtracted e-/s)','Transit difference (e-/s)','Difference / diagonal error']
    for ax,im,title in zip(axes,images,titles):
        limit=np.percentile(abs(im),98) if title!=titles[0] else None
        h=ax.imshow(im,origin='lower',cmap='viridis' if limit is None else 'RdBu_r',vmin=None if limit is None else -limit,vmax=limit)
        for label,p in list(pos.items())[:3]:ax.plot(*p,'x',color='lime' if label=='target' else 'orange');ax.text(p[0]+.15,p[1]+.25,label,fontsize=8,color='white')
        ax.plot(*center,'+',color='magenta',ms=15);ax.set_title(title);fig.colorbar(h,ax=ax,shrink=.7)
    fig.tight_layout();fig.savefig(OUT/'source_localization.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(int(np.ceil((len(ev)+1)/4)),4,figsize=(16,4*int(np.ceil((len(ev)+1)/4))));lim=np.percentile(abs(np.array([e['diff'] for e in ev])),99)
    for ax,e in zip(axes.ravel(),ev):
        ax.imshow(e['diff'],origin='lower',cmap='RdBu_r',vmin=-lim,vmax=lim);ax.plot(*target,'gx');ax.plot(*neighbor,'yx');ax.set_title(f"Cycle {e['cycle']}: {e['center_btjd']:.4f}")
    axes.ravel()[-1].scatter(boot[:,0],boot[:,1],s=3,alpha=.15);axes.ravel()[-1].plot(*target,'gx');axes.ravel()[-1].plot(*neighbor,'yx');axes.ravel()[-1].set_title('Bootstrap locations')
    for ax in axes.ravel()[len(ev):-1]:ax.axis('off')
    fig.tight_layout();fig.savefig(OUT/'per_event_difference.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(9,10))
    for ax,label in zip(axes,['target_r1.5','bright_neighbor_r1.5','intermediate_r1.5']):
        at,af,ae=aparrays[label];ph=(at-T0+P/2)%P-P/2
        ax.errorbar(ph*24,af,yerr=ae,fmt='.',alpha=.5);ax.set_xlim(-5,5);ax.set_title(f"{label}: depth={results[label]['depth']*1e3:.2f} ppt, SNR={results[label]['snr']:.1f}");ax.axvline(0,color='grey',lw=.5)
    axes[-1].set_xlabel('Hours from event');fig.tight_layout();fig.savefig(OUT/'independent_apertures.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,9))
    for ax,label in zip(axes,['target','bright_neighbor','third_star']):
        for mode in ['all_quality_zero','comparison_star_image_outlier_mask']:
            a=np.load(OUT/f'aperture_bls_{label}_{mode}.npz');ax.plot(a['period'],a['power'],lw=.8,label=mode)
        ax.axvline(P,color='grey',ls='--');ax.legend(fontsize=8);ax.set_ylabel(f'{label} BLS SNR (diagonal)')
    axes[-1].set_xlabel('Trial period (days)');fig.tight_layout();fig.savefig(OUT/'regional_aperture_period_checks.png',dpi=160);plt.close(fig)
    print(json.dumps(safe({'difference_modes':modes,'bootstrap':bootstrap,'aperture_results':{k:{j:v[j] for j in ['depth','error_max_diagonal_event_scatter','snr']} for k,v in results.items()}}),indent=2))


if __name__=='__main__':main()
