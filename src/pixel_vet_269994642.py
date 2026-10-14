"""TIC269994642 independent FFI source/recurrence tests in sectors1/27/92.

Reusable fitting primitives come from pixel_vet_49617264; its files/constants
are never changed on disk. Long-baseline sector timing is independently
searched rather than assuming an accurately extrapolated 2018 epoch.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,os,urllib.request,urllib.parse,zipfile,io
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pixel_vet_49617264 as fitting
from archive import ROOT,verify_fits_checksums

TIC=269994642
P=1.7464945200053341
T0=1326.3488754694151
DURATION=.0719571950385378
INGRESS=.41423130172147016
RA,DEC=311.461584,-24.649866
OUT=ROOT/'candidates/broad_s1/tic269994642_pixels'
CACHE=Path('/workspace/tess-research-data/broad/pixels/tic269994642')


def template(t,p,t0,dt,mode):
    offset=np.linspace(-dt/2,dt/2,31)
    ph=(t[:,None]+offset-t0+p/2)%p-p/2
    if mode=='box_core':return np.mean(abs(ph)<.035/2,axis=1)
    if mode=='box_broad':return np.mean(abs(ph)<DURATION/2,axis=1)
    return np.mean(np.clip((DURATION/2-abs(ph))/(DURATION*INGRESS),0,1),axis=1)


def event_fit(t,f,e,bg,dx,dy,p,t0,dt,mode):
    events=[];sh=template(t,p,t0,dt,mode)
    for cy in range(int(np.ceil((t.min()-t0)/p)),int(np.floor((t.max()-t0)/p))+1):
        center=t0+cy*p;keep=abs(t-center)<.25;tx=t[keep]-center;s=sh[keep]
        if np.sum(s>.05)<2 or np.sum((s==0)&(tx<0))<6 or np.sum((s==0)&(tx>0))<6:continue
        cols=[np.ones(len(tx)),tx,tx**2]
        if mode!='quadratic':
            for cov in [bg[keep],dx[keep],dy[keep]]:
                cov=cov-np.mean(cov)
                if np.std(cov)>1e-10:cols.append(cov/np.std(cov))
        a=np.column_stack(cols+[-s]);inv=np.linalg.pinv(a);ff=f[keep].reshape(len(tx),-1);coef=inv@ff;res=ff-a@coef
        sigma=np.maximum(np.std(res,axis=0,ddof=a.shape[1]),np.sqrt(np.mean(e[keep].reshape(len(tx),-1)**2,axis=0)))
        ee=sigma*np.sqrt(np.sum(inv[-1]**2))
        events.append({'cycle':cy,'center_btjd':center,'n_in':int((s>.05).sum()),'n_oot':int((s==0).sum()),
                       'difference':coef[-1].reshape(f.shape[1:]),'error':ee.reshape(f.shape[1:]),'residual':res,'projection':inv[-1],
                       'time_btjd':t[keep],'model_condition':np.linalg.cond(a),'a':a,'coef':coef,'s':s,'flux':ff,'shape':f.shape[1:]})
    return events


def period_consistency(t,values,errors,dt,pseed):
    flat=np.zeros_like(values);trend=np.zeros_like(values)
    boundaries=np.r_[0,np.flatnonzero(np.diff(t)>max(.10,dt*5))+1,len(t)]
    width=max(3,int(.5/dt)//2*2+1)
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        n=min(width,(hi-lo)//2*2-1);trend[lo:hi]=median_filter(values[lo:hi],size=max(3,n),mode='nearest');flat[lo:hi]=values[lo:hi]/trend[lo:hi]
    ef=errors/trend
    noise=float(1.4826*np.median(abs(np.diff(flat)-np.median(np.diff(flat))))/np.sqrt(2));ef=np.maximum(ef,noise)
    finite=np.isfinite(flat)&np.isfinite(ef)&(trend>0)
    # Clipping only upward excursions protects unrecognized downward events.
    finite&=flat<np.median(flat[finite])+8*noise
    model=BoxLeastSquares(t[finite],flat[finite],ef[finite]);periods=np.linspace(pseed-.02,pseed+.02,4001)
    result=model.power(periods,[.035,.05,.07,.09],oversample=10,objective='snr');j=int(np.argmax(result.power))
    peak={'period_days':float(result.period[j]),'epoch_btjd':float(result.transit_time[j]),'duration_days':float(result.duration[j]),
          'depth':float(result.depth[j]),'snr_diagonal':float(result.depth_snr[j]),'error_diagonal':float(result.depth_err[j]),
          'n_used':int(finite.sum()),'noise_floor':noise,'search_period_range_days':[float(periods[0]),float(periods[-1])],
          'limitation':'Targeted all-phase recurrence check; approximate noise, no search-wide false-alarm or discovery probability.'}
    return peak,{'time_btjd':t,'flat_flux':flat,'relative_error':ef,'trend':trend,'keep':finite,'period_grid':result.period,'power':result.power}


def image_shifts(flux,reference,target):
    yy,xx=np.indices(reference.shape);use=(xx-target[0])**2+(yy-target[1])**2>2.5**2
    gy,gx=np.gradient(reference);a=np.column_stack([reference[use],gx[use],gy[use],np.ones(use.sum())]);coef=np.linalg.lstsq(a,flux[:,use].T,rcond=None)[0]
    return -coef[1]/coef[0],-coef[2]/coef[0]


def analyze(sector,gaia,nboot):
    out=OUT/f'sector{sector}';out.mkdir(parents=True,exist_ok=True)
    paths=list(CACHE.glob(f'tess-s{sector:04d}-*.fits'))
    if not paths:
        CACHE.mkdir(parents=True,exist_ok=True)
        url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urllib.parse.urlencode({'ra':RA,'dec':DEC,'x':15,'y':15,'units':'px','sector':sector})
        with urllib.request.urlopen(url,timeout=90) as response:archive=response.read(200000001)
        if len(archive)>200000000:raise RuntimeError('Cutout archive exceeds 200 MB limit.')
        with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
            for name in zipped.namelist():
                if name.endswith('.fits'):
                    dest=CACHE/Path(name).name
                    if dest.exists():raise RuntimeError('Refusing to overwrite existing raw cutout.')
                    dest.write_bytes(zipped.read(name))
        paths=list(CACHE.glob(f'tess-s{sector:04d}-*.fits'))
    if len(paths)!=1:raise RuntimeError(f'Expected one retained sector{sector}cutout,found{paths}')
    raw=paths[0];sha=hashlib.sha256(raw.read_bytes()).hexdigest();manifest=ROOT/'data/broad_pixel_269994642.json'
    records=json.loads(manifest.read_text()) if manifest.exists() else []
    expected=[r for r in records if r['sector']==sector]
    byte_hash_changed=bool(expected and expected[0]['sha256']!=sha)
    with fits.open(raw) as h:
        integrity=verify_fits_checksums(raw,h);header=dict(h[0].header);ph=dict(h[1].header);w=WCS(h[2].header)
        wcs_keywords=dict(w.to_header(relax=True))
        if header['TIMESYS']!='TDB' or header['TIMEREF']!='SOLARSYSTEM' or header['BJDREFI']!=2457000:raise ValueError('Unexpectedtimeframework')
        data=h[1].data;t=np.array(data['TIME']);quality=np.array(data['QUALITY']);flux=np.array(data['FLUX'],float);err=np.array(data['FLUX_ERR'],float)
        alltime=t.copy();ffi=np.array(data['FFI_FILE']).astype(str)
        keep=(quality==0)&np.isfinite(t)&np.all(np.isfinite(flux)&np.isfinite(err),axis=(1,2))
        counts=dict(zip(*np.unique(quality,return_counts=True)));t,flux,err=t[keep],flux[keep],err[keep]
        science=hashlib.sha256()
        for key,dtype in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4'),('FLUX_ERR','<f4')]:science.update(key.encode());science.update(np.asarray(data[key],dtype=dtype).tobytes())
        science.update(np.asarray(ffi,dtype='S38').tobytes())
    if byte_hash_changed and expected[0].get('science_array_sha256')!=science.hexdigest():
        raise ValueError('Retain changed science arrays for diagnosis; no silent expected-manifest replacement.')
    if byte_hash_changed and expected[0].get('wcs_fits_keywords')!=wcs_keywords:
        raise ValueError('Regenerated WCS differs from archived expectation; preserve for source-localization diagnosis.')
    (out/'ffi_product_identifiers.txt').write_text('\n'.join(ffi)+'\n')
    target=w.all_world2pix([[RA,DEC]],0)[0];yy,xx=np.indices(flux.shape[1:]);reference=np.median(flux,axis=0)
    physical={'rawx_min':ph['1CRV5P'],'rawx_max':ph['1CRV5P']+flux.shape[2]-1,'rawy_min':ph['2CRV5P'],'rawy_max':ph['2CRV5P']+flux.shape[1]-1,
              'target_rawx':ph['1CRV5P']+(target[0]+1-ph['1CRP5P'])*ph['1CDL5P'],
              'target_rawy':ph['2CRV5P']+(target[1]+1-ph['2CRP5P'])*ph['2CDL5P']}
    on_science_silicon=44<=physical['target_rawx']<2092 and 0<=physical['target_rawy']<2048
    if not on_science_silicon:
        row={'tic':TIC,'sector':sector,'camera':header['CAMERA'],'ccd':header['CCD'],'raw_path':str(raw),'sha256':sha,'raw_bytes':raw.stat().st_size,
             'science_array_sha256':science.hexdigest(),'archive_url':f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={RA}&dec={DEC}&x=15&y=15&units=px&sector={sector}',
             'data_release':header['DATA_REL'],'processing_version':header['PROCVER'],'cube_date':header['DATE'],'fits_checksums':integrity,
             'time_system':header['TIMESYS'],'time_reference':header['TIMEREF'],'bjdref':header['BJDREFI']+header['BJDREFF'],'cadence_days':ph['TIMEDEL'],
             'date_obs':header['DATE-OBS'],'date_end':header['DATE-END'],'time_min_btjd':np.min(alltime),'time_max_btjd':np.max(alltime),
             'ffi_identifiers':str((out/'ffi_product_identifiers.txt').relative_to(ROOT)),'physical_pixel_coordinates':physical,
             'on_science_silicon':False,'wcs_fits_keywords':wcs_keywords,'exclusion_reason':'TESSCut coverage includes collateral/prescan geometry; projected target is outside illuminated science columns.',
             'checked_utc':datetime.now(timezone.utc).isoformat()}
        records=[r for r in records if r['sector']!=sector]+[row];manifest.write_text(json.dumps(fitting.safe(records),indent=2)+'\n')
        result={'tic':TIC,'sector':sector,'raw':row,'target_xy':target,'status':'unusable_collateral_pixels_no_transit_search',
                'exclusion_reason':row['exclusion_reason'],'physical_pixel_coordinates':physical,'median_frame_flux_e_s':float(np.median(reference)),
                'warning':'Historical pre-guard search outputs are retained under exploratory_unusable_overscan and must not enter candidate tables or ephemeris fitting.'}
        (out/'pixel_vetting.json').write_text(json.dumps(fitting.safe(result),indent=2)+'\n')
        np.savez_compressed(out/'unusable_collateral_reference.npz',reference_image=reference,target_xy=target)
        return result
    stars=gaia.copy() if gaia is not None else pd.DataFrame(columns=['ra','dec','phot_g_mean_mag'])
    if len(stars):stars[['x','y']]=w.all_world2pix(stars[['ra','dec']].to_numpy(),0)
    blank=reference<np.percentile(reference,30)
    bg=np.median(flux[:,blank],axis=1);flux=flux-bg[:,None,None]
    dx,dy=image_shifts(flux,reference-np.median(reference[blank]),target)
    ap=(xx-target[0])**2+(yy-target[1])**2<2**2
    values=flux[:,ap].sum(axis=1)
    skyerr=np.sqrt(np.pi/2*np.median(err[:,blank]**2,axis=1)/blank.sum())
    errors=np.sqrt(np.sum(err[:,ap]**2,axis=1)+ap.sum()**2*skyerr**2)
    peak,curves=period_consistency(t,values,errors,float(ph['TIMEDEL']),P)
    curves.update({'raw_aperture_flux_e_s':values,'raw_aperture_flux_error_e_s':errors,'target_aperture_mask':ap,'sky_e_s_pixel':bg,'sky_error_e_s_pixel':skyerr})
    np.savez_compressed(out/'independent_aperture_periodogram.npz',**curves)
    p,t0=(P,T0) if sector==1 else (peak['period_days'],peak['epoch_btjd'])
    modes={};evs={}
    for mode in ['quadratic','pointing','box_broad','box_core']:
        events=event_fit(t,flux,err,bg,dx,dy,p,t0,float(ph['TIMEDEL']),mode);evs[mode]=events;dd,de=fitting.combine(events)
        sy,sx=np.unravel_index(np.argmax(dd/de),dd.shape);gf,_=fitting.gaussfit(dd,de,[sx,sy]);gf['offset_target_pixels']=np.hypot(gf['x']-target[0],gf['y']-target[1]);gf['sky_ra_dec']=w.all_pix2world([[gf['x'],gf['y']]],0)[0]
        modes[mode]={'n_events':len(events),'events':[{k:e[k] for k in ['cycle','center_btjd','n_in','n_oot','model_condition']} for e in events],'source_gaussian':gf}
        np.savez_compressed(out/f'difference_{mode}.npz',difference=dd,error=de,event_differences=np.array([e['difference'] for e in events]),event_errors=np.array([e['error'] for e in events]))
    gf=modes['pointing']['source_gaussian'];events=evs['pointing'];boot=fitting.boot_centroids(events,[gf['x'],gf['y']],nboot)
    positions={'target':target}
    if len(stars):
        stars['difference_distance_pixels']=np.hypot(stars.x-gf['x'],stars.y-gf['y']);stars['target_distance_pixels']=np.hypot(stars.x-target[0],stars.y-target[1]);stars.to_csv(out/'gaia_positions.csv',index=False)
        for _,s in stars[stars.phot_g_mean_mag<17].iterrows():
            if 0<=s.x<15 and 0<=s.y<15 and np.hypot(s.x-target[0],s.y-target[1])>.5:positions[f'gaia_{s.source_id}']=np.array([s.x,s.y])
        stars['bootstrap_fraction_closer_TIC_than_source']=np.mean(np.linalg.norm(boot[:,None,:]-target,axis=2)<np.linalg.norm(boot[:,None,:]-stars[['x','y']].to_numpy()[None,:,:],axis=2),axis=0)
        stars.loc[stars.source_id=='6806240691223695744','bootstrap_fraction_closer_TIC_than_source']=np.nan
        stars.to_csv(out/'gaia_positions.csv',index=False)
    difference=np.load(out/'difference_pointing.npz');source_tests={}
    for name,position in positions.items():
        hypothesis,_=fitting.gaussfit(difference['difference'],difference['error'],position,fixed_center=True)
        hypothesis['delta_chi2_vs_free_source']=hypothesis['chi2_diagonal']-gf['chi2_diagonal'];source_tests[name]=hypothesis
    sensitivity=[]
    for radius in [2.,2.5,3.]:
        roi=(xx-target[0])**2+(yy-target[1])**2<radius**2
        rg,_=fitting.gaussfit(reference,np.ones(reference.shape),target,roi)
        dg,_=fitting.gaussfit(difference['difference'],difference['error'],target,roi)
        sensitivity.append({'radius_pixels':radius,'reference':rg,'difference':dg,'difference_reference_offset_pixels':np.hypot(rg['x']-dg['x'],rg['y']-dg['y'])})
    photo,arrays=fitting.photometry(events,positions,reference,ap)
    # The helper's reference mask label is specific to native SPOC files;
    # this FFI extraction uses a declared circular aperture instead.
    for entry in photo:
        if entry['aperture']=='spoc_optimal':entry['aperture']='target_r2_reference'
    arrays['target_r2_reference']=arrays.pop('spoc_optimal')
    np.savez_compressed(out/'photometry_arrays.npz',time_btjd=t,reference=reference,blank_sky=blank,background=bg,pointing_dx=dx,pointing_dy=dy,bootstrap_xy=boot,
                        **{f'{key}_{label}':val for key,vs in arrays.items() for label,val in zip(['time_btjd','flux'],vs)})
    pd.DataFrame([{'aperture':a['aperture'],**e} for a in photo for e in a['events']]).to_csv(out/'event_aperture_depths.csv',index=False)
    row={'tic':TIC,'sector':sector,'camera':header['CAMERA'],'ccd':header['CCD'],'raw_path':str(raw),'sha256':sha,'raw_bytes':raw.stat().st_size,
         'science_array_sha256':science.hexdigest(),'sciencehash_recipe':'ASCII column name followed by C-order little-endian TIME float64, QUALITY int32, FLUX float32, FLUX_ERR float32; then S38 FFI_FILE strings.',
         'archive_url':f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={RA}&dec={DEC}&x=15&y=15&units=px&sector={sector}',
         'data_release':header['DATA_REL'],'processing_version':header['PROCVER'],'cube_date':header['DATE'],'fits_checksums':integrity,
         'time_system':header['TIMESYS'],'time_reference':header['TIMEREF'],'bjdref':header['BJDREFI']+header['BJDREFF'],'cadence_days':ph['TIMEDEL'],
         'date_obs':header['DATE-OBS'],'date_end':header['DATE-END'],'time_min_btjd':np.min(alltime),'time_max_btjd':np.max(alltime),
         'ffi_identifiers':str((out/'ffi_product_identifiers.txt').relative_to(ROOT)),'checked_utc':datetime.now(timezone.utc).isoformat()}
    row.update({'physical_pixel_coordinates':physical,'on_science_silicon':True,'wcs_fits_keywords':wcs_keywords})
    row['first_midpoint_utc']=Time(row['time_min_btjd']+row['bjdref'],format='jd',scale='tdb').utc.isot
    row['first_midpoint_minus_dateobs_seconds']=(Time(row['time_min_btjd']+row['bjdref'],format='jd',scale='tdb').utc-Time(row['date_obs'],scale='utc')).sec
    if abs(row['first_midpoint_minus_dateobs_seconds']-row['cadence_days']*86400/2)>5:raise ValueError('Cutout TIME/DATE-OBS midpoint consistency failed.')
    if byte_hash_changed:
        # Astrocut stamps a new FITS checksum creation date on regeneration.
        # Retain the original expectation and independently record the new
        # bytes only after exact scientific-array equality and checksums pass.
        row['original_expected_byte_sha256']=expected[0]['sha256'];row['regenerated_header_with_identical_science_arrays']=True
        records=[r for r in records if r['sector']!=sector]+[{**expected[0],'regenerated_byte_hashes_verified_identical_science':list(set(expected[0].get('regenerated_byte_hashes_verified_identical_science',[])+[sha]))}]
    else:records=[r for r in records if r['sector']!=sector]+[row]
    manifest.write_text(json.dumps(fitting.safe(records),indent=2)+'\n')
    result={'tic':TIC,'sector':sector,'raw':row,'target_xy':target,'cadences_total':len(alltime),'cadences_quality_zero':len(t),'quality_distribution':counts,
            'ephemeris_used':{'period_days':p,'epoch_btjd':t0,'duration_days':DURATION,'ingress_fraction':INGRESS,'source':'S1 fixed QLP template' if sector==1 else f'Independent all-phase sector{sector} pixel-aperture BLS'},
            'independent_aperture_period_check':peak,'difference_modes':modes,'apertures':photo,
            'reference_vs_difference_roi_sensitivity':sensitivity,'fixed_center_gaussian_source_hypotheses':source_tests,
            'bootstrap':{'n':nboot,'median_xy':np.median(boot,axis=0),'percentile16_84_xy':np.percentile(boot,[16,84],axis=0),
                         'offset_target_pixels_16_50_84':np.percentile(np.linalg.norm(boot-target,axis=1),[16,50,84]),'covariance':np.cov(boot.T)},
            'limitations':['Gaussian/PRF systematics, WCS, proper motion and unresolved blends are not covered by bootstrap.',
                           'Empirical pointing estimated from reference gradients outside the target may be noisy; quadratic-only baseline comparison retained.',
                           'BLS period consistency is a selected-signal recurrence test, not a search-wide FAP or confirmation.',
                           'Source localization alone cannot distinguish a planet from a stellar eclipse.']}
    (out/'pixel_vetting.json').write_text(json.dumps(fitting.safe(result),indent=2)+'\n')
    dd=np.load(out/'difference_pointing.npz');fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,image,title in zip(axes,[np.log10(np.maximum(reference-np.median(reference[blank]),1)),dd['difference'],dd['difference']/dd['error']],['Reference log10 sky-subtracted e-/s','Transit difference e-/s','Difference/diagonalerror']):
        lim=None if title.startswith('Reference') else np.percentile(abs(image),99)
        h=ax.imshow(image,origin='lower',cmap='viridis' if lim is None else 'RdBu_r',vmin=None if lim is None else -lim,vmax=lim)
        ax.plot(*target,'gx',ms=10);ax.plot(gf['x'],gf['y'],'m+',ms=13)
        if len(stars):
            for _,s in stars[stars.phot_g_mean_mag<17].iterrows():
                if 0<=s.x<15 and 0<=s.y<15:ax.plot(s.x,s.y,'o',mfc='none',mec='orange')
        ax.set_title(title);fig.colorbar(h,ax=ax,shrink=.65)
    fig.suptitle(f'TIC{TIC}sector{sector},green:TIC,magenta:difference,orangeGaiaG<17');fig.tight_layout();fig.savefig(out/'source_localization.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(9,10))
    for ax,key in zip(axes,['target_r1.0','target_r1.5','target_r2.0']):
        at,af=arrays[key];phase=(at-t0+p/2)%p-p/2;ax.plot(phase*24,af,'.',ms=2,alpha=.25);ax.set_xlim(-4,4)
        measurement=[a for a in photo if a['aperture']==key][0];ax.set_title(f"{key}:depth{measurement['depth']*1e3:.2f}ppt,SNR{measurement['snr']:.1f}")
    axes[-1].set_xlabel('Hoursfromevent');fig.tight_layout();fig.savefig(out/'target_aperture_photometry.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4));ax.plot(curves['period_grid'],curves['power']);ax.axvline(P,color='grey',ls='--');ax.set_xlabel('Period(days)');ax.set_ylabel('BLSdiagonalSNR');fig.tight_layout();fig.savefig(out/'period_consistency.png',dpi=160);plt.close(fig)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--sectors',type=int,nargs='+',default=[1,27,92]);p.add_argument('--bootstrap',type=int,default=500);p.add_argument('--gaia',type=Path);args=p.parse_args()
    gaia=pd.read_csv(args.gaia,dtype={'source_id':str}) if args.gaia else None
    for sector in args.sectors:
        r=analyze(sector,gaia,args.bootstrap)
        if r.get('status')=='unusable_collateral_pixels_no_transit_search':print(json.dumps(fitting.safe(r),indent=2),flush=True);continue
        print(json.dumps(fitting.safe({'sector':sector,'period':r['independent_aperture_period_check'],'source':r['difference_modes']['pointing']['source_gaussian'],
                                      'bootstrap':r['bootstrap'],'photometry':[{k:a[k] for k in ['aperture','depth','error_max_diagonal_event_scatter','snr']} for a in r['apertures']]}),indent=2),flush=True)


if __name__=='__main__':main()
