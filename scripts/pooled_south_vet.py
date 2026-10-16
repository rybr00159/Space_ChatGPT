"""Fixed-ephemeris falsification of selected southern pooled-search leads.

The candidate selection precedes these tests. These are conditional diagnostics,
not calibrated false-alarm probabilities or independent planet validation.
"""
from pathlib import Path
import argparse, hashlib, json, os, sys, traceback
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import pandas as pd
import multisector_followup as f
import numpy as np

ASSIGNED = {382190230: 1.282392, 375090561: 8.136209,
            293268396: 2.596856, 349972227: 8.276215,
            149625254: 10.993242, 238182368: 14.444722}

def selection(tic, requested):
    filename = 'results/multisector_mixed_pool/signals.csv.gz' if tic==238182368 else 'results/multisector/combined_signals.csv.gz'
    table = pd.read_csv(ROOT/filename)
    subset = table[(table.tic == tic) & (table.label == 'unmasked')]
    subset = subset.assign(distance=abs(subset.period_days-requested))
    row = subset.sort_values(['distance', 'width_days']).iloc[0].to_dict()
    ephem = dict(period=float(row['period_days']), epoch_btjd=float(row['epoch_btjd']),
                 duration_days=float(row['duration_days']))
    return row, ephem

def compact(result):
    return dict(sector=result['sector'], source=result['source'],
                primary={key: value['primary'] for key, value in result['column_tests'].items()},
                warnings=result['warnings'], assessment=result['assessment'],
                integrated_trapezoid=result['integrated_trapezoid'],
                quality_context=result['quality_context'])

def explicit_followup_error_policy(record):
    """Apply the existing separately labeled missing-error estimator if eligible."""
    from multisector_screen import read_curve
    import multisector_missing_errors as m
    curve=read_curve(record['path'],record['tic'],record['sector'])
    if curve['sha256']!=record['sha256']:raise ValueError('Follow-up raw checksum mismatch')
    diagnostic=m.diagnostics(curve)
    if not diagnostic['recovery_eligible']:return record
    estimate=m.estimate(curve)
    return dict(record,error_assumption='empirical_adjacent_difference_MAD',
                empirical_sigma_flux=estimate['empirical_sigma_flux'],
                error_estimator_sha256=f.digest(ROOT/'src/multisector_missing_errors.py'),
                original_producer_error_diagnostics=diagnostic,
                empirical_error_estimate=estimate,
                error_policy_note='Explicit separate uncertainty assumption; not a repair or calibrated red-noise significance.')

def incomplete_observation(record,ephem,exception):
    a,metadata=f.read_source(record);good=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a['error'])&(a['error']>0)
    t=a['TIME'][good];p,epoch,duration=[ephem[k]for k in ['period','epoch_btjd','duration_days']]
    tests={col:f.local_quadratic(t,a[col][good],a['error'][good],p,epoch,duration)[0] for col in ['corrected','raw']}
    centers=[]
    if len(t):
        for cycle in range(int(np.floor((t.min()-epoch)/p)),int(np.ceil((t.max()-epoch)/p))+1):
            center=epoch+cycle*p
            if t.min()<=center<=t.max():centers.append(dict(cycle=cycle,center_btjd=center,n_inside=int(np.sum(abs(t-center)<duration/2)),nearest_good_cadence_days=float(np.min(abs(t-center)))))
    return dict(source=record,sector=record['sector'],source_metadata=metadata,tests=tests,n_usable=len(t),predicted_event_coverage=centers,status='frozen_full_diagnostics_incomplete',error=repr(exception),note='Local-only fallback preserves zero-event coverage. An event in a gap is not a physical non-detection; full helper exception is retained.')

def conditional_joint(tic, summary):
    """One common period optimization; these sectors have already been inspected."""
    records=[r['source'] for k in ['discovery','independent','incomplete_observations'] for r in summary.get(k,[])]
    result=f.joint_recurrence(records,summary['ephemeris'],name='pooled_south_conditional_joint')
    result['interpretation']='All included sectors have already been inspected. This targeted common-period optimization is a robustness test, not withheld-sector confirmation.'
    f.save(f.OUT/f'tic{tic}'/'pooled_south_conditional_joint.json',result)

def stable_variability(tic,summary):
    """Low-order baseline with conditioning and extrapolation diagnostics.

    A long phase exclusion can make high-order Fourier extrapolation singular.
    Every fit here reports the condition and full-orbit model amplitude; these
    fits describe an alternative, and cannot exclude an overlapping transit.
    """
    from astropy.timeseries import LombScargle
    results=[]
    ephem=summary['ephemeris'];p,epoch,duration=[ephem[k] for k in ['period','epoch_btjd','duration_days']]
    for item in summary['discovery']:
        a,_=f.read_source(item['source'])
        for col in ['corrected','raw']:
            flat=f.flatten(a,col,1.5);t,y,e=[flat[k]for k in ['t','f','e']]
            oot=abs(f.phase(t,p,epoch))>min(.75*duration,.20*p)
            freq=np.arange(.2,min(30.,.45/np.median(np.diff(t))),1/(8*np.ptp(t)))
            ls=LombScargle(t[oot],y[oot],e[oot]);power=ls.power(freq);j=int(np.argmax(power));frequency=float(freq[j])
            fits=[]
            for harmonics in [1,2]:
                theta=2*np.pi*(t-t.min())*frequency
                design=np.column_stack([np.ones(len(t))]+[fun(k*theta) for k in range(1,harmonics+1) for fun in [np.sin,np.cos]])
                weighted=design[oot]/e[oot,None];condition=float(np.linalg.cond(weighted))
                beta=np.linalg.lstsq(weighted,y[oot]/e[oot],rcond=None)[0];model=design@beta
                full_amplitude=float(np.ptp(model));valid=condition<1e5 and full_amplitude<.2
                corrected=1+y-model
                stats=f.event_stats(t,corrected,e,p,epoch,duration) if valid else None
                fits.append(dict(harmonic_count=harmonics,condition=condition,full_model_peak_to_peak=full_amplitude,valid_for_descriptive_use=valid,coefficients=beta,residual_primary=stats))
            results.append(dict(sector=item['sector'],column=col,dominant_frequency_per_day=frequency,period_days=1/frequency,power=float(power[j]),n_out_of_event_fit=int(oot.sum()),fits=fits))
    f.save(f.OUT/f'tic{tic}'/'pooled_south_conditioned_variability.json',dict(tic=tic,ephemeris=ephem,tests=results,script_sha256=f.digest(__file__),warning='Dominant frequency selected from the same data. Flexible baselines can absorb real transits. Scores are conditional diagonal estimates, not detection significance or planet exclusion.'))

def transit_pixels(tic,summary,sector=95,size=9):
    """Native-cadence local event regression and descriptive difference image."""
    import io,zipfile,requests
    from urllib.parse import urlencode
    from astropy.io import fits
    from astropy.wcs import WCS
    from archive import verify_fits_checksums
    from pixel_vet_49617264 import gaussfit
    import matplotlib.pyplot as plt
    record=next(r['source'] for r in summary['discovery'] if r['sector']==sector)
    ephem=summary['ephemeris'];p,epoch,duration=[ephem[k] for k in ['period','epoch_btjd','duration_days']]
    out=f.OUT/f'tic{tic}'/f'sector{sector}_transit_difference';out.mkdir(parents=True,exist_ok=True)
    cache=f.CACHE/'followup'/f'tic{tic}'/'pixels';cache.mkdir(parents=True,exist_ok=True)
    url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urlencode(dict(ra=record['ra'],dec=record['dec'],x=size,y=size,units='px',sector=sector))
    paths=list(cache.glob(f'tess-s{sector:04d}-*_{size}x{size}_astrocut.fits'))
    if not paths:
        response=requests.get(url,timeout=100);raw=cache/f'sector{sector}_{size}x{size}_response.zip';raw.write_bytes(response.content)
        f.save(out/'http_manifest.json',dict(url=url,status=response.status_code,path=str(raw),sha256=f.digest(raw),bytes=raw.stat().st_size,checked_utc=f.utc()));response.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            for name in archive.namelist():
                if name.endswith('.fits'):
                    path=cache/Path(name).name
                    if path.exists():raise ValueError('Refusing raw pixel overwrite')
                    path.write_bytes(archive.read(name));paths.append(path)
    if len(paths)!=1:raise ValueError('Expected one native pixel cutout')
    path=paths[0];sha=f.digest(path);manifest=out/'raw_manifest.json'
    if manifest.exists() and json.loads(manifest.read_text())['sha256']!=sha:raise ValueError('Pixel checksum mismatch')
    with fits.open(path,memmap=False) as hd:
        hd.verify('exception');checks=verify_fits_checksums(path,hd);h=hd[1].header;d=hd[1].data
        if h.get('TIMESYS')!='TDB' or h.get('BJDREFI',0)+h.get('BJDREFF',0)!=2457000:raise ValueError('Pixel timing not explicitly BTJD TDB')
        wcs=WCS(hd[2].header);target=wcs.all_world2pix([[float(record['ra']),float(record['dec'])]],0)[0]
        good=(d['QUALITY']==0)&np.isfinite(d['TIME'])&np.all(np.isfinite(d['FLUX']),axis=(1,2))
        t=np.asarray(d['TIME'][good],float);cube=np.asarray(d['FLUX'][good],float);errors=np.asarray(d['FLUX_ERR'][good],float);exposure=float(h['TIMEDEL'])
        science=hashlib.sha256()
        for key,typ in [('TIME','<f8'),('QUALITY','<i4'),('FLUX','<f4'),('FLUX_ERR','<f4')]:science.update(key.encode());science.update(np.asarray(d[key],dtype=typ).tobytes())
        ffi=[str(x)for x in d['FFI_FILE']] if 'FFI_FILE' in d.names else []
    source=dict(path=str(path),url=url,sha256=sha,science_sha256=science.hexdigest(),checksums=checks,n_native=len(d),n_quality=len(t),exposure_days=exposure,bytes=path.stat().st_size)
    f.save(manifest,source);(out/'ffi_identifiers.txt').write_text('\n'.join(ffi)+'\n')
    reference=np.median(cube,axis=0);ny,nx=reference.shape;outputs=[];arrays={'reference':reference,'target_xy':target}
    qlp,_=f.read_source(record);near=np.searchsorted(qlp['TIME'],t);near=np.clip(near,1,len(qlp['TIME'])-1)
    near=np.where(abs(t-qlp['TIME'][near])<abs(t-qlp['TIME'][near-1]),near,near-1)
    time_delta=t-qlp['TIME'][near]
    intersection=(qlp['QUALITY'][near]==0)&(abs(time_delta)<exposure/3)
    blank=reference<np.percentile(reference,25)
    sky_series=np.median(cube[:,blank],axis=1)
    sky_error=1.253*np.sqrt(np.sum(errors[:,blank]**2,axis=1))/int(blank.sum())
    for common_quality,sky_subtracted,degree in [(q,s,d) for q in [False,True] for s in [False,True] for d in [1,2]]:
        cadence_mask=intersection if common_quality else np.ones(len(t),bool)
        work_t=t[cadence_mask]
        current=(cube-sky_series[:,None,None] if sky_subtracted else cube)[cadence_mask]
        current_error=(np.sqrt(errors**2+sky_error[:,None,None]**2) if sky_subtracted else errors)[cadence_mask]
        key=('qlp_quality_' if common_quality else '')+('sky_' if sky_subtracted else '')+f'degree{degree}'
        events=[]
        for cycle in np.unique(np.rint((work_t-epoch)/p).astype(int)):
            center=epoch+cycle*p;use=abs(work_t-center)<2.5*duration;tx=work_t[use]-center
            shape=np.mean(np.clip((duration/2-abs(tx[:,None]+np.linspace(-exposure/2,exposure/2,11)))/(.2*duration),0,1),axis=1)
            if np.sum(shape>.5)<3 or np.sum((shape==0)&(tx<-duration))<5 or np.sum((shape==0)&(tx>duration))<5:continue
            design=np.column_stack([tx**k for k in range(degree+1)]+[-shape]);inverse=np.linalg.pinv(design);flux=current[use].reshape(len(tx),-1);coef=inverse@flux;resid=flux-design@coef
            sigma=np.maximum(np.std(resid,axis=0,ddof=design.shape[1]),np.sqrt(np.mean(current_error[use].reshape(len(tx),-1)**2,axis=0)))
            error=sigma*np.linalg.norm(inverse[-1]);events.append(dict(cycle=int(cycle),center_btjd=float(center),depth=coef[-1].reshape(ny,nx),error=error.reshape(ny,nx),n_inside=int(np.sum(shape>.5)),n_used=len(tx),condition=float(np.linalg.cond(design))))
        if not events:raise ValueError('No bracketed pixel events')
        depths=np.array([event['depth']for event in events]);errs=np.array([event['error']for event in events]);weights=1/errs**2
        image=np.sum(weights*depths,axis=0)/np.sum(weights,axis=0);diagonal=1/np.sqrt(np.sum(weights,axis=0));scatter=np.std(depths,axis=0,ddof=1)/np.sqrt(len(events)) if len(events)>1 else diagonal;error=np.maximum(diagonal,scatter)
        seed=np.array(np.unravel_index(np.argmax(image/error),image.shape))[::-1].astype(float)
        free,_=gaussfit(image,error,seed);fixed,_=gaussfit(image,error,target,fixed_center=True)
        center=np.array([free['x'],free['y']]);radec=wcs.all_pix2world([center],0)[0]
        result=dict(baseline_degree=degree,sky_subtracted=sky_subtracted,common_qlp_quality=common_quality,n_cadences=len(work_t),n_events=len(events),events=events,target_xy=target,free_gaussian=free,target_fixed_gaussian=fixed,free_radec=radec,offset_target_pixels=float(np.linalg.norm(center-target)),max_difference_snr=float(np.max(image/error)),fixed_vs_free_delta_chi2=float(fixed['chi2_diagonal']-free['chi2_diagonal']),localization_valid=bool(np.max(image/error)>5 and not free['width_bound_hit']))
        outputs.append(result);arrays.update({key+'_difference':image,key+'_error':error})
        fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained')
        for ax,img,title in zip(axes,[reference,image,image/error],['Native median flux','Local transit difference','Conditional S/N, scatter floor']):
            im=ax.imshow(img,origin='lower');ax.plot(*target,'r+',ms=12);ax.plot(*center,'wo',ms=5,mfc='none');ax.set_title(title,fontsize=9);fig.colorbar(im,ax=ax,shrink=.7)
        fig.suptitle(f'TIC {tic}, S{sector}, P={p:.7f} d, {key}; descriptive Gaussian localization');fig.savefig(out/(key+'.png'),dpi=150);plt.close(fig)
    f.save(out/'localization.json',dict(tic=tic,sector=sector,source=source,ephemeris=ephem,script_sha256=f.digest(__file__),gaussian_helper_sha256=f.digest(ROOT/'src/pixel_vet_49617264.py'),blank_mask_n=int(blank.sum()),cadence_intersection=dict(n_common_quality=int(intersection.sum()),matching_tolerance_seconds=float(exposure/3*86400),pixel_minus_nearest_qlp_seconds_quantiles=np.quantile(time_delta,[0,.5,1])*86400,note='Nearest native cadence within one third exposure; QLP and cutout barycentric timestamps differ by tens of seconds, retained explicitly.'),background_scheme='Median of native reference-faintest quartile pixels, optionally subtracted from every pixel. Approximate median uncertainty added diagonally; induced covariance not fully modeled.',tests=outputs,limitations=['Selected ephemeris; native pixel regressions are conditional, not a calibrated confidence region.','Diagonal errors with event-scatter floor do not model all spacecraft covariance.','Gaussian shape is not a calibrated TESS PRF; no proper-motion correction or crowding deblend.','A displaced difference does not exclude a weaker overlapping on-target transit.','A fitted position with low image S/N or a Gaussian width bound is explicitly invalid for source localization.']))
    np.savez_compressed(out/'maps.npz',**arrays)
    print('TRANSIT_PIXELS',tic,[(r['baseline_degree'],r['offset_target_pixels'],r['free_radec'])for r in outputs],flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tic',type=int)
    parser.add_argument('--independent',type=int,default=0)
    parser.add_argument('--variability',action='store_true')
    parser.add_argument('--stronger-harmonic',action='store_true')
    parser.add_argument('--joint',action='store_true');parser.add_argument('--conditioned-variability',action='store_true');parser.add_argument('--pixels',action='store_true')
    args=parser.parse_args()
    for tic,period in ASSIGNED.items():
        if args.tic and args.tic!=tic:continue
        row,ephem=selection(tic,period)
        if args.stronger_harmonic:
            period = {382190230:.854650,375090561:5.424383,349972227:16.555445}.get(tic,period)
            row,ephem=selection(tic,period)
        tag='pooled_south_stronger' if args.stronger_harmonic else 'pooled_south'
        path=f.OUT/f'tic{tic}'/(tag+'.json')
        if args.joint or args.conditioned_variability or args.pixels:
            summary=json.loads(path.read_text())
            if args.joint:conditional_joint(tic,summary)
            if args.conditioned_variability:stable_variability(tic,summary)
            if args.pixels:transit_pixels(tic,summary)
            continue
        output=dict(tic=tic,checked_utc=f.utc(),selection=row,ephemeris=ephem,
                    script_sha256=f.digest(__file__),frozen_followup_method=f.method(),
                    discovery=[],independent=[],incomplete_observations=[],errors=[],limitations=[
                        'The ephemeris comes from selected pooled discovery data; per-sector tests of those data are not independent.',
                        'An additional sector is evaluated first at the frozen pooled ephemeris, with no per-sector retuning.',
                        'Diagnostic warning thresholds can over-veto shallow known planets and are not a final disposition.'])
        records=f.records_from_db(tic)
        if tic==238182368:
            identity=json.loads((f.CACHE/'combined_mixed'/f'tic{tic}.json').read_text())['identity']
            for record in records:
                source=next(s for s in identity['sources'] if int(s['sector'])==int(record['sector']))
                if record['sha256']!=source['raw_sha256']:raise ValueError('Mixed-pool source checksum mismatch')
                if source.get('error_assumption')=='empirical_adjacent_difference_MAD':
                    record.update(error_assumption=source['error_assumption'],empirical_sigma_flux=source['empirical_sigma_flux'],error_estimator_sha256=source['estimator_sha256'])
        for record in records:
            try:
                result=f.vet_observation(record,row,fixed_ephemeris=ephem)
                output['discovery'].append(compact(result))
                if args.variability:f.stellar_variability_check(record,ephem)
            except Exception as exc:
                output['errors'].append(dict(stage='discovery',sector=record['sector'],error=repr(exc),traceback=traceback.format_exc()))
            f.save(path,output)
        if args.independent and records:
            try:
                retrieved=f.retrieve_independent(records[-1],max_sectors=args.independent,
                                                exclude_sectors=[r['sector'] for r in records])
                output['retrieval']=retrieved
                for obs in retrieved['observations']:
                    if obs.get('status')!='retained':continue
                    try:
                        obs=explicit_followup_error_policy(obs)
                        result=f.vet_observation(obs,row,fixed_ephemeris=ephem)
                        output['independent'].append(compact(result))
                    except Exception as exc:
                        output['errors'].append(dict(stage='independent',sector=obs['sector'],error=repr(exc),traceback=traceback.format_exc()))
                        try:output['incomplete_observations'].append(incomplete_observation(obs,ephem,exc))
                        except Exception as fallback:output['errors'].append(dict(stage='independent_local_fallback',sector=obs['sector'],error=repr(fallback)))
                    f.save(path,output)
            except Exception as exc:output['errors'].append(dict(stage='retrieval',error=repr(exc),traceback=traceback.format_exc()))
        f.save(path,output)
        print('SOUTH',tic,ephem,[(r['sector'],round(r['primary']['corrected']['snr'],2),round(r['primary']['raw']['snr'],2),r['warnings']) for r in output['discovery']+output['independent']], 'errors',output['errors'],flush=True)

if __name__=='__main__':main()
