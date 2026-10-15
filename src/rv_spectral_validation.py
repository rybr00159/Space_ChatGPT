"""Reanalyze published GALAH observed-flux arrays with explicit frame handling.

The downloaded fit spectra have REST-FRAME wavelengths. Their stored fitted
velocity is read from matching fit products to invert the documented
transformation and reconstruct the pipeline's original barycentric grid. This
is a same-data line-consistency check, not independent raw wavelength calibration
or an independent observing epoch. It cannot establish a new stellar companion.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import os
import requests
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize_scalar
os.environ.setdefault('MPLCONFIGDIR','/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
RAW=Path('/workspace/tess-research-data/new_stars/rv/spectra')
META=ROOT/'data/new_stars/rv/spectra'
OUT=ROOT/'candidates/new_stars/rv'
C=299792.458
PAIRS=[('4690783166942834816','140806004701067','230803003801128'),
       ('5437472031746145536','140312002701153','230113003201154')]
REGIONS=[('blue',4720,4850),('green',5660,5860),('red1',6490,6545),('red2',6585,6740)]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def product(sid,kind,download=False):
    suffix='npz'if kind=='covariances'else'fits'
    extra={'spectrum':'','covariances':'_cov','results':'_results'}[kind]
    path=RAW/(sid+extra+'.'+suffix)
    manifest=META/(sid+extra+'.json')
    RAW.mkdir(parents=True,exist_ok=True);META.mkdir(parents=True,exist_ok=True)
    if path.exists():
        record=json.loads(manifest.read_text())
        if sha(path)!=record['sha256']:raise ValueError('Raw identity mismatch '+str(path))
    elif download:
        name=kind+'.'+suffix
        url=f'https://cloud.datacentral.org.au/teamdata/GALAH/public/GALAH_DR4/analysis_products_single/{sid[:6]}/{sid}/{sid}_single_fit_{name}'
        response=requests.get(url,timeout=90);body=response.content
        path.write_bytes(body);record=dict(url=response.url,http_status=response.status_code,path=str(path),
            retrieved_utc=datetime.now(timezone.utc).isoformat(),bytes=len(body),sha256=sha(path))
        manifest.write_text(json.dumps(record,indent=2)+'\n');response.raise_for_status()
    else:raise FileNotFoundError(str(path)+'; use --download')
    if record['http_status']!=200:raise ValueError('Failed product response '+str(path))
    return path,record


def load(sid,download=False):
    path,record=product(sid,'spectrum',download)
    covpath,covrecord=product(sid,'covariances',download)
    with np.load(covpath,allow_pickle=False)as data:
        labels=data['model_labels'].astype(str);parameters=data['model_parameters']
        where=np.flatnonzero(labels=='vrad')
        velocity=float(parameters[where[0]])if len(where)==1 else None
    fitrecord=None;velocity_provenance='covariance_vrad'
    if velocity is None:
        # In the released single-star fit the RV is commonly fixed at the
        # initial Gaussian fit and is therefore absent from the covariance.
        # The pipeline writes that applied value as rv_gauss in fit_results.
        fitpath,fitrecord=product(sid,'results',download)
        with fits.open(fitpath,checksum=True)as hd:
            if 'rv_gauss'not in hd[1].columns.names:raise ValueError('No documented applied velocity')
            velocity=float(hd[1].data['rv_gauss'][0])
            velocity_provenance='rv_gauss_from_matching_results; fixed-RV transformation assumed as in published analysis source'
    with fits.open(path,checksum=True)as hd:
        values=hd[1].data
        wave=np.asarray(values['wave'],float).copy()
        flux=np.asarray(values['sob'],float).copy();error=np.asarray(values['uob'],float).copy()
        mask=np.asarray(values['mob'],bool).copy()
    valid=np.isfinite(wave)&np.isfinite(flux)&np.isfinite(error)&(error>0)
    return dict(wave=wave[valid]*(1+velocity/C),rest_wave=wave[valid],flux=flux[valid],
                error=error[valid],mask=mask[valid],frame_velocity=velocity,
                products=[v for v in [record,covrecord,fitrecord]if v is not None],id=sid,
                velocity_provenance=velocity_provenance)


def profile_shift(first,second,lo,hi,masked=False):
    a=(first['rest_wave']>lo)&(first['rest_wave']<hi)
    b=(second['rest_wave']>lo+1)&(second['rest_wave']<hi-1)
    if masked:a&=first['mask'];b&=second['mask']
    w1=first['wave'][a];f1=first['flux'][a];e1=first['error'][a]
    w2=second['wave'][b];f2=second['flux'][b];e2=second['error'][b]
    if min(len(w1),len(w2))<50:raise ValueError('Insufficient spectral coverage')
    x=(w2-np.median(w2))/np.ptp(w2)
    def evaluate(v,details=False):
        shifted=w2/(1+v/C)
        comparison=np.interp(shifted,w1,f1)
        error=np.hypot(e2,np.interp(shifted,w1,e1))
        design=np.column_stack([np.ones(len(w2)),x,comparison])
        weight=design/error[:,None];target=f2/error
        coefficient=np.linalg.lstsq(weight,target,rcond=None)[0]
        residual=target-weight@coefficient
        statistic=float(residual@residual)
        if coefficient[-1]<0:statistic+=1e9
        return (statistic,coefficient,comparison)if details else statistic
    grid=np.arange(-25,25.00001,.1);stat=np.array([evaluate(v)for v in grid]);j=int(np.argmin(stat))
    if j in [0,len(grid)-1]:raise ValueError('Shift optimum hits search boundary')
    fit=minimize_scalar(evaluate,bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-5})
    best=float(fit.x);step=.03
    curvature=(evaluate(best+step)+evaluate(best-step)-2*evaluate(best))/step**2
    err=float(np.sqrt(2/curvature))if curvature>0 else None
    statistic,beta,comparison=evaluate(best,True)
    return dict(shift_kms=best,formal_profile_error_kms=err,chi2=statistic,n_points=len(w2),
                continuum_and_depth_coefficients=beta.tolist(),masked=masked),dict(grid=grid,chi2=stat,
                    wavelength=w2,observed=f2,old_on_new_grid=comparison)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--download',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);results=[]
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for index,(source,old,new)in enumerate(PAIRS):
        first=load(old,args.download);second=load(new,args.download);measurements=[]
        for name,lo,hi in REGIONS:
            for masked in [False,True]:
                result,arrays=profile_shift(first,second,lo,hi,masked)
                result['region']=name;measurements.append(result)
                if not masked:
                    axes[index,0].plot(arrays['grid'],arrays['chi2']-arrays['chi2'].min(),label=name)
                    axes[index,1].errorbar(len(measurements)//2,result['shift_kms'],yerr=result['formal_profile_error_kms'],fmt='o')
        delta=second['frame_velocity']-first['frame_velocity']
        axes[index,0].axvline(delta,color='black',ls='--',label='Catalogue-frame difference')
        axes[index,0].set(xlabel='New minus old velocity (km/s)',ylabel='Conditional delta chi-square',ylim=(0,60),title=source)
        axes[index,0].legend(fontsize=8)
        axes[index,1].axhline(delta,color='black',ls='--');axes[index,1].set(xlabel='Independent wavelength region index',ylabel='Same-data line shift (km/s)')
        # A shifted copy of an actual observed spectrum validates velocity
        # sign and units on its original sampling; it does not test reality.
        synthetic=dict(first);synthetic['wave']=first['wave']*(1+8/C)
        control,_=profile_shift(first,synthetic,5660,5860)
        assert abs(control['shift_kms']-8)<.02
        results.append(dict(source_id=source,old_id=old,new_id=new,old_frame_velocity=first['frame_velocity'],
            new_frame_velocity=second['frame_velocity'],catalogue_frame_difference=delta,
            velocity_provenance=[first['velocity_provenance'],second['velocity_provenance']],
            line_region_checks=measurements,synthetic_8kms_control=control,
            raw_products=first['products']+second['products']))
        print(source,'catalogue-frame delta',delta,'regions',[(m['region'],m['masked'],m['shift_kms'])for m in measurements],flush=True)
    result=dict(checked_utc=datetime.now(timezone.utc).isoformat(),method_sha256=sha(Path(__file__)),targets=results,
        transformation='wave_original_barycentric = published_rest_wave * (1 + exact_fit_vrad/299792.458)',
        documented_source='GALAH_DR4 fd661a419681e51dfe2d241465aa24ae97c57ad3 spectrum_analysis/galah_dr4_spectrum_analysis_single.py:1096,1364,2110',
        limitations=['Published normalized observed-flux arrays are reused; these are not independent observing epochs.',
            'Rest-frame wavelengths incorporate the original fitted velocity. Inversion uses that value and is not independent calibration proof.',
            'Regional reanalysis can reveal line inconsistency or masking dependence; agreement does not prove a companion or its mass.',
            'Formal profile errors omit continuum, interpolation, calibration and stellar variability systematics.'])
    (OUT/'spectral_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    fig.tight_layout();fig.savefig(OUT/'spectral_validation.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
