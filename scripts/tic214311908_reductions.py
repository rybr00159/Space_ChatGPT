"""Fixed discovery timing in independent public reductions of the same photons.

TARS lacks reported uncertainties/quality columns, so its descriptive local
scores explicitly use adjacent-difference empirical errors and residual floors.
These are not independent observations or calibrated significance estimates.
"""
from pathlib import Path
import os,sys,json
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from astropy.io import fits
import multisector_followup as f
from archive import verify_fits_checksums
OUT=ROOT/'candidates/tic214311908/independent_reductions'

def main():
 OUT.mkdir(parents=True,exist_ok=True);pre=json.loads((ROOT/'candidates/tic214311908/holdout/preregistration.json').read_text());records=json.loads((ROOT/'data/tic214311908/tars_recognition/independent_reduction_manifest.json').read_text());result=[]
 for r in records:
  path=Path(r['path'])
  if f.digest(path)!=r['sha256']:raise ValueError('Independent-reduction original hash differs')
  with fits.open(path,memmap=False)as h:
   checks=verify_fits_checksums(path,h);d=h[1].data;head=h[1].header;h0=h[0].header;ref=head.get('BJDREFI',h0.get('BJDREFI',0))+head.get('BJDREFF',h0.get('BJDREFF',0))
   if ref!=2457000 or head.get('TIMESYS',h0.get('TIMESYS'))!='TDB':raise ValueError('Independent timing not explicit BTJD TDB')
   t=np.asarray(d['TIME'],float);is_tars=r['provenance_name']=='TARS';columns=['FLUX']if is_tars else['RAW_FLUX','CORR_FLUX','PCA_FLUX'];good=np.isfinite(t)&(np.asarray(d['QUALITY'])==0 if'QUALITY'in d.names else True)
   for column in columns:
    flux=np.asarray(d[column],float);use=good&np.isfinite(flux)&(flux>0);tt=t[use];y=flux[use];norm=np.median(y);y=y/norm
    if is_tars:
     dt=np.diff(tt);diff=np.diff(y)[(dt>0)&(dt<2.5*np.median(dt[dt>0]))];sigma=1.4826*np.median(abs(diff-np.median(diff)))/np.sqrt(2);e=np.full(len(y),sigma);error_policy='Explicit constant adjacent-difference MAD; no original QUALITY/FLUX_ERR supplied. Residual floors do not calibrate red noise.'
    else:e=np.asarray(d['FLUX_ERR'],float)[use]/norm;error_policy='Producer FLUX_ERR with original QUALITY==0; approximate errors shared by RAW/CORR/PCA flux, local residual floors retained.'
    positive=np.isfinite(e)&(e>0);tt,y,e=tt[positive],y[positive],e[positive]
    stats,local=f.local_quadratic(tt,y,e,pre['period_days'],pre['epoch_btjd'],pre['duration_days']);result.append(dict(source=r,flux_column=column,checksums=checks,error_policy=error_policy,primary=stats))
    print('REDUCTION',r['provenance_name'],r['sector'],column,{k:stats.get(k)for k in ['n_events','n_positive_events','depth','depth_err','snr']},flush=True)
 f.save(OUT/'assessment.json',dict(script_sha256=f.digest(__file__),fixed_ephemeris=pre,tests=result,limitations=__doc__))

if __name__=='__main__':main()
