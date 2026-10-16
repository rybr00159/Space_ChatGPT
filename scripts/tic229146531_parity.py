"""Signed alternating-event tests at a single common ephemeris for TIC229146531.

All sectors have already been inspected. Sector-dependent normalizations avoid
pooling differing apparent depths as though the aperture dilution were fixed.
These conditional diagonal estimates do not calibrate survey/analysis trials.
"""
from pathlib import Path
import json,os,sys
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import multisector_followup as f
OUT=f.OUT/'tic229146531'

def main():
 summary=json.loads((OUT/'batch02_followup.json').read_text());joint=json.loads((OUT/'batch02_discovery_joint.json').read_text());ep=dict(period=joint['period_days'],epoch_btjd=joint['epoch_btjd'],duration_days=joint['duration_days']);records=[r['source']for branch in ['discovery','other_batch','additional']for r in summary[branch]]
 result=dict(tic=229146531,script_sha256=f.digest(__file__),methods={'src/multisector_followup.py':f.digest(ROOT/'src/multisector_followup.py')},ephemeris=ep,tests=[],limitations=__doc__)
 for duration in [ep['duration_days'],.05]:
  rows=[]
  for record in records:
   a,meta=f.read_source(record);good=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a['error'])&(a['error']>0)
   for col in ['raw','corrected']:
    t,y,e=a['TIME'][good],a[col][good],a['error'][good];p,epoch=ep['period'],ep['epoch_btjd'];depths=[]
    for shift in [0,p]:
     stats,_=f.local_quadratic(t,y,e,2*p,epoch+shift,duration);depths.append(stats)
    normal=.5*(depths[0]['depth']+depths[1]['depth']);difference=depths[0]['depth']-depths[1]['depth'];error=np.hypot(depths[0]['depth_err'],depths[1]['depth_err'])
    secondary=f.local_quadratic(t,y,e,p,epoch+p/2,duration)[0]
    rows.append(dict(sector=record['sector'],source=record,column=col,first_depth=depths[0]['depth'],alternating_depth=depths[1]['depth'],first_error=depths[0]['depth_err'],alternating_error=depths[1]['depth_err'],difference=difference,difference_error=error,signed_sigma=difference/error,mean_depth=normal,fractional_difference=difference/normal,fractional_error=error/normal,secondary={k:v for k,v in secondary.items()if k!='events'}))
  combined=[]
  for col in ['raw','corrected']:
   rr=[r for r in rows if r['column']==col];x=np.array([r['fractional_difference']for r in rr]);e=np.array([r['fractional_error']for r in rr]);w=1/e**2;avg=float(np.sum(w*x)/np.sum(w));err=float(1/np.sqrt(np.sum(w)));chi=float(np.sum(((x-avg)/e)**2));scatter_scale=max(1.,np.sqrt(chi/(len(x)-1)));combined.append(dict(column=col,mean_fractional_alternating_difference=avg,diagonal_error=err,signed_diagonal_sigma=avg/err,heterogeneity_chi2=chi,heterogeneity_dof=len(x)-1,scatter_inflated_error=err*scatter_scale,signed_scatter_inflated_sigma=avg/(err*scatter_scale),warning='Mean depths in denominators have uncertainty; covariance, red noise and analysis trials are not included. This is a descriptive consistency test.'))
  result['tests'].append(dict(duration_days=duration,sectors=rows,combined=combined))
 f.save(OUT/'alternating_depths.json',result)
 for z in result['tests']:
  print('PARITY',z['duration_days'],z['combined'],flush=True)

if __name__=='__main__':main()
