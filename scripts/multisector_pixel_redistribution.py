"""Conditional empirical flux-loss versus image-redistribution diagnostics.

Reference-image gradients approximate pointing/focus artifacts and can also
absorb a blended astrophysical event. These tests do not identify a host or
provide calibrated significance. TIC229146531 is included as a positive control
for a stable target-area difference, not as an independently validated planet.
"""
from pathlib import Path
import os,sys
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import multisector_followup as f

def main():
 for tic,sector in [(183373155,96),(183373155,105),(71274048,105),(229146531,96)]:
  folder=f.OUT/f'tic{tic}'/f'sector{sector}_transit_difference';path=folder/'maps.npz';data=np.load(path);ref=data['reference'];target=data['target_xy'];yy,xx=np.indices(ref.shape);gy,gx=np.gradient(ref);basis=[ref,gx,gy,np.ones_like(ref),xx-target[0],yy-target[1]];tests=[]
  for radius in [2.5,3.5]:
   roi=(xx-target[0])**2+(yy-target[1])**2<radius**2;design=np.column_stack([v[roi]for v in basis])
   for key in data.files:
    if not key.endswith('_difference'):continue
    name=key.removesuffix('_difference');image=data[key];error=data[name+'_error'];a=design/error[roi,None];y=image[roi]/error[roi];beta=np.linalg.lstsq(a,y,rcond=None)[0];chi=float(np.sum((y-a@beta)**2));dof=len(y)-len(beta);cov=np.linalg.pinv(a.T@a)*max(1,chi/dof);se=np.sqrt(np.maximum(0,np.diag(cov)));compare={}
    for label,cols in [('no_gradients',[0,3,4,5]),('no_flux_scaling',[1,2,3,4,5])]:
     small=a[:,cols];b=np.linalg.lstsq(small,y,rcond=None)[0];compare[label]=float(np.sum((y-small@b)**2))
    tests.append(dict(mode=name,radius_pixels=radius,reference_fractional_decrement=float(beta[0]),scaled_diagonal_error=float(se[0]),conditional_snr=float(beta[0]/se[0]),gradient_coefficients_pixels=beta[1:3],chi2=chi,dof=dof,reduced_model_chi2=compare,condition=float(np.linalg.cond(a)),coefficients=beta))
  f.save(folder/'redistribution.json',dict(tic=tic,sector=sector,script_sha256=f.digest(__file__),input_map=dict(path=str(path.relative_to(ROOT)),sha256=f.digest(path)),parameters=['reference_fractional_decrement','x_gradient_pixels','y_gradient_pixels','constant_background','x_background_gradient','y_background_gradient'],tests=tests,limitations=__doc__))
  print('REDISTRIBUTION',tic,sector,[(r['mode'],r['radius_pixels'],round(r['reference_fractional_decrement']*1e6),round(r['scaled_diagonal_error']*1e6),round(r['conditional_snr'],2),round(r['chi2']/r['dof'],2))for r in tests if r['mode']in ['qlp_quality_degree2','qlp_quality_sky_degree2']],flush=True)

if __name__=='__main__':main()
