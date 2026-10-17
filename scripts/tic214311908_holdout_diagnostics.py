"""Descriptive phase controls and plots after the frozen native-FFI holdout test.

Offsets/apertures share data. These finite diagnostic controls are not a formal
null ensemble, independent trial count or calibrated survey false-alarm rate.
"""
from pathlib import Path
import os,sys,json
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib.pyplot as plt
import multisector_followup as f
OUT=ROOT/'candidates/tic214311908/holdout'

def main():
 pre=json.loads((OUT/'preregistration.json').read_text());p,e0,duration=[pre[k]for k in ['period_days','epoch_btjd','duration_days']];results=[];fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
 for row,sector in enumerate([103,104]):
  path=OUT/f'sector{sector}/native_aperture_arrays.npz';z=np.load(path);t,y,e=z['t'],z['r1.5_f'],z['r1.5_e'];stats,loc=f.local_quadratic(t,y,e,p,e0,duration);phase_tests=[]
  for offset in np.arange(0,1,.1):
   q,_=f.local_quadratic(t,y,e,p,e0+offset*p,duration);phase_tests.append(dict(phase_offset=float(offset),statistics=q))
  results.append(dict(sector=sector,input=dict(path=str(path.relative_to(ROOT)),sha256=f.digest(path)),primary_aperture_radius_pixels=1.5,duration_days=duration,phase_controls=phase_tests))
  axes[row,0].plot(t,(y-1)*1000,'.',ms=1,color='.45');axes[row,0].set(xlabel='BTJD TDB',ylabel='Native aperture flux (ppt)',title=f'Sector{sector}, fixed 1.5pixel aperture')
  x=f.phase(loc['t'],p,e0)*24;fy=(loc['f']-1)*1000;er=loc['e']*1000;axes[row,1].plot(x,fy,'.',ms=1,color='.7');edges=np.linspace(-3,3,31)
  for lo,hi in zip(edges[:-1],edges[1:]):
   use=(x>=lo)&(x<hi)
   if use.sum()<2:continue
   w=1/er[use]**2;mean=np.sum(w*fy[use])/w.sum();err=max(1/np.sqrt(w.sum()),np.std(fy[use],ddof=1)/np.sqrt(use.sum()));axes[row,1].errorbar(np.mean(x[use]),mean,yerr=err,fmt='o',ms=3,color='#246a73')
  axes[row,1].set(xlim=(-3,3),xlabel='Hours from frozen discovery ephemeris',ylabel='Local native flux (ppt)',title=f'Sector{sector}: selected window depth {stats["depth"]*1e6:.0f}±{stats["depth_err"]*1e6:.0f} ppm')
 fig.suptitle('Previously uninspected S103/104 original FFIs; fixed timing, no QLP timing applied');fig.savefig(OUT/'native_holdout.png',dpi=150);plt.close(fig)
 f.save(OUT/'phase_controls.json',dict(script_sha256=f.digest(__file__),preregistration=dict(path=str((OUT/'preregistration.json').relative_to(ROOT)),sha256=f.digest(OUT/'preregistration.json')),results=results,limitations=__doc__))
 for r in results:print('PHASE_CONTROLS',r['sector'],[(round(q['phase_offset'],1),round(q['statistics']['snr'],2))for q in r['phase_controls']],flush=True)

if __name__=='__main__':main()
