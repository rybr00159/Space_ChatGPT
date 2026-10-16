"""Six-sector source-hashed raw local folds at a single discovery ephemeris."""
from pathlib import Path
import os,sys,json
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib.pyplot as plt
import multisector_followup as f
OUT=f.OUT/'tic229146531'

def main():
 s=json.loads((OUT/'batch02_followup.json').read_text());j=json.loads((OUT/'batch02_discovery_joint.json').read_text());p,epoch=j['period_days'],j['epoch_btjd'];records=sorted([r['source']for group in ['discovery','other_batch','additional']for r in s[group]],key=lambda r:r['sector']);fig,axes=plt.subplots(2,3,figsize=(12,6),layout='constrained');arrays={};rows=[]
 for ax,r in zip(axes.flat,records):
  a,m=f.read_source(r);g=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a['error'])&(a['error']>0);stats,loc=f.local_quadratic(a['TIME'][g],a['raw'][g],a['error'][g],p,epoch,.05);x=f.phase(loc['t'],p,epoch)*24;y=(loc['f']-1)*1000;e=loc['e']*1000;use=abs(x)<2
  ax.plot(x[use],y[use],'.',ms=1,color='.75');edges=np.linspace(-2,2,41)
  for lo,hi in zip(edges[:-1],edges[1:]):
   use=(x>=lo)&(x<hi)
   if use.sum()<2:continue
   w=1/e[use]**2;mean=np.sum(w*y[use])/sum(w);err=max(1/np.sqrt(sum(w)),np.std(y[use],ddof=1)/np.sqrt(use.sum()));ax.errorbar(np.mean(x[use]),mean,yerr=err,fmt='o',ms=3,color='#20639b',lw=.7)
  ax.set(xlim=(-2,2),xlabel='Hours from common ephemeris',ylabel='Raw local flux (ppt)',title=f'Sector {r["sector"]}');arrays.update({f's{r["sector"]}_{k}':loc[k]for k in ['t','f','e','cycle']});rows.append(dict(sector=r['sector'],source=r,metadata=m,local_statistics={k:v for k,v in stats.items()if k!='events'}))
 fig.suptitle('TIC 229146531: repeating target-area dips; planetary origin unresolved');fig.savefig(OUT/'six_sector_raw_folds.png',dpi=150);plt.close(fig);np.savez_compressed(OUT/'six_sector_raw_arrays.npz',**arrays);f.save(OUT/'six_sector_raw_folds.json',dict(script_sha256=f.digest(__file__),period_days=p,epoch_btjd=epoch,local_window_duration_days=.05,sectors=rows,warning='Common ephemeris was selected using the discovery sectors. Local polynomials, different apertures and correlated errors limit physical interpretation. No independent significance or novelty claim.'))

if __name__=='__main__':main()
