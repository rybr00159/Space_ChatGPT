"""Consolidate original independent checks without double-counting photons."""
from pathlib import Path
import os,sys,json,csv
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib.pyplot as plt
import multisector_followup as f
OUT=ROOT/'candidates/tic214311908'

def main():
 discovery=json.loads((OUT/'raw_local_vetting.json').read_text());pre=json.loads((OUT/'holdout/preregistration.json').read_text());rows=[];fig,axes=plt.subplots(3,2,figsize=(12,10),layout='constrained')
 for i,r in enumerate(discovery['observations']):
  a,m=f.read_source(r['source']);good=(a['QUALITY']==0)&np.isfinite(a['TIME'])&np.isfinite(a['error'])&(a['error']>0);stats,loc=f.local_quadratic(a['TIME'][good],a['raw'][good],a['error'][good],pre['period_days'],pre['epoch_btjd'],pre['duration_days']);ax=axes.flat[i]
  rows.append(dict(sector=r['sector'],product='QLP',flux_normalization='TIC target-reference flux; assumes varying source is target',n_local_events=stats['n_events'],n_positive=stats['n_positive_events'],local_depth_ppm=stats['depth']*1e6,local_error_ppm=stats['depth_err']*1e6,conditional_local_snr=stats['snr'],assessment_file=f'raw_local_vetting.json'))
  draw(ax,loc,pre,f'Sector{r["sector"]}: QLP raw target-reference flux')
 for i,sector in enumerate([103,104],4):
  z=np.load(OUT/f'holdout/sector{sector}/native_aperture_arrays.npz');stats,loc=f.local_quadratic(z['t'],z['r1.5_f'],z['r1.5_e'],pre['period_days'],pre['epoch_btjd'],pre['duration_days']);rows.append(dict(sector=sector,product='native_FFI_holdout',flux_normalization='Total sky-subtracted1.5pixel aperture flux; blended, unlikeQLP',n_local_events=stats['n_events'],n_positive=stats['n_positive_events'],local_depth_ppm=stats['depth']*1e6,local_error_ppm=stats['depth_err']*1e6,conditional_local_snr=stats['snr'],assessment_file=f'holdout/sector{sector}/assessment.json'));draw(axes.flat[i],loc,pre,f'Sector{sector}: native blended aperture (held out)')
 fig.suptitle('TIC214311908 field: fixed5.691498-day timing; QLP and native apertures have different flux normalization');fig.savefig(OUT/'six_sector_local_folds.png',dpi=150);plt.close(fig)
 with(OUT/'measurements.csv').open('w',newline='')as stream:w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 frozen=json.loads((OUT/'independent_variability/ephemeris_frozen_pre_holdout.json').read_text())
 # Preserve the original worker-local record while attaching the global chronology correction.
 frozen['global_chronology_audit']='holdout/timeline_audit.json'
 frozen['global_chronology_note']='Discovery-only worker had no holdout access, but saved this snapshot after root retrieval. Predictive tests used the earlier preregistered BLS timing.'
 f.save(OUT/'signal_summary.json',dict(tic=214311908,review_date='2026-10-17',script_sha256=f.digest(__file__),selected_discovery_ephemeris=pre,independent_frozen_timing=frozen,primary_checks=rows,observed_sectors=[12,39,66,93,103,104],primary_sampled_positive_event_windows=sum(r['n_positive']for r in rows),unique_stars_nominated=1,classification='promising_recurrent_transit_or_eclipsing_signal_unresolved_host',planet_established=False,novelty_established=False,openrouter_research_cost_usd=0,notes=['All21sampled primary windows are positive, but not all individually significant.','TARS/ELEANOR reductions share originalphotons and are not newsector observations.','Conditional scores omit fullrednoise/trials; no falsepositiveprobability measured.','ClosestfaintGaia/TIC1515709327at3.15arcsec and comovingTIC214311910 remain viablehostalternatives.']))

def draw(ax,loc,pre,title):
 x=f.phase(loc['t'],pre['period_days'],pre['epoch_btjd'])*24;y=(loc['f']-1)*1000;e=loc['e']*1000;ax.plot(x,y,'.',ms=1,color='.72');edges=np.linspace(-3,3,25)
 for lo,hi in zip(edges[:-1],edges[1:]):
  g=(x>=lo)&(x<hi)
  if g.sum()<2:continue
  w=1/e[g]**2;mean=np.sum(w*y[g])/w.sum();err=max(1/np.sqrt(w.sum()),np.std(y[g],ddof=1)/np.sqrt(g.sum()));ax.errorbar(np.mean(x[g]),mean,yerr=err,fmt='o',ms=3,color='#234e70')
 ax.axvline(0,lw=.6,color='.5');ax.set(xlim=(-3,3),xlabel='Hours from fixed discovery ephemeris',ylabel='Local flux (ppt)',title=title)

if __name__=='__main__':main()
