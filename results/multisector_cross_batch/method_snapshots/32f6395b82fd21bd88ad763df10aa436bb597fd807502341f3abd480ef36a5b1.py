"""Blind pooled searches of every repeated target, including weak single sectors.

Each sector is independently prepared before pooling; no median crosses an
inter-sector gap. Cross-sector statistics at pooled-fitted ephemerides are
conditional consistency tests, not independent confirmation or false-alarm rates.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import argparse,csv,hashlib,json,sqlite3,time

import numpy as np
from astropy.timeseries import BoxLeastSquares

from fast_screen import FastScreenConfig,period_grid,_peak_seeds,event_stats
from multisector_screen import WIDTHS,seed_reasons
from multisector_recognition import FreshCatalogIndex

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/workspace/tess-research-data/multisector')
OUT=ROOT/'results/multisector'
CONFIG=FastScreenConfig(min_period=.5,max_period=40.,durations=(.035,.05,.075,.1,.15,.225),
                        n_peaks=8,n_refine=32,grid_drift_durations=.8,compute_full_stats=True)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def search_prepared(t,f,e):
    baseline=float(np.ptp(t));periods=period_grid(baseline,CONFIG)
    model=BoxLeastSquares(t,f,e)
    coarse=model.power(periods,CONFIG.durations,oversample=5,objective='snr')
    seeds=_peak_seeds(periods,coarse.power,baseline,CONFIG.n_refine)
    step=min(CONFIG.durations)*CONFIG.grid_drift_durations/baseline
    refined_periods=np.unique(np.concatenate([periods[j]*np.exp(np.linspace(-step,step,31))for j in seeds]))
    refined_periods=refined_periods[(refined_periods>=CONFIG.min_period)&(refined_periods<=periods[-1])]
    refined=model.power(refined_periods,CONFIG.durations,oversample=12,objective='snr')
    peaks=[]
    for j in _peak_seeds(refined_periods,refined.power,baseline,CONFIG.n_peaks):
        p,epoch,duration=(float(x[j])for x in [refined.period,refined.transit_time,refined.duration])
        stats=event_stats(t,f,e,p,epoch,duration)
        guards=[]
        if stats['n_events']<2:guards.append('fewer_than_two_sampled_events')
        if stats['n_positive_events']<2:guards.append('fewer_than_two_positive_local_depths')
        if duration/p>.15:guards.append('excessive_duty_cycle')
        if (stats['event_snr_max_fraction']or 0)>.85:guards.append('single_event_dominates')
        peaks.append(dict(period=p,t0_btjd=epoch,duration_days=duration,depth=float(refined.depth[j]),
                          depth_err=float(refined.depth_err[j]),snr=float(refined.depth_snr[j]),
                          screening_pass=not guards,guard_failures=guards,**stats))
    return dict(baseline_days=baseline,grid_periods=len(periods),maximum_period_searched=float(periods[-1]),peaks=peaks)


def run(limit=0):
    with sqlite3.connect((CACHE/'survey.sqlite').as_uri()+'?mode=ro',uri=True)as conn:
        records=[json.loads(r[0])for r in conn.execute("SELECT result_json FROM observations WHERE status='screened' ORDER BY tic,sector")]
    grouped=defaultdict(list)
    for record in records:grouped[record['tic']].append(record)
    repeated={tic:rows for tic,rows in grouped.items()if len(rows)>=2}
    if limit:repeated=dict(list(repeated.items())[:limit])
    index=FreshCatalogIndex();destination=CACHE/'combined';destination.mkdir(parents=True,exist_ok=True)
    methods={n:sha(ROOT/n)for n in ['src/multisector_combine.py','src/multisector_screen.py','src/fast_screen.py']}
    summaries=[];table=[]
    for i,(tic,rows)in enumerate(sorted(repeated.items()),1):
        path=destination/f'tic{tic}.json';sources=[dict(sector=r['sector'],path=r['path'],sha256=r['sha256'])for r in rows]
        identity=dict(methods=methods,sources=sources,catalog_snapshot=index.metadata()['fresh_snapshot'])
        if path.exists():
            result=json.loads(path.read_text())
            if result['identity']!=identity:raise ValueError('Pooled checkpoint differs; preserve and fork rather than replace '+str(path))
        else:
            scans=[]
            for width in WIDTHS:
                chunks=[]
                for record in rows:
                    p=CACHE/'prepared'/f'tic{tic}_s{record["sector"]:04d}_w{width:g}.npz'
                    with np.load(p)as data:chunk={k:data[k]for k in ['t','f','e']}
                    chunk['sector']=np.full(len(chunk['t']),record['sector']);chunks.append(chunk)
                joined={k:np.concatenate([c[k]for c in chunks])for k in ['t','f','e','sector']}
                order=np.argsort(joined['t']);joined={k:v[order]for k,v in joined.items()}
                known_keep=index.transit_mask(tic,joined['t'],float(rows[0]['ra']),float(rows[0]['dec']))
                for label,mask in [('unmasked',np.ones(len(joined['t']),bool)),('known_masked_residual',known_keep)]:
                    if label!='unmasked'and mask.all():continue
                    if mask.sum()<50:continue
                    t,f,e=(joined[k][mask]for k in ['t','f','e']);ss=joined['sector'][mask]
                    scan=search_prepared(t,f,e);scan.update(width_days=width,label=label,n_known_masked=int((~mask).sum()))
                    for peak in scan['peaks']:
                        peak['per_sector']={}
                        for sector in sorted(set(ss)):
                            use=ss==sector
                            peak['per_sector'][str(sector)]=event_stats(t[use],f[use],e[use],peak['period'],peak['t0_btjd'],peak['duration_days'])
                        peak['positive_support_sectors']=sum(v['n_positive_events']>=1 and v['local_snr']>=2.5 for v in peak['per_sector'].values())
                    scans.append(scan)
            result=dict(tic=tic,sectors=[r['sector']for r in rows],identity=identity,scans=scans,
                        checked_utc=datetime.now(timezone.utc).isoformat(),
                        limitations='Pooled ephemeris optimized on all sectors; per-sector scores are conditional, not holdouts. Known masks applied to independently detrended arrays.')
            path.write_text(json.dumps(result,indent=2)+'\n')
        summaries.append(dict(tic=tic,sectors=result['sectors'],path=str(path),sha256=sha(path)))
        for scan in result['scans']:
            for rank,peak in enumerate(scan['peaks'],1):
                match=index.match_signal(tic,peak['period'],float(rows[0]['ra']),float(rows[0]['dec']))
                table.append(dict(tic=tic,sectors=';'.join(map(str,result['sectors'])),width_days=scan['width_days'],label=scan['label'],rank=rank,
                                  period_days=peak['period'],epoch_btjd=peak['t0_btjd'],duration_days=peak['duration_days'],depth_ppm=peak['depth']*1e6,
                                  local_depth_ppm=(peak['local_depth']or 0)*1e6,local_snr=peak['local_snr'],n_events=peak['n_events'],
                                  n_positive_events=peak['n_positive_events'],positive_support_sectors=peak['positive_support_sectors'],
                                  seed_retained=not seed_reasons(peak),known_period_match=match['known_period_match'],known_harmonic_match=match['known_harmonic_match'],
                                  recognition_status=match['recognition_status'],evidence_path=str(path)))
        if i%16==0:print('POOLED',i,'/',len(repeated),flush=True)
    OUT.mkdir(parents=True,exist_ok=True)
    if table:
        with (OUT/'combined_signals.csv').open('w',newline='')as stream:
            writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    report=dict(checked_utc=datetime.now(timezone.utc).isoformat(),methods=methods,repeated_stars_searched=len(summaries),
                sector_observations_in_pooled_search=sum(len(r['sectors'])for r in summaries),signal_rows=len(table),
                both_sector_retained_rows=sum(r['seed_retained']and r['positive_support_sectors']>=2 for r in table),sources=summaries,
                limitation='Selection includes all repeated searched targets, regardless of single-sector SNR. Counts overlap the primary survey.')
    (OUT/'combined_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print('POOLED COMPLETE',len(summaries),'stars;',len(table),'rows',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int,default=0);args=parser.parse_args();run(args.limit)
