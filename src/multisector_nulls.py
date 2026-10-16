"""Small moving-block surrogate diagnostic for the expanded TESS screen.

This is not a distribution-free null or survey false-alarm calibration. Fixed
pre-flux targets can contain astrophysical events. Block joins alter their time
correlations; signals longer than a block and nonstationarity are not preserved.
Original data, inputs and existing control experiments are never modified.
"""
from pathlib import Path
from datetime import datetime,timezone
import csv
import hashlib
import json
import sys
import numpy as np
from scipy.ndimage import median_filter

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from multisector_controls import read_observation
from multisector_screen import screen_arrays,seed_reasons

OUT=ROOT/'results/multisector_controls/nulls'
CACHE=Path('/workspace/tess-research-data/multisector/controls/nulls')
SEED=20261016


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def surrogate(arrays,rng,block_days=.4):
    t,f,e,q=arrays
    use=np.isfinite(t)&np.isfinite(f)&np.isfinite(e)&(f>0)&(e>0)&(q==0)
    order=np.argsort(t[use]);t,f,e=[a[use][order] for a in [t,f,e]]
    cadence=float(np.median(np.diff(t)[np.diff(t)>0]))
    gap=max(1/48,2.5*cadence)
    starts=np.r_[0,np.flatnonzero(np.diff(t)>=gap)+1,len(t)]
    ff,ee=f.copy(),e.copy()
    source_index=np.arange(len(t))
    segment_rows=[]
    for a,b in zip(starts[:-1],starts[1:]):
        n=b-a
        if n<3: continue
        local_dt=float(np.median(np.diff(t[a:b])))
        width=max(3,int(1./local_dt)//2*2+1)
        trend=median_filter(f[a:b],size=width,mode='reflect')
        if np.any(~np.isfinite(trend)) or np.any(trend<=0):
            raise ValueError('Invalid positive-flux surrogate trend')
        normalized=f[a:b]/trend
        normalized_error=e[a:b]/trend
        block_n=max(2,int(round(block_days/local_dt)))
        pieces=[]
        while sum(len(x) for x in pieces)<n:
            start=int(rng.integers(n))
            pieces.append((start+np.arange(min(block_n,n)))%n)
        idx=np.concatenate(pieces)[:n]
        ff[a:b]=normalized[idx]*trend
        ee[a:b]=normalized_error[idx]*trend
        source_index[a:b]=a+idx
        segment_rows.append(dict(start=int(a),stop=int(b),n=int(n),
                                 duration_days=float(t[b-1]-t[a]),block_cadences=int(block_n),
                                 nominal_block_days=float(block_n*local_dt),
                                 n_block_joins=int(np.ceil(n/min(block_n,n))-1)))
    return (t,ff,ee,np.zeros(len(t),dtype=int)),source_index,segment_rows


def main():
    OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
    selection=ROOT/'data/multisector/selection.csv'
    selected=[r for r in csv.DictReader(selection.open()) if int(r['sector'])==83][:8]
    if len(selected)!=8: raise ValueError('Need fixed eight pre-flux selected targets')
    hashes={p:sha(ROOT/p) for p in ['src/multisector_nulls.py','src/multisector_controls.py',
                                  'src/multisector_screen.py','src/fast_screen.py']}
    config=dict(seed=SEED,selection_sha256=sha(selection),targets=[int(r['tic']) for r in selected],
                sector=83,iterations_per_star=4,nominal_block_days=.4,methods_sha256=hashes,
                construction='Quality0 finite positive flux/error. Gap-defined segments. Normalize each segment by 1-day reflected running median; resample flux/error pairs in random circular 0.4-day blocks within the same segment; restore the original smooth baseline. Screen with the frozen two-width method.',
                limitations=['Not a valid survey-wide FAP calibration; only eight selected observed stars.',
                             'Stars were not selected to lack true events; astrophysical signal fragments remain.',
                             'Block joins introduce artifacts; correlations longer than blocks change.',
                             'Variability, gaps, heterogeneous noise and short segments limit stationarity.',
                             'Repeated surrogates of a star and the two detrenders are correlated.',
                             'Maximum is across reported refined peaks and two widths, not every possible model.'])
    manifest=OUT/'config.json'
    if manifest.exists() and json.loads(manifest.read_text())!=config:
        raise ValueError('Existing surrogate experiment differs; preserve it, use a new output run')
    manifest.write_text(json.dumps(config,indent=2)+'\n')
    dest=OUT/'iterations.jsonl'
    done={json.loads(line)['case_id'] for line in dest.read_text().splitlines()} if dest.exists() else set()
    observations=[]
    for row in selected:
        path=Path('/workspace/tess-research-data/multisector/raw/s0083')/Path(row['url']).name
        record=json.loads(path.with_suffix('.manifest.json').read_text())
        arrays,provenance=read_observation(record);observations.append(provenance)
        unchanged=screen_arrays(*arrays[:3],quality=arrays[3])
        (OUT/f'unchanged_tic{row["tic"]}.json').write_text(json.dumps(dict(observation=provenance,**unchanged),indent=2)+'\n')
        for iteration in range(4):
            case=f'{row["tic"]}:83:{iteration}'
            if case in done: continue
            seed=int(hashlib.sha256(f'{SEED}:{case}'.encode()).hexdigest()[:16],16)
            generated,index,segments=surrogate(arrays,np.random.default_rng(seed))
            prepared=CACHE/f'tic{row["tic"]}_i{iteration}.npz'
            np.savez_compressed(prepared,t=generated[0],f=generated[1],e=generated[2],q=generated[3],source_index=index)
            result=screen_arrays(*generated[:3],quality=generated[3])
            peaks=[dict(width_index=i,**p) for i,s in enumerate(result['scans']) for p in s['peaks']]
            retained=[p for p in peaks if not seed_reasons(p)]
            stats=dict(case_id=case,tic=int(row['tic']),sector=83,iteration=iteration,seed=seed,
                       surrogate_path=str(prepared),surrogate_sha256=sha(prepared),n_input=len(index),
                       segments=segments,max_local_snr=max((p['local_snr'] for p in peaks),default=0),
                       max_retained_local_snr=max((p['local_snr'] for p in retained),default=None),
                       any_retained=bool(retained),retained_peak_records=len(retained),scans=result['scans'])
            with dest.open('a')as stream: stream.write(json.dumps(stats,allow_nan=False)+'\n')
            print('SURROGATE',case,round(stats['max_local_snr'],2),len(retained),flush=True)
    (OUT/'observations.json').write_text(json.dumps(observations,indent=2)+'\n')
    results=[json.loads(s) for s in dest.read_text().splitlines()]
    snr=np.array([r['max_local_snr'] for r in results])
    summary=dict(generated_utc=datetime.now(timezone.utc).isoformat(),surrogate_searches=len(results),
                 distinct_original_stars=len({r['tic'] for r in results}),
                 searches_with_retained_seed=sum(r['any_retained'] for r in results),
                 max_local_snr_quantiles={str(q):float(np.quantile(snr,q)) for q in [0,.25,.5,.75,.9,.95,1]},
                 searches_with_max_local_snr_ge={str(c):int(np.sum(snr>=c)) for c in [5.5,6,7,8,10]},
                 per_star=[dict(tic=tic,n=sum(r['tic']==tic for r in results),
                                searches_with_retained_seed=sum(r['any_retained'] for r in results if r['tic']==tic),
                                max_local_snrs=[r['max_local_snr'] for r in results if r['tic']==tic])
                           for tic in sorted({r['tic'] for r in results})],
                 limitations=config['limitations'])
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('SUMMARY',json.dumps({k:summary[k] for k in ['surrogate_searches','distinct_original_stars','searches_with_retained_seed','max_local_snr_quantiles']}),flush=True)


if __name__=='__main__': main()
