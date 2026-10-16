"""Conditional raw-flux checks of the final blind pooled-search queues.

Selection is recorded before opening these targets' raw flux: top eight distinct
stars per queue with retained seeds, >=2 supported sectors and depth <5000ppm.
Already characterized TIC229146531 and TIC238182368 remain in their own dossiers.
A pooled score or a warning is not a calibrated discovery/false-positive label.
"""
from pathlib import Path
import hashlib,json,os,sys,traceback
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import pandas as pd
import multisector_followup as f
from pooled_south_vet import incomplete_observation
OUT=ROOT/'results/multisector_final_pool_followup'
QUEUES=['results/multisector_cross_batch/signals.csv.gz','results/multisector_batches/b02_096_105/pooled/signals.csv.gz']

def main():
    OUT.mkdir(parents=True,exist_ok=True);selection=[]
    for queue in QUEUES:
        d=pd.read_csv(ROOT/queue);d=d[d.seed_retained.eq(True)&d.positive_support_sectors.ge(2)&d.depth_ppm.lt(5000)]
        d=d.sort_values('local_snr',ascending=False).drop_duplicates('tic').head(8)
        for row in d.to_dict('records'):
            if row['tic']not in [229146531,238182368]:selection.append(dict(queue=queue,**row))
    selected=OUT/'selection.json'
    if not selected.exists():f.save(selected,dict(checked_utc=f.utc(),script_sha256=f.digest(__file__),selection=selection,criteria=__doc__))
    else:
        prior=json.loads(selected.read_text())
        if prior['script_sha256']!=f.digest(__file__):raise ValueError('Driver changed after selection')
        selection=prior['selection']
    rows=[]
    for peak in selection:
        tic=int(peak['tic']);path=OUT/f'tic{tic}.json';source=Path(peak['evidence_path'])
        if path.exists():
            output=json.loads(path.read_text())
            if output['pooled_input_sha256']!=f.digest(source):raise ValueError('Pooled input changed')
        else:
            pool=json.loads(source.read_text());ephem=dict(period=peak['period_days'],epoch_btjd=peak['epoch_btjd'],duration_days=peak['duration_days'])
            output=dict(tic=tic,checked_utc=f.utc(),selection=peak,ephemeris=ephem,pooled_input_path=str(source),pooled_input_sha256=f.digest(source),script_sha256=f.digest(__file__),observations=[],errors=[],limitations=['Selection and ephemeris used these sectors; these are conditional checks, not independent confirmation.','No calibrated false-alarm probability, planet classification or novelty inference.'])
            for s in pool['identity']['sources']:
                record=dict(tic=tic,sector=s['sector'],path=s['path'],sha256=s['raw_sha256'],error_assumption=s['error_assumption'])
                if s['error_assumption']=='empirical_adjacent_difference_MAD':record.update(empirical_sigma_flux=s['empirical_sigma_flux'],error_estimator_sha256=s['estimator_sha256'])
                try:r=f.vet_observation(record,peak,output=False,fixed_ephemeris=ephem)
                except Exception as exc:
                    output['errors'].append(dict(sector=record['sector'],error=repr(exc),traceback=traceback.format_exc()))
                    try:r=incomplete_observation(record,ephem,exc)
                    except Exception as fallback:
                        output['errors'].append(dict(sector=record['sector'],fallback_error=repr(fallback)));continue
                output['observations'].append(r)
            f.save(path,output)
        for r in output['observations']:
            cols=r.get('column_tests',{});corr=cols.get('corrected',{}).get('primary',{});raw=cols.get('raw',{}).get('primary',{})
            row=dict(tic=tic,sector=r['sector'],period_days=peak['period_days'],pooled_local_snr=peak['local_snr'],pooled_depth_ppm=peak['depth_ppm'],corrected_local_snr=corr.get('snr'),raw_local_snr=raw.get('snr'),raw_depth_ppm=raw.get('depth',0)*1e6 if raw.get('depth')is not None else None,n_events=raw.get('n_events'),warnings=';'.join(r.get('warnings',[])),error_assumption=r.get('source',{}).get('error_assumption'),evidence=str(path.relative_to(ROOT)))
            rows.append(row);print('FINAL_POOL',row,flush=True)
        pd.DataFrame(rows).to_csv(OUT/'observations.csv',index=False)
    f.save(OUT/'summary.json',dict(checked_utc=f.utc(),selected_stars=len(selection),observations=len(rows),incomplete_stars=sum(bool(json.loads((OUT/f'tic{int(p["tic"])}.json').read_text())['errors'])for p in selection),script_sha256=f.digest(__file__),selection_sha256=f.digest(selected),scope=__doc__))

if __name__=='__main__':main()
