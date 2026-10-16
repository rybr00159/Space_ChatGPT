"""Selected shallow-lead follow-up for the new b02_096_105 survey batch.

The discovery and empirical-error databases are opened read-only. Additional
sectors first receive the original discovery ephemeris; no sector retuning is
allowed to count as confirmation. Diagnostic scores remain conditional.
"""
from pathlib import Path
import argparse, json, os, sqlite3, sys, traceback
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import pandas as pd
import multisector_followup as f
from pooled_south_vet import explicit_followup_error_policy, incomplete_observation, compact, transit_pixels

BATCH='b02_096_105'
DBDIR=f.CACHE/BATCH
TABLES=ROOT/'results/multisector_batches'/BATCH
DEFAULT={229146531:.934911,71274048:1.359244,147271811:.574451,
         389752955:10.458554,80204331:6.081661,80427865:2.165521}

def records(tic):
    result={}
    for name,empirical in [('survey.sqlite',False),('error_recovery.sqlite',True)]:
        con=sqlite3.connect('file:'+str(DBDIR/name)+'?mode=ro',uri=True)
        for raw, in con.execute('SELECT result_json FROM observations WHERE tic=? AND result_json IS NOT NULL',(tic,)):
            row=json.loads(raw)
            if empirical:
                if row.get('status')!='screened':continue
                sha=row['source_sha256']['src/multisector_missing_errors.py']
                if sha!=f.digest(ROOT/'src/multisector_missing_errors.py'):raise ValueError('Empirical estimator source changed')
                row.update(error_assumption='empirical_adjacent_difference_MAD',empirical_sigma_flux=row['empirical_error']['empirical_sigma_flux'],error_estimator_sha256=sha)
            result[int(row['sector'])]=row
        con.close()
    return result

def select(tic,period):
    producer=pd.read_csv(TABLES/'signals.csv');empirical=pd.read_csv(TABLES/'error_recovery/signals.csv')
    producer['analysis_assumption']='producer_errors';empirical['analysis_assumption']='empirical_adjacent_difference_MAD'
    table=pd.concat([producer,empirical],ignore_index=True);table=table[table.tic.eq(tic)]
    table=table.assign(distance=abs(table.period_days-period))
    row=table.sort_values(['distance','local_snr'],ascending=[True,False]).iloc[0].to_dict()
    ephem=dict(period=float(row['period_days']),epoch_btjd=float(row['epoch_btjd']),duration_days=float(row['duration_days']))
    return row,ephem

def run_one(tic,period,followup):
    peak,ephem=select(tic,period);available=records(tic);source=available[int(peak['sector'])]
    output=dict(tic=tic,batch=BATCH,status='in_progress',checked_utc=f.utc(),selection=peak,ephemeris=ephem,
                methods={str(p):f.digest(ROOT/p) for p in ['scripts/batch02_followup.py','scripts/pooled_south_vet.py','src/multisector_followup.py','src/multisector_missing_errors.py']},
                selection_rationale='Nominal R<1.2 solar, selected depth<5000ppm, at least3 events, positive fraction>=.7, odd/even score<4, max-event fraction<.6, and relatively short fractional event duration. Catalog status remains separate.',
                discovery=[],other_batch=[],additional=[],incomplete_observations=[],errors=[],
                limitations=['Single discovery-sector search selected these ephemerides; no calibrated survey false-alarm probability.','All additional sectors first receive the unchanged discovery ephemeris.','Warning counts are not a physical classification or strict automatic planet veto.','Same data through different detrenders are not independent evidence.'])
    path=f.OUT/f'tic{tic}'/'batch02_followup.json'
    def evaluate(record,group):
        try:
            result=f.vet_observation(record,peak,fixed_ephemeris=ephem)
            output[group].append(compact(result))
        except Exception as exc:
            output['errors'].append(dict(stage=group,sector=record['sector'],error=repr(exc),traceback=traceback.format_exc()))
            try:output['incomplete_observations'].append(incomplete_observation(record,ephem,exc))
            except Exception as fallback:output['errors'].append(dict(stage=group+'_local_fallback',sector=record['sector'],error=repr(fallback)))
        f.save(path,output)
    evaluate(source,'discovery')
    for sector,record in available.items():
        if sector!=int(source['sector']):evaluate(record,'other_batch')
    if followup:
        try:
            retrieved=f.retrieve_independent(source,max_sectors=followup,exclude_sectors=list(available))
            output['retrieval']=retrieved
            for record in retrieved['observations']:
                if record.get('status')=='retained':evaluate(explicit_followup_error_policy(record),'additional')
        except Exception as exc:output['errors'].append(dict(stage='retrieval',error=repr(exc),traceback=traceback.format_exc()))
    output['status']='complete_with_preserved_errors' if output['errors'] else 'complete'
    f.save(path,output)
    print('BATCH02',tic,ephem,[(r['sector'],round(r['primary']['corrected']['snr'],2),round(r['primary']['raw']['snr'],2),r['warnings'])for k in ['discovery','other_batch','additional']for r in output[k]],'errors',output['errors'],flush=True)
    return output

def joint_retest(tic,summary):
    selected=[r['source']for group in ['discovery','other_batch']for r in summary[group]]
    independent=len(selected)>=2
    if not independent:selected += [r['source']for r in summary['additional']]
    if len(selected)<2:
        print('JOINT_UNAVAILABLE',tic,'fewer than two observations',flush=True);return
    name='batch02_discovery_joint' if independent else 'batch02_conditional_all_sector_joint'
    joint=f.joint_recurrence(selected,summary['ephemeris'],name=name)
    ephem=dict(period=joint['period_days'],epoch_btjd=joint['epoch_btjd'],duration_days=joint['duration_days'])
    peak=dict(summary['selection'],period_days=ephem['period'])
    output=dict(tic=tic,ephemeris=ephem,original_ephemeris=summary['ephemeris'],joint_source=name+'.json',retests=[],errors=[],script_sha256=f.digest(__file__),interpretation='One common discovery-period refinement followed by unchanged additional-sector ephemeris. Additional observations were already inspected at the original ephemeris; this is conditional follow-up, not blind confirmation.' if independent else 'Common optimization includes all inspected observations. Neither independent confirmation nor calibrated significance.')
    for item in summary['additional']:
        try:output['retests'].append(compact(f.vet_observation(item['source'],peak,fixed_ephemeris=ephem)))
        except Exception as exc:output['errors'].append(dict(sector=item['sector'],error=repr(exc),traceback=traceback.format_exc()))
    f.save(f.OUT/f'tic{tic}'/'batch02_joint_retest.json',output)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tic',type=int);parser.add_argument('--period',type=float)
    parser.add_argument('--followup',type=int,default=2)
    parser.add_argument('--joint',action='store_true');parser.add_argument('--pixels',action='store_true');args=parser.parse_args()
    targets={args.tic:args.period or DEFAULT[args.tic]}if args.tic else DEFAULT
    for tic,period in targets.items():
        if args.joint or args.pixels:
            summary=json.loads((f.OUT/f'tic{tic}'/'batch02_followup.json').read_text())
            if summary.get('status')not in ['complete','complete_with_preserved_errors']:
                raise ValueError('Complete the source follow-up checkpoint before dependent joint/pixel analysis')
            if args.joint:joint_retest(tic,summary)
            if args.pixels:transit_pixels(tic,summary,sector=summary['discovery'][0]['sector'])
        else:run_one(tic,period,args.followup)

if __name__=='__main__':main()
