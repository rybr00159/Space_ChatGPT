#!/usr/bin/env python3
"""Execution-only all-observation pooled search of stars shared by two batches.

Reuses the frozen batch-pool initializer, compute and export unchanged. Source
preparation changes only its parent cache directory per observation. Scores and
period statistics are selected, conditional diagnostics, not discoveries/FAPs.
"""
from pathlib import Path
from collections import defaultdict,Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,gzip,json,sqlite3,sys,time,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import pool_multisector_batch as frozen
import multisector_mixed_pool as preparation
from multisector_combine_parallel import atomic_json,atomic_bytes,sha,utc
from multisector_recognition import FreshCatalogIndex
CACHE=Path('/workspace/tess-research-data/multisector')
DEST=CACHE/'cross_batch_combined'
OUTPUT=ROOT/'results/multisector_cross_batch'
BATCHES={'b01':CACHE,'b02_096_105':CACHE/'b02_096_105'}
METHOD_NAMES=['scripts/pool_multisector_batch.py','src/multisector_combine.py','src/multisector_mixed_pool.py','src/fast_screen.py','src/multisector_screen.py','src/multisector_missing_errors.py','scripts/pool_cross_batch.py']
SOURCE=Path(__file__).read_bytes()

def main(workers=2):
    if not 1<=workers<=2:raise ValueError('At most two workers are authorized for this cross-batch stage')
    started=time.monotonic();OUTPUT.mkdir(parents=True,exist_ok=True);DEST.mkdir(parents=True,exist_ok=True)
    methods={name:sha(ROOT/name)for name in METHOD_NAMES};snapdir=OUTPUT/'method_snapshots';snapdir.mkdir(exist_ok=True)
    for name,digest in methods.items():
        snapshot=snapdir/(digest+'.py');body=ROOT.joinpath(name).read_bytes()
        if snapshot.exists():
            if snapshot.read_bytes()!=body:raise ValueError('Method snapshot mismatch')
        else:snapshot.write_bytes(body)
    grouped=defaultdict(dict);batch_counts={};records_exact=[]
    for batch,cache in BATCHES.items():
        counts=Counter();tics=set()
        for filename,assumption in [('survey.sqlite','producer_reported'),('error_recovery.sqlite','empirical_adjacent_difference_MAD')]:
            db=cache/filename
            with sqlite3.connect(db.as_uri()+'?mode=ro',uri=True)as conn:
                rows=conn.execute("SELECT result_json FROM observations WHERE status='screened' ORDER BY tic,sector").fetchall()
            for raw, in rows:
                record=json.loads(raw);tic=int(record['tic']);sector=int(record['sector'])
                if sector in grouped[tic]:raise ValueError('Duplicate successful TIC/sector across batches/assumptions')
                record.update(error_assumption=assumption,cross_batch_name=batch,cross_batch_cache=str(cache),source_database=str(db),source_database_record_sha256=hashlib.sha256(raw.encode()).hexdigest())
                grouped[tic][sector]=record;tics.add(tic);counts[assumption]+=1
                records_exact.append(dict(batch=batch,database=str(db),raw_result_json=raw,raw_result_json_sha256=hashlib.sha256(raw.encode()).hexdigest(),tic=tic,sector=sector,assumption=assumption))
        batch_counts[batch]=dict(successful_observations=sum(counts.values()),unique_stars=len(tics),error_assumptions=dict(counts))
    shared={tic:[values[s]for s in sorted(values)]for tic,values in grouped.items()if len({r['cross_batch_name']for r in values.values()})==len(BATCHES)and len(values)>=2}
    index=FreshCatalogIndex();snapshot=index.metadata()['fresh_snapshot']
    plans=[];jobs=[];completed=set();raw_verified=set();prepared_verified=set()
    for tic,rows in sorted(shared.items()):
        sources=[]
        for record in rows:
            rawpath=Path(record['path'])
            if sha(rawpath)!=record['sha256']:raise ValueError('Source raw observation changed')
            raw_verified.add(str(rawpath));preparation.CACHE=Path(record['cross_batch_cache'])
            source=preparation.prepared_sources(tic,[record])[0]
            source.update(batch=record['cross_batch_name'],source_database=record['source_database'],source_database_record_sha256=record['source_database_record_sha256'])
            for width,prepared in source['prepared'].items():
                path=Path(prepared['path'])
                if sha(path)!=prepared['sha256']:raise ValueError('Prepared source changed')
                with np.load(path)as a:
                    if not all(np.all(np.isfinite(a[k]))for k in ['t','f','e'])or not np.all(a['e']>0):raise ValueError('Invalid prepared scientific arrays')
                prepared_verified.add(str(path))
            sources.append(source)
        identity=dict(methods=methods,sources=sources,catalog_snapshot=snapshot,cross_batch_scope='All successful observations of stars present in both batches; no score-based selection.')
        plan=dict(tic=tic,sectors=[r['sector']for r in rows],batches=[r['cross_batch_name']for r in rows],sources=sources,identity=identity);plans.append(plan)
        path=DEST/f'tic{tic}.json'
        if path.exists():
            if json.loads(path.read_text())['identity']!=identity:raise ValueError('Existing cross-batch identity differs; preserve and fork')
            completed.add(tic)
        else:jobs.append((tic,rows,identity))
    chosen=set(shared);original=[r for r in records_exact if r['tic']in chosen]
    body=''.join(json.dumps(r,separators=(',',':'))+'\n'for r in original).encode();atomic_bytes(OUTPUT/'input_database_rows.jsonl.gz',gzip.compress(body,mtime=0))
    scope=dict(created_utc=utc(),batch_input_counts=batch_counts,shared_unique_stars=len(shared),shared_observations=sum(len(v)for v in shared.values()),sector_counts=dict(Counter(r['sector']for rows in shared.values()for r in rows)),error_assumption_counts=dict(Counter(r['error_assumption']for rows in shared.values()for r in rows)),observations_per_star=dict(Counter(len(rows)for rows in shared.values())),raw_sha_verified=len(raw_verified),prepared_sha_verified=len(prepared_verified),methods=methods,catalog_snapshot=snapshot,selected_tics=sorted(shared),plans=plans,limitations=['This overlaps previously counted stars and observations; do not add to survey unique-star totals.','All available successful observations are combined irrespective of individual-sector or earlier pooled scores.','The frozen search remains period0.5..40days and inherits all its sensitivity and systematic limits.','Reported and reconstructed errors are tagged separately; empirical adjacent differences do not calibrate red noise.','Pooled selection, detrenders, masks and per-sector support do not provide independent confirmation or a calibrated false-alarm probability.'])
    atomic_json(OUTPUT/'input_scope.json',scope)
    frozen.export(completed,DEST,OUTPUT,methods,snapshot,len(shared))
    execution=dict(started_utc=utc(),workers=workers,initial_completed=len(completed),pending=len(jobs),method=methods,cache_destination=str(DEST),status='running',scope=dict(unique_stars=len(shared),observations=scope['shared_observations']))
    atomic_json(OUTPUT/'execution.json',execution)
    failures=[]
    with ProcessPoolExecutor(max_workers=workers,initializer=frozen.initialize,initargs=(snapshot,str(DEST)))as pool:
        futures={pool.submit(frozen.compute,job):job[0]for job in jobs}
        for i,future in enumerate(as_completed(futures),1):
            tic=futures[future]
            try:completed.add(future.result());print('CROSS_DONE',tic,len(completed),'/',len(shared),flush=True)
            except Exception as exc:
                failures.append(dict(tic=tic,error=repr(exc),checked_utc=utc()));print('CROSS_FAILED',tic,repr(exc),flush=True)
            if i%4==0 or i==len(jobs):
                frozen.export(completed,DEST,OUTPUT,methods,snapshot,len(shared));execution.update(completed=len(completed),failures=failures,elapsed_seconds=time.monotonic()-started,last_checked_utc=utc());atomic_json(OUTPUT/'execution.json',execution)
    for name,digest in methods.items():
        if sha(ROOT/name)!=digest:raise ValueError('Scientific or execution source changed while stage ran')
    execution.update(status='complete'if len(completed)==len(shared)and not failures else'partial_failed',completed=len(completed),failures=failures,elapsed_seconds=time.monotonic()-started,finished_utc=utc());atomic_json(OUTPUT/'execution.json',execution)
    frozen.export(completed,DEST,OUTPUT,methods,snapshot,len(shared));print('CROSS_FINISHED',execution,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);args=parser.parse_args();main(args.workers)
