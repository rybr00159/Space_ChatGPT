#!/usr/bin/env python3
"""Parallel executor for the unchanged, explicitly empirical-error alternative.

One parent writes SQLite; workers use the frozen estimator and screening
functions. The first complete record is compared with serial execution before
the remaining jobs run. Producer files and earlier accepted rows are immutable.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,hashlib,json,sqlite3,sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts'))
from run_multisector_batch import configure
import multisector_missing_errors as method
from multisector_combine_parallel import atomic_json,compare_nested


def analyze(job):
    tic,sector,raw,source_methods=job
    original=json.loads(raw)
    base={k:original[k]for k in ['tic','sector','ra','dec','path','url','sha256','raw_bytes']}
    base.update(original_status='insufficient_data',original_record_sha256=hashlib.sha256(raw.encode()).hexdigest(),
                original_main_methods=source_methods,source_sha256=method.signature(),checked_utc=method.utc())
    try:
        curve=method.read_curve(original['path'],tic,sector)
        if curve['sha256']!=original['sha256']:raise ValueError('Original raw product differs')
        base['producer_error_diagnostics']=method.diagnostics(curve)
        base.update(metadata=curve['metadata'],fits_checksums=curve['fits_checksums'],flux_column=curve['column'])
        if not base['producer_error_diagnostics']['recovery_eligible']:return tic,sector,'ineligible',base,None
        noise=method.estimate(curve);base['empirical_error']=noise
        result=method.screen_arrays(curve['t'],curve['f'],np.full(curve['e'].shape,noise['empirical_sigma_flux']),curve['q'])
        result.update(base)
        status='screened'if any(s['status']=='screened'for s in result['scans'])else'insufficient_data_after_recovery'
        result['status']=status
        return tic,sector,status,result,None
    except Exception as exc:return tic,sector,'error',base,type(exc).__name__+': '+str(exc)


def worker_init(name):configure(name)


def scientific(value):
    if isinstance(value,dict):return {k:scientific(v)for k,v in value.items()if k not in ['checked_utc','runtime_seconds']}
    if isinstance(value,(list,tuple)):return [scientific(v)for v in value]
    return value


def main(name,workers):
    configure(name)
    with sqlite3.connect(method.MAIN_DB.as_uri()+'?mode=ro',uri=True)as source:
        source_methods=dict(source.execute('SELECT key,value FROM metadata'))
        rows=source.execute("SELECT tic,sector,result_json FROM observations WHERE status='insufficient_data' AND result_json IS NOT NULL ORDER BY sector,tic").fetchall()
    conn=method.connect();done=set(conn.execute('SELECT tic,sector FROM observations'))
    jobs=[(tic,sector,raw,source_methods)for tic,sector,raw in rows if(tic,sector)not in done]
    execution=dict(batch=name,driver_sha256=method.sha(__file__),frozen_methods=method.signature(),
                   formula=method.ASSUMPTIONS,started_utc=method.utc(),workers=workers,planned_new_records=len(jobs))
    method.OUT.mkdir(parents=True,exist_ok=True)
    execution_path=method.OUT/'parallel_execution.json'
    if execution_path.exists():
        old=json.loads(execution_path.read_text())
        if old['driver_sha256']!=execution['driver_sha256']or old['frozen_methods']!=execution['frozen_methods']:
            raise ValueError('Existing execution driver differs; preserve and fork')
    if jobs:
        first=jobs.pop(0);serial=analyze(first)
        with ProcessPoolExecutor(max_workers=workers,initializer=worker_init,initargs=(name,))as pool:
            parallel=pool.submit(analyze,first).result()
            compared=compare_nested(scientific(serial),scientific(parallel),'first_full_recovery_record')
            execution['serial_parallel_exact_equivalence']=dict(tic=first[0],sector=first[1],numeric_values_compared=compared)
            method.save(conn,*serial[:3],serial[3],serial[4]);method.export(conn)
            atomic_json(execution_path,execution)
            futures=[pool.submit(analyze,job)for job in jobs]
            for i,future in enumerate(as_completed(futures),1):
                tic,sector,status,result,error=future.result();method.save(conn,tic,sector,status,result,error)
                if i%32==0 or i==len(jobs):method.export(conn)
    else:method.export(conn)
    execution['finished_utc']=method.utc();execution['final_counts']=dict(conn.execute('SELECT status,count(*) FROM observations GROUP BY status'))
    atomic_json(execution_path,execution);conn.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--name',required=True);parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args();main(args.name,args.workers)
