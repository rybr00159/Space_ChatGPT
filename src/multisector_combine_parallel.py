"""Parallel execution driver for the unchanged, separately frozen pooled search.

The scientific function, period grid, data joins, masks and statistics match
multisector_combine.py exactly. Driver identity is additional provenance; it
does not replace or relax the original per-result scientific identity.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import asdict
import argparse,csv,gzip,hashlib,io,json,os,sqlite3,time
import numpy as np

from multisector_combine import search_prepared,CONFIG
from fast_screen import event_stats
from multisector_screen import WIDTHS,seed_reasons
from multisector_recognition import FreshCatalogIndex

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/workspace/tess-research-data/multisector')
DEST=CACHE/'combined'
OUT=ROOT/'results/multisector'
METHODS=['src/multisector_combine.py','src/multisector_screen.py','src/fast_screen.py']
INDEX=None


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def utc():return datetime.now(timezone.utc).isoformat()


def atomic_bytes(path,body):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+f'.tmp.{os.getpid()}')
    with temporary.open('wb')as stream:
        stream.write(body);stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)


def atomic_json(path,value):
    atomic_bytes(path,(json.dumps(value,indent=2,allow_nan=False)+'\n').encode())


def original_methods():return {p:sha(ROOT/p)for p in METHODS}


def snapshot_records():
    with sqlite3.connect((CACHE/'survey.sqlite').as_uri()+'?mode=ro',uri=True)as conn:
        conn.execute('BEGIN')
        counts=dict(conn.execute('SELECT status,COUNT(*) FROM observations GROUP BY status'))
        records=[json.loads(r[0])for r in conn.execute("SELECT result_json FROM observations WHERE status='screened' ORDER BY tic,sector")]
    grouped=defaultdict(list)
    for record in records:grouped[int(record['tic'])].append(record)
    return {tic:rows for tic,rows in grouped.items()if len(rows)>=2},counts


def expected_identity(rows,index,methods):
    return dict(methods=methods,sources=[dict(sector=r['sector'],path=r['path'],sha256=r['sha256'])for r in rows],
                catalog_snapshot=index.metadata()['fresh_snapshot'])


def initialize_worker(snapshot):
    global INDEX
    INDEX=FreshCatalogIndex(snapshot=ROOT/snapshot)


def compute(tic,rows,identity,index):
    """Same scientific operations and ordering as frozen base run(); no writes."""
    scans=[]
    for width in WIDTHS:
        chunks=[]
        for record in rows:
            path=CACHE/'prepared'/f'tic{tic}_s{record["sector"]:04d}_w{width:g}.npz'
            with np.load(path)as data:chunk={k:data[k]for k in ['t','f','e']}
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
    return dict(tic=tic,sectors=[r['sector']for r in rows],identity=identity,scans=scans,
                checked_utc=utc(),limitations='Pooled ephemeris optimized on all sectors; per-sector scores are conditional, not holdouts. Known masks applied to independently detrended arrays.')


def worker(job):
    tic,rows,identity=job
    started=time.monotonic()
    try:
        result=compute(tic,rows,identity,INDEX)
        path=DEST/f'tic{tic}.json'
        # The main process checked existing accepted results before scheduling.
        # If another process creates one, stop instead of replacing its result.
        if path.exists():raise FileExistsError('Accepted result appeared during computation: '+str(path))
        atomic_json(path,result)
        return dict(tic=tic,path=str(path),sha256=sha(path),seconds=time.monotonic()-started,status='completed')
    except Exception as exc:
        failure=dict(tic=tic,identity=identity,error=type(exc).__name__+': '+str(exc),checked_utc=utc())
        atomic_json(DEST/f'tic{tic}.error.json',failure)
        return dict(tic=tic,error=failure['error'],status='error',seconds=time.monotonic()-started)


def accepted(path,identity,tic):
    if not path.exists():return None
    try:
        result=json.loads(path.read_text())
    except json.JSONDecodeError:
        interrupted=path.with_name(path.name+'.interrupted_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        path.rename(interrupted)
        print('PRESERVED PARTIAL',str(interrupted),flush=True)
        return None
    if not all(k in result for k in ['tic','sectors','identity','scans','checked_utc','limitations']):
        raise ValueError('Valid JSON with incomplete result structure; preserve for review: '+str(path))
    if result['tic']!=tic or result['identity']!=identity:
        raise ValueError('Accepted scientific identity differs; preserve and fork rather than replace: '+str(path))
    return result


def compare_nested(a,b,path='result'):
    """Exact numeric comparison, with no tolerance capable of hiding drift."""
    if type(a)!=type(b):raise AssertionError('Changed type at '+path)
    if isinstance(a,dict):
        if a.keys()!=b.keys():raise AssertionError('Changed keys at '+path)
        return sum(compare_nested(a[k],b[k],path+'.'+k)for k in a)
    if isinstance(a,list):
        if len(a)!=len(b):raise AssertionError('Changed length at '+path)
        return sum(compare_nested(x,y,path+f'[{i}]')for i,(x,y)in enumerate(zip(a,b)))
    if a!=b:raise AssertionError('Scientific mismatch at '+path+': '+repr(a)+' != '+repr(b))
    return int(isinstance(a,(int,float))and not isinstance(a,bool))


def verify_one():
    repeated,counts=snapshot_records();index=FreshCatalogIndex();methods=original_methods()
    candidates=[]
    for tic,rows in sorted(repeated.items()):
        path=DEST/f'tic{tic}.json'
        if path.exists():candidates.append((tic,rows,path))
    if not candidates:raise RuntimeError('No accepted base result available for equivalence test')
    tic,rows,path=candidates[0]
    identity=expected_identity(rows,index,methods)
    before=accepted(path,identity,tic)
    original_hash=sha(path);started=time.monotonic()
    after=compute(tic,rows,identity,index)
    numeric=compare_nested({k:v for k,v in before.items()if k!='checked_utc'},
                           {k:v for k,v in after.items()if k!='checked_utc'})
    if sha(path)!=original_hash:raise AssertionError('Verification changed accepted result')
    report=dict(checked_utc=utc(),passed=True,tic=tic,base_result_path=str(path),base_result_sha256=original_hash,
                exact_numeric_values_compared=numeric,scans=len(before['scans']),
                peaks=sum(len(s['peaks'])for s in before['scans']),methods=methods,
                execution_driver_sha256=sha(Path(__file__)),seconds=time.monotonic()-started,
                scope='Every scientific JSON field compared exactly, including periods/epochs/depths/SNR/event statistics and per-sector statistics; only timestamp excluded.')
    atomic_json(OUT/'combined_parallel_equivalence.json',report)
    print('EXACT EQUIVALENCE PASSED',tic,numeric,'numeric values',flush=True)
    return report


def export(repeated,completed,index,methods,execution,errors):
    summaries=[];table=[]
    for tic in sorted(completed):
        rows=repeated[tic];path=DEST/f'tic{tic}.json'
        result=accepted(path,expected_identity(rows,index,methods),tic)
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
    if table:
        text=io.StringIO(newline='');writer=csv.DictWriter(text,fieldnames=list(table[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(table)
        # The uncompressed table is retained for existing follow-up consumers.
        body=text.getvalue().encode()
        atomic_bytes(OUT/'combined_signals.csv',body)
        atomic_bytes(OUT/'combined_signals.csv.gz',gzip.compress(body,mtime=0))
    report=dict(checked_utc=utc(),methods=methods,repeated_stars_searched=len(summaries),
                repeated_stars_in_snapshot=len(repeated),pending_stars=len(repeated)-len(completed)-len(errors),
                error_stars=len(errors),errors=errors,sector_observations_in_pooled_search=sum(len(r['sectors'])for r in summaries),
                signal_rows=len(table),both_sector_retained_rows=sum(r['seed_retained']and r['positive_support_sectors']>=2 for r in table),
                sources=summaries,execution=execution,
                limitation='Selection includes all targets with at least two screened sectors in the fixed input snapshot, regardless of single-sector SNR. Insufficient-data observations are not called searched. Pooled counts overlap the primary survey.')
    atomic_json(OUT/'combined_summary.json',report)


def run(workers=3,require_finalized=4096):
    repeated,counts=snapshot_records()
    if sum(counts.values())<require_finalized:
        raise RuntimeError(f'Main survey is still running: {sum(counts.values())}/{require_finalized} finalized records; do not launch bulk pool yet')
    methods=original_methods();driver=sha(Path(__file__))
    verification=json.loads((OUT/'combined_parallel_equivalence.json').read_text())
    if not verification['passed']or verification['methods']!=methods or verification['execution_driver_sha256']!=driver:
        raise ValueError('Current exact-equivalence verification is required before parallel run')
    index=FreshCatalogIndex();snapshot=index.metadata()['fresh_snapshot']
    DEST.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    prepared=[];jobs=[];completed=set();reused=0
    for tic,rows in sorted(repeated.items()):
        identity=expected_identity(rows,index,methods)
        if accepted(DEST/f'tic{tic}.json',identity,tic)is not None:
            completed.add(tic);reused+=1
        else:
            slim=[{k:r[k]for k in ['tic','sector','path','sha256','ra','dec']}for r in rows]
            jobs.append((tic,slim,identity))
        for row in rows:
            if sha(row['path'])!=row['sha256']:raise ValueError('Raw observation hash mismatch')
            for width in WIDTHS:
                p=CACHE/'prepared'/f'tic{tic}_s{row["sector"]:04d}_w{width:g}.npz'
                prepared.append(dict(tic=tic,sector=row['sector'],width_days=width,path=str(p),sha256=sha(p)))
    identity_manifest=dict(methods=methods,catalog_snapshot=snapshot,
                           prepared=prepared,sources=[dict(tic=tic,**expected_identity(rows,index,methods))for tic,rows in sorted(repeated.items())])
    manifest_path=DEST/'parallel_input_snapshot.json'
    if manifest_path.exists()and json.loads(manifest_path.read_text())!=identity_manifest:
        raise ValueError('Parallel input snapshot differs; preserve and fork run')
    atomic_json(manifest_path,identity_manifest)
    execution=dict(driver='src/multisector_combine_parallel.py',driver_sha256=driver,workers=workers,
                   input_snapshot_path=str(manifest_path),input_snapshot_sha256=sha(manifest_path),
                   finalized_primary_observation_counts=counts,reused_matching_base_results=reused,
                   equivalence_report=str(OUT/'combined_parallel_equivalence.json'),
                   separate_science_and_execution_identity=True)
    atomic_json(OUT/'combined_execution.json',execution)
    errors=[];export(repeated,completed,index,methods,execution,errors)
    print('POOLED PARALLEL START',len(repeated),'stars;',reused,'exact-identity checkpoints reused;',len(jobs),'pending',flush=True)
    with ProcessPoolExecutor(max_workers=workers,initializer=initialize_worker,initargs=(snapshot,))as pool:
        futures={pool.submit(worker,job):job[0]for job in jobs}
        for i,future in enumerate(as_completed(futures),1):
            result=future.result()
            if result['status']=='completed':completed.add(result['tic'])
            else:errors.append(result)
            if i%16==0 or i==len(jobs):
                export(repeated,completed,index,methods,execution,errors)
                print('POOLED PARALLEL',len(completed),'/',len(repeated),'errors',len(errors),flush=True)
    export(repeated,completed,index,methods,execution,errors)
    print('POOLED PARALLEL COMPLETE',len(completed),'stars;',len(errors),'errors',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['verify','run']);p.add_argument('--workers',type=int,default=3)
    p.add_argument('--require-finalized',type=int,default=4096);args=p.parse_args()
    if args.stage=='verify':verify_one()
    else:run(args.workers,args.require_finalized)
