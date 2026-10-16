#!/usr/bin/env python3
"""Blind repeated-star searches for a named batch, with explicit input errors.

Uses the same frozen pooled BLS function as the first tranche. Prepared source
maps separate reported and estimated errors and preserve every exact input hash.
"""
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,csv,gzip,io,json,sqlite3,sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from run_multisector_batch import configure
from multisector_combine import search_prepared
from multisector_combine_parallel import atomic_json,atomic_bytes,sha,utc
from multisector_screen import WIDTHS,seed_reasons
from fast_screen import event_stats
from multisector_recognition import FreshCatalogIndex
import multisector_mixed_pool as preparation

INDEX=None
DEST=None


def initialize(snapshot,destination):
    global INDEX,DEST
    INDEX=FreshCatalogIndex(snapshot=ROOT/snapshot);DEST=Path(destination)


def compute(job):
    tic,rows,identity=job;scans=[]
    for width in WIDTHS:
        chunks=[]
        for source in identity['sources']:
            record=source['prepared'][str(width)];path=Path(record['path'])
            if sha(path)!=record['sha256']:raise ValueError('Prepared input differs')
            with np.load(path)as arrays:chunk={k:arrays[k]for k in ['t','f','e']}
            chunk['sector']=np.full(len(chunk['t']),source['sector']);chunks.append(chunk)
        joined={k:np.concatenate([c[k]for c in chunks])for k in ['t','f','e','sector']}
        order=np.argsort(joined['t']);joined={k:v[order]for k,v in joined.items()}
        known_keep=INDEX.transit_mask(tic,joined['t'],float(rows[0]['ra']),float(rows[0]['dec']))
        for label,mask in [('unmasked',np.ones(len(joined['t']),bool)),('known_masked_residual',known_keep)]:
            if label!='unmasked'and mask.all():continue
            if mask.sum()<50:continue
            t,f,e=(joined[k][mask]for k in ['t','f','e']);sectors=joined['sector'][mask]
            scan=search_prepared(t,f,e);scan.update(width_days=width,label=label,n_known_masked=int((~mask).sum()))
            for peak in scan['peaks']:
                peak['per_sector']={}
                for sector in sorted(set(sectors)):
                    take=sectors==sector
                    peak['per_sector'][str(sector)]=event_stats(t[take],f[take],e[take],peak['period'],peak['t0_btjd'],peak['duration_days'])
                peak['positive_support_sectors']=sum(s['n_positive_events']>=1 and s['local_snr']>=2.5 for s in peak['per_sector'].values())
            scans.append(scan)
    result=dict(tic=tic,sectors=[r['sector']for r in rows],identity=identity,scans=scans,checked_utc=utc(),
                limitations='Pooled-fitted statistics are conditional, not independent confirmations or FAP. Each source explicitly labels reported or reconstructed uncertainty.')
    path=DEST/f'tic{tic}.json'
    if path.exists():raise FileExistsError('Refusing to replace a completed result')
    atomic_json(path,result);return tic


def export(completed,destination,output,methods,snapshot,planned):
    index=FreshCatalogIndex(snapshot=ROOT/snapshot);table=[];sources=[]
    for tic in sorted(completed):
        path=destination/f'tic{tic}.json';result=json.loads(path.read_text())
        sources.append(dict(tic=tic,path=str(path),sha256=sha(path)))
        for scan in result['scans']:
            for rank,peak in enumerate(scan['peaks'],1):
                match=index.match_signal(tic,peak['period'])
                table.append(dict(tic=tic,sectors=';'.join(map(str,result['sectors'])),width_days=scan['width_days'],label=scan['label'],rank=rank,
                                  period_days=peak['period'],epoch_btjd=peak['t0_btjd'],duration_days=peak['duration_days'],depth_ppm=peak['depth']*1e6,
                                  local_snr=peak['local_snr'],n_events=peak['n_events'],n_positive_events=peak['n_positive_events'],
                                  positive_support_sectors=peak['positive_support_sectors'],seed_retained=not seed_reasons(peak),
                                  known_period_match=match['known_period_match'],known_harmonic_match=match['known_harmonic_match'],
                                  recognition_status=match['recognition_status'],evidence_path=str(path)))
    output.mkdir(parents=True,exist_ok=True)
    if table:
        stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
        atomic_bytes(output/'signals.csv.gz',gzip.compress(stream.getvalue().encode(),mtime=0))
    atomic_json(output/'summary.json',dict(checked_utc=utc(),methods=methods,catalog_snapshot=snapshot,planned_repeated_tics=planned,
                                          completed_repeated_tics=len(completed),signal_rows=len(table),sources=sources,
                                          both_sector_retained_rows=sum(r['seed_retained']and r['positive_support_sectors']>=2 for r in table),
                                          limitation='Repeated-star subset of this batch, searched regardless of individual-sector scores; not additional unique stars.'))
    print('BATCH POOLED',len(completed),'/',planned,flush=True)


def main(name,workers):
    cache,data,out=configure(name);preparation.CACHE=cache
    grouped=defaultdict(dict)
    for filename,assumption in [('survey.sqlite','producer_reported'),('error_recovery.sqlite','empirical_adjacent_difference_MAD')]:
        with sqlite3.connect((cache/filename).as_uri()+'?mode=ro',uri=True)as conn:
            for raw, in conn.execute("SELECT result_json FROM observations WHERE status='screened' ORDER BY tic,sector"):
                record=json.loads(raw);record['error_assumption']=assumption
                if record['sector']in grouped[record['tic']]:raise ValueError('Duplicate successful observation')
                grouped[record['tic']][record['sector']]=record
    repeated={tic:[rows[s]for s in sorted(rows)]for tic,rows in grouped.items()if len(rows)>=2}
    index=FreshCatalogIndex();snapshot=index.metadata()['fresh_snapshot']
    methods={p:sha(ROOT/p)for p in ['scripts/pool_multisector_batch.py','src/multisector_combine.py','src/multisector_mixed_pool.py','src/fast_screen.py','src/multisector_screen.py','src/multisector_missing_errors.py']}
    destination=cache/'combined';destination.mkdir(parents=True,exist_ok=True);output=out/'pooled'
    jobs=[];completed=set()
    for tic,rows in sorted(repeated.items()):
        identity=dict(methods=methods,sources=preparation.prepared_sources(tic,rows),catalog_snapshot=snapshot)
        path=destination/f'tic{tic}.json'
        if path.exists():
            if json.loads(path.read_text())['identity']!=identity:raise ValueError('Accepted result identity differs; preserve and fork')
            completed.add(tic)
        else:jobs.append((tic,rows,identity))
    export(completed,destination,output,methods,snapshot,len(repeated))
    with ProcessPoolExecutor(max_workers=workers,initializer=initialize,initargs=(snapshot,str(destination)))as pool:
        futures=[pool.submit(compute,job)for job in jobs]
        for i,future in enumerate(as_completed(futures),1):
            completed.add(future.result())
            if i%8==0 or i==len(jobs):export(completed,destination,output,methods,snapshot,len(repeated))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--name',required=True);parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args();main(args.name,args.workers)
