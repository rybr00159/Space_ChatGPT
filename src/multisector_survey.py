"""Checkpointed bright-dwarf survey across multiple independently selected sectors.

Never writes to the original sector-1 database. Preserves every fetched raw file,
including errors, and source hashes prohibit mixed scientific methods in a run.
"""
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path
from collections import Counter
import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import threading
import time

import requests

from multisector_screen import analyze_file, seed_reasons, CONFIG, WIDTHS

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path('/workspace/tess-research-data/multisector')
DATA = ROOT/'data/multisector'
OUT = ROOT/'results/multisector'
DB = CACHE/'survey.sqlite'
LOCAL = threading.local()


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def signature():
    return {name:sha(ROOT/name) for name in ['src/fast_screen.py','src/multisector_screen.py']}


def url_for(tic, sector):
    s = str(tic).zfill(16)
    sub = '/'.join(s[i:i+4] for i in range(0,16,4))
    return f'https://archive.stsci.edu/hlsps/qlp/s{sector:04d}/{sub}/hlsp_qlp_tess_ffi_s{sector:04d}-{s}_tess_v01_llc.fits'


def plan(metadata, per_sector):
    """Pre-flux selection; hash-order metadata-eligible rows, never flux-ranked."""
    DATA.mkdir(parents=True, exist_ok=True)
    destination = DATA/'selection.csv'
    manifest = DATA/'selection.json'
    if destination.exists():
        m = json.loads(manifest.read_text())
        if m['per_sector'] != per_sector or sha(metadata) != m['metadata_sha256']:
            raise ValueError('Existing fixed selection differs; use a new named run for expansion')
        if sha(destination) != m['selection_sha256']:
            raise ValueError('Selection checksum mismatch')
        return list(csv.DictReader(destination.open()))
    rows = list(csv.DictReader(metadata.open()))
    groups = {}
    for row in rows:
        tic,sector = int(row['tic']),int(row['sector'])
        # Metadata are retrieved before any search flux. No missing values pass.
        try:
            eligible = 0<float(row['tmag'])<=12.5 and 3000<=float(row['teff'])<=6500 and .15<=float(row['radius'])<=2.5
        except (ValueError,TypeError):
            eligible = False
        if eligible:
            groups.setdefault(sector,[]).append(row)
    selected = []
    for sector in sorted(groups):
        ranked = sorted(groups[sector], key=lambda r:hashlib.sha256(('multisector-v1:'+str(int(r['tic']))).encode()).hexdigest())
        for row in ranked[:per_sector]:
            selected.append(dict(order=len(selected),sector=sector,tic=int(row['tic']),ra=row['ra'],dec=row['dec'],
                                 tmag=row['tmag'],teff=row['teff'],radius=row['radius'],url=url_for(int(row['tic']),sector)))
    with destination.open('w',newline='')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(selected[0]));writer.writeheader();writer.writerows(selected)
    m = dict(created_utc=utc(),metadata_path=str(metadata),metadata_sha256=sha(metadata),
             per_sector=per_sector,eligible_counts={str(k):len(v)for k,v in groups.items()},
             selected_counts=dict(Counter(str(r['sector'])for r in selected)),
             unique_selected_stars=len({r['tic']for r in selected}),selected_observations=len(selected),
             selection_sha256=sha(destination),config=CONFIG,widths=WIDTHS,
             rule='SHA256 multisector-v1:TIC order within each sector; Tmag<=12.5,Teff3000..6500K,R.15..2.5Rsun; pre-flux selection',
             note='Requested count is not actual completed search count. Repeated TICs are counted separately only as sector observations.')
    manifest.write_text(json.dumps(m,indent=2)+'\n')
    return selected


def download(row):
    record = {**row, 'tic':int(row['tic']), 'sector':int(row['sector'])}
    path = CACHE/'raw'/f's{record["sector"]:04d}'/Path(record['url']).name
    meta = path.with_suffix('.manifest.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        m = json.loads(meta.read_text())
        if sha(path) != m['sha256']:
            raise ValueError('Retained raw data hash mismatch')
    else:
        if shutil.disk_usage(CACHE).free < 5*1024**3:
            raise RuntimeError('Storage reserve reached; no observations deleted')
        if not hasattr(LOCAL,'session'):
            LOCAL.session=requests.Session()
            LOCAL.session.headers.update({'User-Agent':'Space_ChatGPT-public-TESS-multisector/1.0'})
        for attempt in range(3):
            response=LOCAL.session.get(record['url'],timeout=(15,75))
            body=response.content
            m=dict(url=record['url'],path=str(path),tic=record['tic'],sector=record['sector'],
                   http_status=response.status_code,bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),retrieved_utc=utc())
            if response.status_code==200 and body.startswith(b'SIMPLE  '):
                path.write_bytes(body);meta.write_text(json.dumps(m,indent=2)+'\n');break
            failure=path.with_name(path.name+f'.failed{attempt}')
            failure.write_bytes(body)
            failure.with_suffix(failure.suffix+'.json').write_text(json.dumps(m,indent=2)+'\n')
            if response.status_code not in [429,500,502,503,504] or attempt==2:
                raise RuntimeError('Archive HTTP '+str(response.status_code)+'; response preserved')
            time.sleep(2**attempt)
    record.update(path=str(path),sha256=m['sha256'],raw_bytes=m['bytes'],prepared_dir=str(CACHE/'prepared'))
    return record


def connect():
    CACHE.mkdir(parents=True,exist_ok=True)
    conn=sqlite3.connect(DB)
    conn.execute('PRAGMA journal_mode=WAL')
    conn.execute('CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY,value TEXT)')
    conn.execute('CREATE TABLE IF NOT EXISTS observations (tic INTEGER,sector INTEGER,status TEXT,result_json TEXT,error TEXT,updated_utc TEXT,PRIMARY KEY(tic,sector))')
    methods=json.dumps(signature(),sort_keys=True)
    previous=conn.execute("SELECT value FROM metadata WHERE key='methods'").fetchone()
    if previous and previous[0]!=methods:
        raise ValueError('Method changed after completed observations; preserve existing run and fork a new database')
    conn.execute("INSERT OR IGNORE INTO metadata VALUES ('methods',?)",(methods,));conn.commit()
    return conn


def record_result(conn,tic,sector,result=None,error=None):
    status=result['status']if result else'error'
    conn.execute('INSERT OR REPLACE INTO observations VALUES (?,?,?,?,?,?)',
                 (tic,sector,status,json.dumps(result)if result else None,str(error)if error else None,utc()))
    conn.commit()


def export(conn):
    OUT.mkdir(parents=True,exist_ok=True)
    targets=[];signals=[];manifest=[]
    for tic,sector,status,raw,error in conn.execute('SELECT tic,sector,status,result_json,error FROM observations ORDER BY sector,tic'):
        if not raw:
            targets.append(dict(tic=tic,sector=sector,status=status,error=error));continue
        result=json.loads(raw);scans=result['scans']
        targets.append(dict(tic=tic,sector=sector,status=status,error='',ra=result['ra'],dec=result['dec'],
                            tmag=result['metadata'].get('TESSMAG'),teff=result['metadata'].get('TEFF'),radius=result['metadata'].get('RADIUS'),
                            n_input=max(s['n_input']for s in scans),n_clean=max(s['n_clean']for s in scans),
                            n_searched=max(s['n_searched']for s in scans),path=result['path'],sha256=result['sha256']))
        manifest.append({k:result[k]for k in ['tic','sector','path','url','sha256','raw_bytes']})
        for scan in scans:
            for rank,peak in enumerate(scan.get('peaks',[]),1):
                reasons=seed_reasons(peak)
                signals.append(dict(tic=tic,sector=sector,ra=result['ra'],dec=result['dec'],width_days=scan['width_days'],label=scan['label'],rank=rank,
                                    period_days=peak['period'],epoch_btjd=peak['t0_btjd'],duration_days=peak['duration_days'],
                                    depth_ppm=peak['depth']*1e6,local_depth_ppm=(peak.get('local_depth')or 0)*1e6,
                                    local_depth_error_ppm=(peak.get('local_depth_err')or 0)*1e6,
                                    snr_diagonal=peak['snr'],local_snr=peak['local_snr'],
                                    n_events=peak['n_events'],n_positive_events=peak['n_positive_events'],
                                    event_snr_max_fraction=peak.get('event_snr_max_fraction'),odd_even_sigma=peak.get('odd_even_sigma'),
                                    seed_retained=not reasons,rejection_reasons=';'.join(reasons),raw_path=result['path'],raw_sha256=result['sha256']))
    for name,rows in [('targets.csv',targets),('signals.csv',signals)]:
        columns=list(dict.fromkeys(k for row in rows for k in row))
        with (OUT/name).open('w',newline='')as stream:
            writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    with gzip.open(OUT/'observations.jsonl.gz','wt')as stream:
        for row in manifest:stream.write(json.dumps(row)+'\n')
    good=[r for r in targets if r['status']=='screened']
    sector_counts={}
    for sector in sorted({r['sector']for r in targets}):
        rows=[r for r in targets if r['sector']==sector]
        sector_counts[str(sector)]=dict(Counter(r['status']for r in rows))
    summary=dict(updated_utc=utc(),methods=signature(),config=CONFIG,widths=WIDTHS,
                 sector_counts=sector_counts,completed_observations=len(targets),
                 actually_period_searched_observations=len(good),unique_period_searched_stars=len({r['tic']for r in good}),
                 stored_signal_rows=len(signals),retained_seed_rows=sum(r['seed_retained']for r in signals),
                 retained_seed_stars=len({r['tic']for r in signals if r['seed_retained']}),
                 raw_files=len(manifest),raw_bytes=sum(r['raw_bytes']for r in manifest),
                 paid_calls=0,discovery_claim=False,
                 limitations='Correlated detrenders and repeated sectors do not add independent stars; retention scores are not false-alarm probabilities.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('EXPORTED',summary['actually_period_searched_observations'],'observations /',summary['unique_period_searched_stars'],'stars;',summary['retained_seed_stars'],'seed stars',flush=True)


def run(limit,workers,downloads,retry_errors):
    selected=list(csv.DictReader((DATA/'selection.csv').open()))
    conn=connect()
    done={(tic,sector)for tic,sector,status in conn.execute('SELECT tic,sector,status FROM observations')if status!='error' or not retry_errors}
    todo=[r for r in selected if (int(r['tic']),int(r['sector']))not in done]
    if limit:todo=todo[:limit]
    # Bounded I/O and CPU futures avoid a large memory queue of raw FITS arrays.
    with ThreadPoolExecutor(downloads)as io,ProcessPoolExecutor(workers)as cpu:
        pending_io={};pending_cpu={};iterator=iter(todo);finished=0
        def fill():
            while len(pending_io)+len(pending_cpu)<downloads+workers*2:
                row=next(iterator,None)
                if row is None:break
                pending_io[io.submit(download,row)]=row
        fill()
        while pending_io or pending_cpu:
            completed,_=wait(list(pending_io)+list(pending_cpu),timeout=30,return_when=FIRST_COMPLETED)
            for future in completed:
                if future in pending_io:
                    row=pending_io.pop(future)
                    try:record=future.result();pending_cpu[cpu.submit(analyze_file,record)]=row
                    except Exception as exc:
                        record_result(conn,int(row['tic']),int(row['sector']),error=exc);finished+=1
                        print('ERROR',row['sector'],row['tic'],str(exc),flush=True)
                else:
                    row=pending_cpu.pop(future)
                    try:record_result(conn,int(row['tic']),int(row['sector']),result=future.result())
                    except Exception as exc:
                        record_result(conn,int(row['tic']),int(row['sector']),error=exc)
                        print('ERROR',row['sector'],row['tic'],str(exc),flush=True)
                    finished+=1
                    if finished%32==0:
                        print('CHECKPOINT',finished,'/',len(todo),flush=True)
                    if finished%128==0:export(conn)
            fill()
    export(conn);conn.close()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['plan','run','export'])
    parser.add_argument('--metadata',type=Path)
    parser.add_argument('--per-sector',type=int,default=1024)
    parser.add_argument('--limit',type=int,default=0)
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--downloads',type=int,default=12)
    parser.add_argument('--retry-errors',action='store_true')
    args=parser.parse_args()
    if args.mode=='plan':print('SELECTED',len(plan(args.metadata,args.per_sector)),'sector observations')
    elif args.mode=='export':
        conn=connect();export(conn);conn.close()
    else:run(args.limit,args.workers,args.downloads,args.retry_errors)
