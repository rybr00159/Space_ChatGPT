#!/usr/bin/env python3
"""Consistent SQLite backup and exact offline restore for the new survey.

Freeze may run alongside the survey writer: SQLite backup establishes the
snapshot. Raw observations remain cached and are recoverable from exact URLs
and expected SHA-256, without repeating completed period searches.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,gzip,hashlib,json,sqlite3,sys
from contextlib import closing

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/workspace/tess-research-data/multisector')
OUT=ROOT/'data/multisector/checkpoint'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb')as stream:
        while block:=stream.read(1024*1024):h.update(block)
    return h.hexdigest()


def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    local=CACHE/'checkpoint/survey_snapshot.sqlite';local.parent.mkdir(parents=True,exist_ok=True)
    with closing(sqlite3.connect((CACHE/'survey.sqlite').as_uri()+'?mode=ro',uri=True))as source:
        with closing(sqlite3.connect(local))as destination:
            source.backup(destination)
            # A backup inherits WAL mode. Close/checkpoint the destination
            # before hashing its main file; otherwise later rows can reside
            # solely in the local -wal and be omitted from the gzip archive.
            destination.execute('PRAGMA wal_checkpoint(TRUNCATE)')
            mode=destination.execute('PRAGMA journal_mode=DELETE').fetchone()[0]
            if mode.lower()!='delete':raise ValueError('Snapshot must be self-contained, without WAL')
            destination.commit()
    with closing(sqlite3.connect(local))as conn:
        if conn.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('Snapshot database failed integrity check')
        rows=conn.execute('SELECT tic,sector,status,result_json FROM observations').fetchall()
        methods=json.loads(conn.execute("SELECT value FROM metadata WHERE key='methods'").fetchone()[0])
    compressed=OUT/'survey.sqlite.gz'
    with local.open('rb')as source,compressed.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0,compresslevel=6)as target:
            while block:=source.read(1024*1024):target.write(block)
    manifest=[]
    for tic,sector,status,body in rows:
        if not body:continue
        record=json.loads(body)
        manifest.append({k:record[k]for k in ['tic','sector','path','url','sha256','raw_bytes']})
    raw_manifest=OUT/'raw_observations.jsonl.gz'
    with raw_manifest.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',mtime=0)as zipped:
            for record in manifest:zipped.write((json.dumps(record)+'\n').encode())
    metadata=dict(checked_utc=datetime.now(timezone.utc).isoformat(),methods=methods,
                  database_uncompressed_sha256=sha(local),database_uncompressed_bytes=local.stat().st_size,
                  database_compressed_sha256=sha(compressed),database_compressed_bytes=compressed.stat().st_size,
                  raw_manifest_sha256=sha(raw_manifest),completed_observations=len(rows),
                  searched_observations=sum(r[2]=='screened'for r in rows),
                  searched_unique_stars=len({r[0]for r in rows if r[2]=='screened'}),
                  status_counts={status:sum(r[2]==status for r in rows)for status in sorted({r[2]for r in rows})},
                  raw_files=len(manifest),raw_bytes=sum(r['raw_bytes']for r in manifest),
                  caveat='Snapshot counts can lag a running writer. Storage persistence across cloud tasks is unverified; this tracked backup restores completed period-search results.')
    (OUT/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('CHECKPOINT',metadata['searched_observations'],'searched observations;',metadata['database_compressed_bytes'],'compressed bytes')


def restore(destination):
    m=json.loads((OUT/'manifest.json').read_text());archive=OUT/'survey.sqlite.gz'
    if sha(archive)!=m['database_compressed_sha256']:raise ValueError('Compressed database differs')
    destination=Path(destination)
    if destination.exists():raise FileExistsError('Refusing to overwrite a database: '+str(destination))
    destination.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(archive,'rb')as source,destination.open('wb')as target:
        while block:=source.read(1024*1024):target.write(block)
    if sha(destination)!=m['database_uncompressed_sha256']:raise ValueError('Restored database differs')
    with sqlite3.connect(destination)as conn:
        assert conn.execute('PRAGMA quick_check').fetchone()[0]=='ok'
        assert conn.execute('SELECT count(*) FROM observations').fetchone()[0]==m['completed_observations']
    print('EXACT DATABASE RESTORE VERIFIED',destination,m['completed_observations'],'observations')


def verify():
    m=json.loads((OUT/'manifest.json').read_text())
    for name,digest in m['methods'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Scientific method hash differs: '+name)
    if sha(OUT/'survey.sqlite.gz')!=m['database_compressed_sha256']:raise ValueError('Checkpoint archive differs')
    if sha(OUT/'raw_observations.jsonl.gz')!=m['raw_manifest_sha256']:raise ValueError('Raw observation manifest differs')
    count=0;missing=[]
    with gzip.open(OUT/'raw_observations.jsonl.gz','rt')as stream:
        for line in stream:
            record=json.loads(line);path=Path(record['path'])
            if not path.exists():missing.append(record);continue
            if sha(path)!=record['sha256']:raise ValueError('Original observation differs: '+str(path))
            count+=1
    print('VERIFIED',count,'original observations;',len(missing),'missing; scientific hashes unchanged')
    if missing:raise SystemExit('Restore missing products using raw_observations.jsonl.gz URLs and expected hashes; do not overwrite altered bytes.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['freeze','verify','restore'])
    parser.add_argument('--destination',type=Path);args=parser.parse_args()
    if args.mode=='restore':
        if args.destination is None:parser.error('--destination is required for restore')
        restore(args.destination)
    elif args.mode=='freeze':freeze()
    else:verify()
