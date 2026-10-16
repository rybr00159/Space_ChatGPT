#!/usr/bin/env python3
"""Freeze supplemental multi-sector science without touching active writers.

Completed pooled scans, consistent SQLite backups and small cached follow-up
responses are archived separately from large raw observations. Expected input
hashes come from existing scientific records, never from replacement values.
"""
from collections import Counter
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path
import argparse
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sqlite3
import tarfile

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/workspace/tess-research-data/multisector')
OUT=ROOT/'data/multisector/supplemental_checkpoint'
MAX_ARCHIVE_BYTES=30*1024**2
DATABASES=['error_recovery.sqlite','b02_096_105/survey.sqlite','b02_096_105/error_recovery.sqlite']
POOL_DIRS=['combined','combined_mixed','b02_096_105/combined','cross_batch_combined']


def utc():return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb')as stream:
        while block:=stream.read(1024*1024):digest.update(block)
    return digest.hexdigest()


def atomic_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


class Inputs:
    def __init__(self):
        self.expected={};self.verified={};self.urls={}
        preservation=ROOT/'data/multisector/diagnostic_preservation.json'
        self.relocated={}
        if preservation.exists():
            for record in json.loads(preservation.read_text())['files']:
                self.relocated[(ROOT/record['original_path']).resolve()]=record

    def check(self,path,digest,label):
        path=Path(path)
        if not path.is_absolute():path=ROOT/path
        path=path.resolve()
        if not path.is_file()and path in self.relocated:
            record=self.relocated[path]
            if record['sha256']!=digest:raise ValueError('Preservation manifest differs from original expected input: '+str(path))
            path=Path(record['retained_path']).resolve()
        if not isinstance(digest,str)or not re.fullmatch('[0-9a-f]{64}',digest):
            raise ValueError('Malformed recorded hash: '+label)
        if path in self.expected and self.expected[path]!=digest:
            raise ValueError('Conflicting existing hashes for '+str(path))
        self.expected[path]=digest
        if path not in self.verified:
            if not path.is_file()or sha(path)!=digest:
                raise ValueError('Recorded scientific input is missing or changed: '+str(path)+' from '+label)
            self.verified[path]=digest
    def walk(self,value,label):
        if isinstance(value,dict):
            for key in ['methods','source_sha256','method_sha256','source_hashes']:
                mapping=value.get(key)
                if isinstance(mapping,dict):
                    for name,digest in mapping.items():
                        if isinstance(name,str)and name.endswith('.py')and isinstance(digest,str):self.check(name,digest,label+':'+key)
            path=value.get('path');digest=value.get('sha256')or value.get('raw_sha256')
            if isinstance(path,str)and isinstance(digest,str)and not path.startswith(('http:','https:')):
                self.check(path,digest,label)
                if value.get('url'):self.urls[Path(path).resolve()]=value['url']
            if isinstance(value.get('raw_path'),str)and isinstance(value.get('raw_sha256'),str):
                self.check(value['raw_path'],value['raw_sha256'],label)
            for key,child in value.items():
                if isinstance(child,(dict,list)):self.walk(child,label+':'+key)
        elif isinstance(value,list):
            for child in value:
                if isinstance(child,(dict,list)):self.walk(child,label)


def backup_database(relative,work,inputs):
    source=CACHE/relative;target=work/relative;target.parent.mkdir(parents=True,exist_ok=True)
    with closing(sqlite3.connect(source.resolve().as_uri()+'?mode=ro',uri=True))as src:
        src.execute('PRAGMA query_only=ON')
        with closing(sqlite3.connect(target))as dst:
            src.backup(dst)
            dst.execute('PRAGMA wal_checkpoint(TRUNCATE)')
            mode=dst.execute('PRAGMA journal_mode=DELETE').fetchone()[0]
            if mode.lower()!='delete':raise ValueError('Backup is not self-contained')
            dst.commit()
    if target.with_name(target.name+'-wal').exists():raise ValueError('Closed backup still has a WAL sidecar')
    with closing(sqlite3.connect(target.resolve().as_uri()+'?mode=ro',uri=True))as conn:
        if conn.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('Backup failed SQLite integrity check')
        rows=conn.execute('SELECT tic,sector,status,result_json FROM observations ORDER BY sector,tic').fetchall()
        methods=json.loads(conn.execute("SELECT value FROM metadata WHERE key='methods'").fetchone()[0])
        inputs.walk(dict(methods=methods),'database:'+relative)
        inputs.walk(methods,'database:'+relative)
        for tic,sector,status,body in rows:
            if body:inputs.walk(json.loads(body),'database:'+relative+':'+str(tic))
    return dict(path=relative,snapshot_path=str(target),bytes=target.stat().st_size,sha256=sha(target),
                completed_observations=len(rows),status_counts=dict(Counter(r[2]for r in rows)),
                searched_observations=sum(r[2]=='screened'for r in rows),
                searched_unique_stars=len({r[0]for r in rows if r[2]=='screened'}),
                methods=methods,sqlite_quick_check='ok',journal_mode='delete',closed_before_hash=True)


def check_original_records(backups):
    byname={r['path']:r for r in backups};checked=0
    for recovery,main in [('error_recovery.sqlite','survey.sqlite'),('b02_096_105/error_recovery.sqlite','b02_096_105/survey.sqlite')]:
        with closing(sqlite3.connect(Path(byname[recovery]['snapshot_path']).as_uri()+'?mode=ro',uri=True))as recovered:
            original_path=Path(byname[main]['snapshot_path'])if main in byname else CACHE/main
            with closing(sqlite3.connect(original_path.resolve().as_uri()+'?mode=ro',uri=True))as original:
                for tic,sector,body in recovered.execute('SELECT tic,sector,result_json FROM observations WHERE result_json IS NOT NULL'):
                    record=json.loads(body);expected=record.get('original_record_sha256')
                    if expected:
                        prior=original.execute('SELECT result_json FROM observations WHERE tic=? AND sector=?',(tic,sector)).fetchone()
                        if not prior or hashlib.sha256(prior[0].encode()).hexdigest()!=expected:
                            raise ValueError('Recovery original-record hash differs: '+str((tic,sector)))
                        checked+=1
    return checked


def write_archive(name,members):
    """Read every immutable member once; hash exactly the bytes added to tar."""
    path=OUT/name;temporary=path.with_suffix(path.suffix+'.tmp');records=[]
    with temporary.open('wb')as file:
        with gzip.GzipFile(filename='',fileobj=file,mode='wb',mtime=0,compresslevel=6)as zipped:
            with tarfile.open(fileobj=zipped,mode='w|')as archive:
                for relative,source,kind in sorted(members):
                    body=Path(source).read_bytes();digest=hashlib.sha256(body).hexdigest()
                    info=tarfile.TarInfo(relative);info.size=len(body);info.mode=0o644;info.mtime=0
                    archive.addfile(info,io.BytesIO(body))
                    records.append(dict(path=relative,bytes=len(body),sha256=digest,kind=kind))
    if temporary.stat().st_size>MAX_ARCHIVE_BYTES:raise ValueError('Archive exceeds configured 30MiB limit: '+str(temporary))
    temporary.replace(path)
    return dict(path=name,bytes=path.stat().st_size,sha256=sha(path),members=records,
                uncompressed_member_bytes=sum(r['bytes']for r in records))


def query_files():
    paths=set()
    allowed={'.json','.txt','.csv','.xml','.html','.bin','.http404','.failed0','.failed1','.failed2'}
    for relative in ['followup','probe','latest_probe','raw','b02_096_105/raw']:
        for path in (CACHE/relative).rglob('*'):
            if path.is_file()and path.suffix.lower()in allowed and path.stat().st_size<=2*1024**2:paths.add(path)
    for relative in ['combined/run_identity.json','b02_096_105/combined/run_identity.json','cross_batch_combined/run_identity.json']:
        path=CACHE/relative
        if path.is_file():paths.add(path)
    return sorted(paths)


def raw_inventory(inputs):
    records=[]
    for relative in ['raw','probe','followup','b02_096_105/raw']:
        for path in sorted((CACHE/relative).rglob('*')):
            if not path.is_file()or path.suffix.lower()not in {'.fits','.zip'}:continue
            path=path.resolve();digest=sha(path)
            expected=inputs.expected.get(path)
            manifest=path.with_suffix('.manifest.json')
            url=inputs.urls.get(path)
            if manifest.is_file():
                m=json.loads(manifest.read_text());prior=m.get('sha256')
                if prior:
                    if prior!=digest:raise ValueError('Raw observation differs from existing download manifest: '+str(path))
                    expected=prior;inputs.expected[path]=prior;inputs.verified[path]=prior
                url=m.get('url')or url
            if expected and expected!=digest:raise ValueError('Raw observation differs from existing scientific record: '+str(path))
            if not url and path.name.startswith('hlsp_qlp_'):
                found=re.search(r'_s(\d+)-(\d{16})_',path.name)
                if found:
                    sector,tic=found.groups();sub='/'.join(tic[i:i+4]for i in range(0,16,4))
                    url=f'https://archive.stsci.edu/hlsps/qlp/s{sector}/{sub}/{path.name}'
            records.append(dict(path=str(path.relative_to(CACHE)),bytes=path.stat().st_size,sha256=digest,
                expected_sha256=expected,validated_against_recorded_input=bool(expected),url=url,
                kind='fits'if path.suffix.lower()=='.fits'else'tesscut_zip'))
    target=OUT/'raw_cache_inventory.jsonl.gz'
    with target.open('wb')as stream:
        with gzip.GzipFile(filename='',fileobj=stream,mode='wb',mtime=0)as zipped:
            for row in records:zipped.write((json.dumps(row,sort_keys=True)+'\n').encode())
    return dict(path=target.name,sha256=sha(target),bytes=target.stat().st_size,files=len(records),
        fits_files=sum(r['kind']=='fits'for r in records),zip_files=sum(r['kind']=='tesscut_zip'for r in records),
        retained_raw_bytes=sum(r['bytes']for r in records),
        fits_bytes=sum(r['bytes']for r in records if r['kind']=='fits'),
        files_with_recorded_expected_hash=sum(r['validated_against_recorded_input']for r in records),
        files_without_public_url=sum(not r['url']for r in records))


def preserve_previous(stamp):
    if not (OUT/'manifest.json').exists():return
    destination=CACHE/'supplemental_checkpoint/archive_history'/stamp;destination.mkdir(parents=True,exist_ok=False)
    for path in OUT.iterdir():
        if path.is_file()and path.suffix!='.tmp':shutil.copy2(path,destination/path.name)


def freeze():
    OUT.mkdir(parents=True,exist_ok=True);stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    preserve_previous(stamp)
    work=CACHE/'supplemental_checkpoint/snapshots'/stamp;work.mkdir(parents=True,exist_ok=False)
    inputs=Inputs();backups=[backup_database(p,work,inputs)for p in DATABASES]
    original_checks=check_original_records(backups)
    pooled=[];pool_counts={}
    for directory in POOL_DIRS:
        paths=sorted((CACHE/directory).glob('tic*.json'));pool_counts[directory]=len(paths)
        for path in paths:
            record=json.loads(path.read_text());inputs.walk(record,'pooled:'+str(path))
            pooled.append((str(path.relative_to(CACHE)),path,'completed_pooled_scans'))
    # Catalog and selection checks reuse the original recorded digests.
    for path in sorted((ROOT/'data/multisector').glob('inventory_s*_manifest.json')):
        inputs.walk(json.loads(path.read_text()),'inventory:'+str(path))
    for path in sorted((ROOT/'data/multisector/catalogs').rglob('*.manifest.json')):
        m=json.loads(path.read_text());body=path.with_name(path.name.removesuffix('.manifest.json'))
        if body.is_file()and m.get('sha256'):inputs.check(body,m['sha256'],'catalog:'+str(path))
    for directory in ['candidates/multisector','results/multisector_final_pool_followup']:
        for path in sorted((ROOT/directory).rglob('*.json')):
            inputs.walk(json.loads(path.read_text()),'scientific_output:'+str(path))
    queries=query_files()
    for path in queries:
        if path.name.endswith('.manifest.json'):
            m=json.loads(path.read_text());body=path.with_name(path.name.removesuffix('.manifest.json')+'.json')
            # Coverage manifests use stem.manifest.json -> stem.json, whereas
            # download manifests identify their raw FITS via an explicit path.
            if body.is_file()and m.get('sha256'):inputs.check(body,m['sha256'],'query:'+str(path))
    raw=raw_inventory(inputs)
    archives=[write_archive('pooled_scans.tar.gz',pooled),
        write_archive('database_snapshots.tar.gz',[(r['path'],Path(r['snapshot_path']),'consistent_sqlite_backup')for r in backups]),
        write_archive('query_responses.tar.gz',[(str(p.relative_to(CACHE)),p,'public_query_or_download_manifest')for p in queries])]
    checks=OUT/'validated_inputs.jsonl.gz'
    with checks.open('wb')as stream:
        with gzip.GzipFile(filename='',fileobj=stream,mode='wb',mtime=0)as zipped:
            for path,digest in sorted(inputs.verified.items()):zipped.write((json.dumps(dict(path=str(path),expected_sha256=digest))+'\n').encode())
    manifest=dict(created_utc=utc(),checkpoint_driver_sha256=sha(Path(__file__)),
        archives=archives,databases=backups,pooled_counts=pool_counts,query_files=len(queries),
        original_recovery_record_hashes_checked=original_checks,
        validated_scientific_inputs=dict(path=checks.name,sha256=sha(checks),files=len(inputs.verified)),
        raw_cache_inventory=raw,all_recorded_input_checks_passed=True,
        scope='Supplemental phase-1 empirical-error DB, named phase-2 DBs, completed pooled JSON and small follow-up responses. Original producer-error DB has its separate main checkpoint.',
        limitations='Different archives and live DB snapshots can have different completion times. Counts can lag running writers; repeated sectors and pooled searches overlap existing stars. External storage persistence is unverified. Raw FITS/zip files remain local with hashes and available public URLs; they are not embedded in these archives. Prepared NPZ files are hash-checked where recorded, but are not duplicated in Git archives.')
    atomic_json(OUT/'manifest.json',manifest)
    print('SUPPLEMENTAL_FROZEN',json.dumps(dict(pooled=pool_counts,databases={r['path']:r['completed_observations']for r in backups},raw=raw,archives={r['path']:r['bytes']for r in archives})),flush=True)


def verify(destination=None,restore=False,check_inputs=True):
    m=json.loads((OUT/'manifest.json').read_text())
    if sha(Path(__file__))!=m['checkpoint_driver_sha256']:raise ValueError('Checkpoint driver differs')
    count=0;dbcounts={}
    for descriptor in m['archives']:
        archive=OUT/descriptor['path']
        if sha(archive)!=descriptor['sha256']:raise ValueError('Compressed archive differs')
        expected={r['path']:r for r in descriptor['members']}
        with tarfile.open(archive,'r:gz')as stream:
            for member in stream:
                if not member.isfile()or member.name not in expected:raise ValueError('Unexpected archive member')
                record=expected.pop(member.name);body=stream.extractfile(member).read()
                if len(body)!=record['bytes']or hashlib.sha256(body).hexdigest()!=record['sha256']:raise ValueError('Archived member differs')
                if destination is not None:
                    target=(destination/member.name).resolve()
                    if not target.is_relative_to(destination.resolve()):raise ValueError('Unsafe member path')
                    if target.exists():
                        if sha(target)!=record['sha256']:raise ValueError('Refusing different existing destination: '+str(target))
                    elif restore:
                        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(body)
                    else:raise ValueError('Verification target missing: '+str(target))
                    if record['kind']=='consistent_sqlite_backup':
                        with closing(sqlite3.connect(target.as_uri()+'?mode=ro',uri=True))as conn:
                            if conn.execute('PRAGMA quick_check').fetchone()[0]!='ok':raise ValueError('Restored DB integrity failure')
                            n=conn.execute('SELECT count(*) FROM observations').fetchone()[0]
                            descriptor_db=next(r for r in m['databases']if r['path']==member.name)
                            if n!=descriptor_db['completed_observations']:raise ValueError('Restored DB row count differs')
                            if conn.execute('PRAGMA journal_mode').fetchone()[0].lower()!='delete':raise ValueError('Restored DB requires WAL')
                            dbcounts[member.name]=n
                count+=1
        if expected:raise ValueError('Archive is incomplete')
    for record in [m['validated_scientific_inputs'],m['raw_cache_inventory']]:
        if sha(OUT/record['path'])!=record['sha256']:raise ValueError('Input/raw manifest differs')
    checked=0
    if check_inputs:
        with gzip.open(OUT/m['validated_scientific_inputs']['path'],'rt')as stream:
            for line in stream:
                record=json.loads(line);path=Path(record['path'])
                if not path.is_file()or sha(path)!=record['expected_sha256']:raise ValueError('Recorded scientific input differs: '+str(path))
                checked+=1
    result=dict(verified_utc=utc(),archive_members_verified=count,recorded_inputs_verified=checked,
                restored_database_rows=dbcounts,offline_destination=str(destination)if destination else None,
                original_manifest_sha256=sha(OUT/'manifest.json'),external_persistence_verified=False)
    atomic_json(OUT/('offline_restore_verification.json'if restore else'verification.json'),result)
    print('SUPPLEMENTAL_VERIFIED',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['freeze','verify','restore'])
    parser.add_argument('--destination',type=Path);parser.add_argument('--archive-only',action='store_true',help='Skip local original-input checks; useful for an offline restored workspace.')
    args=parser.parse_args()
    if args.mode=='freeze':freeze()
    else:
        if args.mode=='restore'and args.destination is None:parser.error('restore requires --destination')
        verify(args.destination,args.mode=='restore',not args.archive_only)
