#!/usr/bin/env python3
"""Preserve exact query snapshots and verify this research phase's checkpoint.

Large original images, spectra and all intermediate plots remain in the cache.
Small HTTP query snapshots are additionally archived in Git for exact offline
restoration. Missing static FITS/NPZ products can be retrieved from their own
manifests; changed bytes are preserved separately and never silently accepted.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,tarfile
import requests

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/workspace/tess-research-data/new_stars')
OUT=ROOT/'data/new_stars/checkpoint'
METHODS=['src/uncatalogued_imaging.py','src/new_companion_search.py','src/infrared_star_scout.py',
 'src/spectroscopic_companion_scout.py','src/rv_spectral_validation.py','src/rv_photometric_veto.py',
 'scripts/rv_lead_recognition.py','scripts/run_uncatalogued_imaging.py','scripts/summarize_new_stars.py',
 'scripts/new_star_checkpoint.py']


def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as stream:
        while block:=stream.read(1024*1024):h.update(block)
    return h.hexdigest()


def verify_input_manifests():
    checked={}
    def check(value,context):
        if isinstance(value,dict):
            digest=value.get('sha256')
            loc=next((value[k]for k in ['retained_path','local_path','raw_path','path','file']if value.get(k)),None)
            if digest and loc:
                path=Path(loc)
                if not path.is_absolute():
                    path=ROOT/path if (ROOT/path).exists() else context.parent/path
                if not path.exists():raise FileNotFoundError('Manifest input missing: '+str(path))
                if path not in checked:checked[path]=sha(path)
                if checked[path]!=digest:raise ValueError('Expected scientific input hash mismatch: '+str(path))
            for child in value.values():check(child,context)
        elif isinstance(value,list):
            for child in value:check(child,context)
    for context in (ROOT/'data/new_stars').rglob('*.json'):
        if context.is_relative_to(OUT):continue
        value=json.loads(context.read_text())
        if context.name.endswith('.manifest.json')and isinstance(value,dict)and value.get('sha256'):
            adjacent=context.with_name(context.name.removesuffix('.manifest.json'))
            if adjacent.exists():check(dict(path=str(adjacent),sha256=value['sha256']),context)
        check(value,context)
    print('VERIFIED expected scientific-input manifests:',len(checked),'unique retained files')
    return len(checked)


def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    expected_inputs=verify_input_manifests()
    inventory=[dict(path=str(p.relative_to(CACHE)),bytes=p.stat().st_size,sha256=sha(p))
               for p in sorted(CACHE.rglob('*'))if p.is_file()]
    (OUT/'cache_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    small=[p for p in sorted(CACHE.rglob('*'))if p.is_file()and p.suffix in ['.json','.csv','.txt','.html','.xml']]
    archive=OUT/'query_snapshots.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=6)as tar:
        for p in small:tar.add(p,arcname=str(p.relative_to(CACHE)),recursive=False)
    archived=dict(path=str(archive.relative_to(ROOT)),sha256=sha(archive),bytes=archive.stat().st_size,
                  files=len(small),uncompressed_bytes=sum(p.stat().st_size for p in small))
    (OUT/'query_archive.json').write_text(json.dumps(archived,indent=2)+'\n')
    tracked=[]
    for directory in [ROOT/'data/new_stars',ROOT/'candidates/new_stars',ROOT/'results/new_stars']:
        for p in sorted(directory.rglob('*')):
            if not p.is_file()or p.name=='verification.json':continue
            tracked.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
    report=dict(checked_utc=datetime.now(timezone.utc).isoformat(),methods={n:sha(ROOT/n)for n in METHODS},
                tracked=tracked,cache_files=len(inventory),cache_bytes=sum(v['bytes']for v in inventory),expected_inputs=expected_inputs,
                caveat='Byte identity is not scientific validation, discovery priority or proof that this cloud cache persists into a fresh task.')
    (OUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('FROZEN',len(tracked),'tracked artifacts;',len(inventory),'cache files;',archived['bytes'],'query-archive bytes')


def restore_queries():
    meta=json.loads((OUT/'query_archive.json').read_text());archive=ROOT/meta['path']
    assert sha(archive)==meta['sha256']
    count=0
    with tarfile.open(archive,'r:gz')as tar:
        for member in tar:
            if not member.isfile():raise ValueError('Unexpected archive member')
            target=(CACHE/member.name).resolve()
            if not target.is_relative_to(CACHE):raise ValueError('Unsafe archive path')
            body=tar.extractfile(member).read()
            if target.exists():
                if hashlib.sha256(body).hexdigest()!=sha(target):raise ValueError('Existing cache differs: '+str(target))
                continue
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(body);count+=1
    print('RESTORED',count,'exact query snapshot files')


def verify():
    report=json.loads((OUT/'verification.json').read_text())
    for name,digest in report['methods'].items():assert sha(ROOT/name)==digest,name
    for r in report['tracked']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    missing=[]
    for r in json.loads((OUT/'cache_inventory.json').read_text()):
        path=CACHE/r['path']
        if not path.exists():missing.append(r['path']);continue
        assert sha(path)==r['sha256'],str(path)
    print('VERIFIED',len(report['methods']),'methods,',len(report['tracked']),'tracked artifacts;',
          report['cache_files']-len(missing),'cache files; missing',len(missing))
    if missing:
        print('Restore query snapshots first; static observations use --download --max-downloads N.')
        raise SystemExit(2)


def download(limit):
    """Read-only retrieval of missing, explicitly identified static products."""
    candidates=[]
    # Reconstructed image metadata live next to their missing raw product.
    for p in CACHE.rglob('*.json'):
        try:m=json.loads(p.read_text())
        except (ValueError,OSError):continue
        if not isinstance(m,dict):continue
        target=Path(m.get('path',str(p.with_suffix(''))))
        if m.get('status',m.get('http_status'))==200 and m.get('url')and m.get('sha256'):
            candidates.append((target,m))
    def collect(value):
        if isinstance(value,dict):
            path=value.get('path')or value.get('raw_path')
            if path and value.get('sha256')and value.get('url')and value.get('status',value.get('http_status'))==200:
                candidates.append((Path(path),value))
            for child in value.values():collect(child)
        elif isinstance(value,list):
            for child in value:collect(child)
    for p in (ROOT/'data/new_stars').rglob('*.json'):
        if p.is_relative_to(OUT):continue
        try:collect(json.loads(p.read_text()))
        except ValueError:continue
    done=set();n=0
    for path,record in candidates:
        if path in done or path.exists():continue
        done.add(path)
        if not path.resolve().is_relative_to(CACHE):continue
        if n>=limit:break
        r=requests.get(record['url'],timeout=120);r.raise_for_status();body=r.content
        path.parent.mkdir(parents=True,exist_ok=True)
        if hashlib.sha256(body).hexdigest()!=record['sha256']:
            altered=path.with_name(path.name+'.unexpected_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S'))
            altered.write_bytes(body);raise ValueError('Changed product preserved separately: '+str(altered))
        path.write_bytes(body);n+=1;print('RESTORED',path.name,flush=True)
    print('Downloaded missing static products',n)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['freeze','verify','restore-queries','download'])
    parser.add_argument('--max-downloads',type=int,default=20);a=parser.parse_args()
    if a.mode=='download':download(a.max_downloads)
    else:globals()[a.mode.replace('-','_')]()
