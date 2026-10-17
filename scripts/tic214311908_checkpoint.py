"""Exact offline FITS/query preservation for the independent TIC214311908 review.

Large redundant download ZIPs and the136MBpublished TARSparent table remain
in local storage with immutable hashes/URLs; every scientific TESS FITS used
here is embedded in bounded Git archives. No original file is deleted.
"""
from pathlib import Path
import argparse,gzip,hashlib,io,json,tarfile
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];CACHE=Path('/workspace/tess-research-data/tic214311908');OUT=ROOT/'data/tic214311908/checkpoint';CAND=ROOT/'candidates/tic214311908'

def sha(path):
 d=hashlib.sha256()
 with Path(path).open('rb')as f:
  while b:=f.read(1024*1024):d.update(b)
 return d.hexdigest()

def save(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def recorded_sources():
 expected={};urls={}
 def register(path,digest,url):
  if not isinstance(path,str)or not isinstance(digest,str)or len(digest)!=64:return
  p=Path(path)
  if not p.is_absolute():p=ROOT/p
  p=p.resolve()
  if not p.is_relative_to(CACHE):return
  if p in expected and expected[p]!=digest:raise ValueError('Conflicting original scientific hashes: '+str(p))
  expected[p]=digest
  if url:urls[p]=url
 def walk(d):
  if isinstance(d,list):
   for v in d:walk(v)
  elif isinstance(d,dict):
   url=d.get('url')or d.get('retrieval_url')
   for key,digest in [('path','sha256'),('local_path','sha256'),('fits_path','fits_sha256'),('zip_path','zip_sha256')]:register(d.get(key),d.get(digest),url)
   for v in d.values():
    if isinstance(v,(dict,list)):walk(v)
 files=list((ROOT/'data/tic214311908').rglob('*.json'))+list(CAND.rglob('*.json'))+list(CACHE.rglob('*.manifest.json'))+list(CACHE.rglob('download_manifest.json'))
 for p in files:
  if p.is_relative_to(OUT):continue
  try:walk(json.loads(p.read_text()))
  except json.JSONDecodeError:continue
 return expected,urls

def archive(name,paths):
 target=OUT/name;temporary=target.with_suffix(target.suffix+'.tmp');members=[]
 with temporary.open('wb')as f:
  with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0,compresslevel=6)as z:
   with tarfile.open(fileobj=z,mode='w|')as tar:
    for p in sorted(paths):
     body=p.read_bytes();relative=str(p.relative_to(CACHE));info=tarfile.TarInfo(relative);info.mode=0o644;info.mtime=0;info.size=len(body);tar.addfile(info,io.BytesIO(body));members.append(dict(path=relative,sha256=hashlib.sha256(body).hexdigest(),bytes=len(body)))
 if temporary.stat().st_size>30*1024**2:raise ValueError('Archive exceeds30MiB; preserve and split, do not delete originals')
 temporary.replace(target);return dict(path=name,sha256=sha(target),bytes=target.stat().st_size,members=members)

def freeze():
 OUT.mkdir(parents=True,exist_ok=True);expected,urls=recorded_sources();fits_files=sorted(CACHE.rglob('*.fits'))
 if len(fits_files)!=13:raise ValueError('Expected thirteen original TESSscience FITS; inspect scope before changing this assertion')
 for p in fits_files:
  if p.resolve()not in expected or sha(p)!=expected[p.resolve()]:raise ValueError('Scientific FITS lacks matching previously recorded hash: '+str(p))
 native=[p for p in fits_files if'astrocut'in p.name];static=[p for p in fits_files if p not in native];archives=[archive('static_lightcurves.tar.gz',static)]
 for p in native:
  sector=int(p.name.split('-')[1][1:]);archives.append(archive(f'native_sector{sector}.tar.gz',[p]))
 queries=[p for p in CACHE.rglob('*')if p.is_file()and p.suffix not in {'.fits','.zip','.npz'}and p.stat().st_size<2*1024**2];archives.append(archive('small_query_responses.tar.gz',queries));inventory=[]
 for p in sorted(CACHE.rglob('*')):
  if not p.is_file():continue
  digest=sha(p);prior=expected.get(p.resolve())
  if prior and digest!=prior:raise ValueError('Previously recorded cached bytes changed: '+str(p))
  inventory.append(dict(path=str(p.relative_to(CACHE)),sha256=digest,bytes=p.stat().st_size,recorded_expected_sha256=prior,url=urls.get(p.resolve())))
 raw=OUT/'retained_cache_inventory.json.gz'
 with raw.open('wb')as f:
  with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0)as z:z.write(json.dumps(inventory,indent=2).encode())
 methods=set(ROOT.glob('scripts/tic214311908_*.py'))|set(CAND.rglob('*.py'))
 for name in ['src/fast_screen.py','src/broad_vet.py','src/archive.py','src/multisector_screen.py','src/multisector_combine.py','src/multisector_followup.py','scripts/pooled_south_vet.py']:methods.add(ROOT/name)
 methods|=set((ROOT/'candidates/method_snapshots').glob('*.py'))
 artifacts=[p for p in CAND.rglob('*')if p.is_file()and p.suffix not in {'.pyc','.py'}and 'verification'not in p.name]
 artifacts +=[p for p in(ROOT/'data/tic214311908').rglob('*')if p.is_file()and not p.is_relative_to(OUT)and'verification'not in p.name]
 m=dict(created_utc=datetime.now(timezone.utc).isoformat(),driver_sha256=sha(__file__),archives=archives,original_tess_fits=len(fits_files),native_cubes=len(native),static_lightcurves=len(static),all_scientific_fits_match_original_hashes=True,methods={str(p.relative_to(ROOT)):sha(p)for p in sorted(methods)},artifacts={str(p.relative_to(ROOT)):sha(p)for p in sorted(artifacts)},retained_cache_inventory=dict(path=raw.name,sha256=sha(raw),files=len(inventory),bytes=sum(r['bytes']for r in inventory),files_matching_recorded_input_hash=sum(bool(r['recorded_expected_sha256'])for r in inventory)),scope=__doc__,limitations=['A prior cloud task cache was absent on this instance; do not assume local cache persistence.','Archival/software byte identities do not prove planetary origin, host assignment or novelty.','Original compressed download ZIPs and parent catalogues remain local; exact scientific FITS and compact query responses are embedded in these archives.'])
 save(OUT/'manifest.json',m);print('FROZEN',dict(fits=len(fits_files),native=len(native),static=len(static),archives={r['path']:r['bytes']for r in archives},cache_files=len(inventory),cache_bytes=sum(r['bytes']for r in inventory)),flush=True)

def verify(destination=None,science_only=False):
 m=json.loads((OUT/'manifest.json').read_text())
 if sha(__file__)!=m['driver_sha256']:raise ValueError('Checkpoint driver changed')
 members=0
 for descriptor in m['archives']:
  path=OUT/descriptor['path']
  if sha(path)!=descriptor['sha256']:raise ValueError('Compressed archive differs')
  expected={r['path']:r for r in descriptor['members']}
  with tarfile.open(path,'r:gz')as tar:
   for entry in tar:
    if not entry.isfile()or entry.name not in expected:raise ValueError('Unexpected archive member')
    r=expected.pop(entry.name);body=tar.extractfile(entry).read()
    if len(body)!=r['bytes']or hashlib.sha256(body).hexdigest()!=r['sha256']:raise ValueError('Archived scientific bytes differ')
    if destination:
     target=(destination/entry.name).resolve()
     if not target.is_relative_to(destination.resolve()):raise ValueError('Unsafe archive path')
     if target.exists():
      if sha(target)!=r['sha256']:raise ValueError('Refusing existing mismatching destination')
     else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(body)
     if sha(target)!=r['sha256']:raise ValueError('Offline restoration hash differs')
    members+=1
  if expected:raise ValueError('Missing archive member')
 if sha(OUT/m['retained_cache_inventory']['path'])!=m['retained_cache_inventory']['sha256']:raise ValueError('Raw inventory differs')
 checked_cache_files=0
 if destination is None:
  for kind in ['methods','artifacts']:
   for relative,digest in m[kind].items():
    if sha(ROOT/relative)!=digest:raise ValueError('Current scientific artifact/source differs: '+relative)
  with gzip.open(OUT/m['retained_cache_inventory']['path'],'rt')as f:inventory=json.load(f)
  for r in inventory:
   if science_only and Path(r['path']).suffix!='.fits':continue
   if sha(CACHE/r['path'])!=r['sha256']:raise ValueError('Retained cached file differs')
   checked_cache_files+=1
 result=dict(verified_utc=datetime.now(timezone.utc).isoformat(),archive_members=members,original_tess_fits=m['original_tess_fits'],scientific_methods_verified=len(m['methods'])if destination is None else 0,scientific_artifacts_verified=len(m['artifacts'])if destination is None else 0,retained_cache_files_verified=checked_cache_files,science_only=bool(science_only),offline_destination=str(destination)if destination else None,manifest_sha256=sha(OUT/'manifest.json'))
 save(OUT/('offline_restore_verification.json'if destination else'verification.json'),result);print('VERIFIED',result,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','verify','restore']);p.add_argument('--destination',type=Path);p.add_argument('--science-only',action='store_true',help='Verify sources, artifacts and all thirteen science FITS without requiring redundant cache ZIPs or the parent catalogue.');a=p.parse_args()
 if a.mode=='freeze':freeze()
 else:
  if a.mode=='restore'and a.destination is None:p.error('restore requires --destination')
  verify(a.destination if a.mode=='restore'else None,a.science_only)
