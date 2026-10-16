"""Exact offline preservation of TIC229146531's seven original scientific FITS.

The primary target's dynamic TESSCut container must not depend solely on a
future API response reproducing identical headers. This small archive retains
the actual native cube and six original QLP observations, without redundant ZIPs.
"""
from pathlib import Path
import argparse,gzip,hashlib,io,json,tarfile
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];CACHE=Path('/workspace/tess-research-data/multisector');OUT=ROOT/'data/multisector/candidate_checkpoint';CAND=ROOT/'candidates/multisector/tic229146531'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def records():
 d=json.loads((CAND/'batch02_followup.json').read_text());rows=[r['source']for group in ['discovery','other_batch','additional']for r in d[group]];rows.append(json.loads((CAND/'sector96_transit_difference/raw_manifest.json').read_text()));return rows

def freeze():
 OUT.mkdir(parents=True,exist_ok=True);members=[];archive=OUT/'tic229146531_original_fits.tar.gz';temp=archive.with_suffix('.tmp')
 with temp.open('wb')as stream:
  with gzip.GzipFile(filename='',fileobj=stream,mode='wb',mtime=0,compresslevel=6)as zipped:
   with tarfile.open(fileobj=zipped,mode='w|')as tf:
    for r in records():
     path=Path(r['path']);expected=r['sha256']
     if sha(path)!=expected:raise ValueError('Original FITS hash differs')
     body=path.read_bytes();name=str(path.relative_to(CACHE));info=tarfile.TarInfo(name);info.size=len(body);info.mtime=0;info.mode=0o644;tf.addfile(info,io.BytesIO(body));members.append(dict(path=name,sha256=expected,bytes=len(body),url=r['url']))
 if temp.stat().st_size>30*1024**2:raise ValueError('Important-candidate archive exceeds30MiB')
 temp.replace(archive);m=dict(created_utc=datetime.now(timezone.utc).isoformat(),script_sha256=sha(__file__),archive=archive.name,sha256=sha(archive),bytes=archive.stat().st_size,members=members,scope=__doc__);(OUT/'manifest.json').write_text(json.dumps(m,indent=2)+'\n');print('CANDIDATE_FROZEN',len(members),'FITS',m['bytes'],'bytes')

def verify(destination=None):
 m=json.loads((OUT/'manifest.json').read_text());archive=OUT/m['archive']
 if sha(__file__)!=m['script_sha256']or sha(archive)!=m['sha256']:raise ValueError('Candidate checkpoint changed')
 expected={r['path']:r for r in m['members']};count=0
 with tarfile.open(archive,'r:gz')as tf:
  for member in tf:
   if member.name not in expected or not member.isfile():raise ValueError('Unexpected archive member')
   r=expected.pop(member.name);body=tf.extractfile(member).read()
   if len(body)!=r['bytes']or hashlib.sha256(body).hexdigest()!=r['sha256']:raise ValueError('Archived scientific bytes differ')
   if destination:
    target=(destination/member.name).resolve()
    if not target.is_relative_to(destination.resolve()):raise ValueError('Unsafe member path')
    if target.exists():
     if sha(target)!=r['sha256']:raise ValueError('Refusing existing different file')
    else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(body)
    if sha(target)!=r['sha256']:raise ValueError('Restore identity failed')
   count+=1
 if expected:raise ValueError('Missing archive members')
 result=dict(verified_utc=datetime.now(timezone.utc).isoformat(),exact_original_fits=count,destination=str(destination)if destination else None,manifest_sha256=sha(OUT/'manifest.json'))
 (OUT/('offline_restore_verification.json'if destination else'verification.json')).write_text(json.dumps(result,indent=2)+'\n');print('CANDIDATE_VERIFIED',result)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','verify','restore']);p.add_argument('--destination',type=Path);a=p.parse_args()
 if a.mode=='freeze':freeze()
 else:
  if a.mode=='restore'and a.destination is None:p.error('restore requires --destination')
  verify(a.destination if a.mode=='restore'else None)
