"""Bounded restoration of missing static QLP FITS from recorded public URLs.

Existing files are hash checked and never replaced. Mismatching fresh responses
are retained separately and never accepted by changing an expected hash. Native
TESSCut containers may change on regeneration; use the important-candidate exact
archive for TIC229146531. Other native cube URLs/FFI identities remain recorded.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,re
from datetime import datetime,timezone
import requests
ROOT=Path(__file__).resolve().parents[1];CACHE=Path('/workspace/tess-research-data/multisector');INVENTORY=ROOT/'data/multisector/supplemental_checkpoint/raw_cache_inventory.jsonl.gz'

def main():
 p=argparse.ArgumentParser();p.add_argument('--tic',type=int);p.add_argument('--max-downloads',type=int,default=20);a=p.parse_args();done=0;present=0;errors=[]
 if a.max_downloads<0:p.error('nonnegative max-downloads required')
 manifest=json.loads((INVENTORY.parent/'manifest.json').read_text())['raw_cache_inventory']
 if hashlib.sha256(INVENTORY.read_bytes()).hexdigest()!=manifest['sha256']:raise ValueError('Raw inventory changed')
 with gzip.open(INVENTORY,'rt')as f:
  for line in f:
   r=json.loads(line);path=CACHE/r['path']
   if r['kind']!='fits'or not path.name.startswith('hlsp_qlp_'):continue
   match=re.search(r'_s\d+-(\d{16})_',path.name)
   if not match or(a.tic is not None and int(match.group(1))!=a.tic):continue
   if path.exists():
    if hashlib.sha256(path.read_bytes()).hexdigest()!=r['sha256']:raise ValueError('Existing raw FITS differs: '+str(path))
    present+=1;continue
   if done>=a.max_downloads:continue
   response=requests.get(r['url'],timeout=(20,180));body=response.content;digest=hashlib.sha256(body).hexdigest();path.parent.mkdir(parents=True,exist_ok=True);done+=1
   if response.status_code!=200 or digest!=r['sha256']:
    retained=path.with_name(path.name+'.retrieval_unmatched_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));retained.write_bytes(body);errors.append(dict(url=r['url'],http_status=response.status_code,expected_sha256=r['sha256'],retrieved_sha256=digest,response_path=str(retained)));continue
   path.write_bytes(body)
 print(json.dumps(dict(existing_verified=present,network_requests=done,errors=errors),indent=2))
 if errors:raise SystemExit(1)

if __name__=='__main__':main()
