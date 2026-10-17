"""Cached, immutable public observation retrieval for TIC214311908."""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,urllib.parse,urllib.request
ROOT=Path(__file__).resolve().parents[1];CACHE=Path('/workspace/tess-research-data/tic214311908');OUT=ROOT/'data/tic214311908';TIC=214311908;RA=260.96969903043;DEC=-51.32383313824

def fetch(name,url,payload=None):
 folder=CACHE/'queries';folder.mkdir(parents=True,exist_ok=True);path=folder/name;manifest=path.with_suffix(path.suffix+'.manifest.json')
 if path.exists():
  m=json.loads(manifest.read_text())
  if hashlib.sha256(path.read_bytes()).hexdigest()!=m['sha256']:raise ValueError('Cached query hash differs')
  if m['url']!=url or m.get('request_body')!=(payload.decode()if payload else None):raise ValueError('Cached query request differs; use a separate snapshot name')
  return path.read_bytes(),m
 request=urllib.request.Request(url,data=payload,headers={'User-Agent':'Space_ChatGPT-public-research/1.0'})
 try:
  with urllib.request.urlopen(request,timeout=90)as response:body=response.read();status=response.status
 except urllib.error.HTTPError as exc:body=exc.read();status=exc.code
 path.write_bytes(body);m=dict(url=url,http_status=status,path=str(path),sha256=hashlib.sha256(body).hexdigest(),bytes=len(body),retrieved_utc=datetime.now(timezone.utc).isoformat(),request_body=payload.decode()if payload else None);manifest.write_text(json.dumps(m,indent=2)+'\n');return body,m

def probe():
 query={'service':'Mast.Caom.Filtered','params':{'columns':'*','filters':[{'paramName':'obs_collection','values':['TESS']},{'paramName':'target_name','values':[str(TIC)]}]},'format':'json','pagesize':10000,'page':1}
 queryurl='https://mast.stsci.edu/api/v0/invoke?'+urllib.parse.urlencode({'request':json.dumps(query)})
 coverage='https://mast.stsci.edu/tesscut/api/v0.1/sector?'+urllib.parse.urlencode(dict(ra=RA,dec=DEC,units='d'))
 with ThreadPoolExecutor(max_workers=2)as pool:
  future=[pool.submit(fetch,'mast_target.json',queryurl),pool.submit(fetch,'tesscut_coverage.json',coverage)]
  results=[f.result()for f in future]
 output=dict(tic=TIC,queries=[r[1]for r in results]);OUT.mkdir(parents=True,exist_ok=True);(OUT/'query_manifest.json').write_text(json.dumps(output,indent=2)+'\n')
 for (body,m)in results:
  print(m['path'],'HTTP',m['http_status'],flush=True)
  try:
   d=json.loads(body)
   if 'results'in d:print('Geometric sectors',[r['sector']for r in d['results']],flush=True)
   elif 'data'in d:
    print('MAST status',d.get('status'),'rows',len(d['data']),flush=True)
    for r in d['data']:print({k:r.get(k)for k in ['obsid','obs_id','provenance_name','sequence_number','t_exptime','dataURL','target_name']},flush=True)
  except Exception as e:print(type(e).__name__,str(e)[:150],flush=True)

def download_qlp():
    expected_path=OUT/'observations.json'
    if not expected_path.exists():raise FileNotFoundError('Restore the tracked original observation manifest first')
    expected=json.loads(expected_path.read_text())
    for row in expected:
        if row['http_status']!=200:continue
        path=Path(row['path']);path.parent.mkdir(parents=True,exist_ok=True)
        if not path.exists():
            request=urllib.request.Request(row['url'],headers={'User-Agent':'Space_ChatGPT-public-research/1.0'})
            with urllib.request.urlopen(request,timeout=120)as response:body=response.read()
            temporary=path.with_suffix(path.suffix+'.retrieved_part');temporary.write_bytes(body)
            if hashlib.sha256(body).hexdigest()!=row['sha256']:raise ValueError('New archive response differs; retained temporary bytes without replacing expected hash')
            temporary.replace(path)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Existing original differs; refusing overwrite')
        print('ORIGINAL_QLP_VERIFIED',row['sector'],path.stat().st_size,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--download-qlp',action='store_true');args=parser.parse_args()
    download_qlp()if args.download_qlp else probe()
