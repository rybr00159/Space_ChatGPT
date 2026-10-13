"""Retrieve public Crossref metadata for initial host recognition checks (not exhaustive)."""
from concurrent.futures import ThreadPoolExecutor
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from datetime import datetime,timezone
from pathlib import Path
import csv,json,hashlib
ROOT=Path(__file__).resolve().parent
HOSTS=['TOI-776','TOI-1064','TOI-1468','TOI-1470','TOI-2095','TOI-4342','TOI-1803','TOI-904','L 98-59']

def get(host):
 name='literature_'+host.replace(' ','_')+'.json';path=ROOT/name
 url='https://api.crossref.org/works?'+urlencode({'query.bibliographic':host,'rows':8})
 if path.exists(): body=path.read_bytes()
 else:
  with urlopen(Request(url,headers={'User-Agent':'Space_ChatGPT-public-research/0.1'}),timeout=60) as r:body=r.read()
  path.write_bytes(body)
  path.with_suffix('.json.manifest.json').write_text(json.dumps({'retrieved_utc':datetime.now(timezone.utc).isoformat(),'url':url,'sha256':hashlib.sha256(body).hexdigest(),'notes':'Public bibliographic search by host label, not an exhaustive novelty search.'},indent=2)+'\n')
 data=json.loads(body)
 items=data['message']['items'];out=[]
 for i in items:
  title=' '.join(i.get('title',[]))
  if host.lower().replace(' ','') in title.lower().replace(' ',''):
   out.append({'host':host,'title':title,'doi':i.get('DOI',''),'url':i.get('URL',''),'date':str(i.get('published',{}).get('date-parts','')),'abstract':i.get('abstract',''),'type':i.get('type',''),'score':i.get('score','')})
 return host,out
if __name__=='__main__':
 rows=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for h,result in pool.map(get,HOSTS):
   rows.extend(result);print(h,[(r['title'],r['doi']) for r in result])
 p=ROOT/'pilot_literature.csv'
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=['host','title','doi','url','date','abstract','type','score']);w.writeheader();w.writerows(rows)
