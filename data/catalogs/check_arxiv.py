"""Reproducible public arXiv host-label search; metadata is not an exhaustive literature review."""
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from pathlib import Path
from datetime import datetime,timezone
import xml.etree.ElementTree as ET
import json,hashlib,csv,re
ROOT=Path(__file__).resolve().parent
HOSTS=['TOI-776','TOI-1064','TOI-1468','TOI-1470','TOI-2095','TOI-4342','TOI-1803','TOI-904','L 98-59']
NS={'a':'http://www.w3.org/2005/Atom','op':'http://a9.com/-/spec/opensearch/1.1/'}
if __name__=='__main__':
 query=' OR '.join('all:"'+h+'"' for h in HOSTS)
 url='https://export.arxiv.org/api/query?'+urlencode({'search_query':query,'max_results':100,'sortBy':'submittedDate','sortOrder':'descending'})
 p=ROOT/'pilot_arxiv.xml'
 if p.exists():body=p.read_bytes()
 else:
  with urlopen(Request(url,headers={'User-Agent':'Space_ChatGPT-public-research/0.1'}),timeout=120) as r:body=r.read()
  p.write_bytes(body)
  p.with_suffix('.xml.manifest.json').write_text(json.dumps({'retrieved_utc':datetime.now(timezone.utc).isoformat(),'url':url,'query':query,'sha256':hashlib.sha256(body).hexdigest(),'notes':'Host-label search; abstracts only. This does not cover all aliases/full text or journals without an arXiv deposit.'},indent=2)+'\n')
 root=ET.fromstring(body);out=[]
 for e in root.findall('a:entry',NS):
  title=' '.join(e.findtext('a:title','',NS).split());summary=' '.join(e.findtext('a:summary','',NS).split())
  text=title+' '+summary
  hosts=[]
  for h in HOSTS:
   pattern=re.escape(h).replace(r'\-',r'[-\s]?').replace(r'\ ',r'\s*')
   if re.search(pattern,text,re.I):hosts.append(h)
  out.append({'hosts':';'.join(hosts),'title':title,'arxiv_id':e.findtext('a:id','',NS),'published':e.findtext('a:published','',NS),'updated':e.findtext('a:updated','',NS),'abstract':summary})
 with (ROOT/'pilot_arxiv.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=['hosts','title','arxiv_id','published','updated','abstract']);w.writeheader();w.writerows(out)
 print('reported total',root.findtext('op:totalResults','?',NS),'returned',len(out))
 for r in out:print(r['hosts'],r['arxiv_id'],r['title'])
