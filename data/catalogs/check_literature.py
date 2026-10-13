"""Optional Crossref title searches, with irrelevant results filtered and failures recorded.

Broad hyphenated-name bibliographic searches produced unrelated documents during the
pilot. Use title terms plus 'planet'. For the successful main literature inventory,
run check_arxiv.py instead. Neither service alone establishes a novelty claim.
"""
import argparse,csv,json,re
from urllib.parse import urlencode
from pathlib import Path
from check_field_external import fetch
ROOT=Path(__file__).resolve().parent
HOSTS=['TOI-776','TOI-1064','TOI-1468','TOI-1470','TOI-2095','TOI-4342','TOI-1803','TOI-904','L 98-59']

def normalized(s):return re.sub(r'[^a-z0-9]','',s.lower())
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--host',action='append');args=p.parse_args()
 rows=[];errors=[]
 for host in args.host or HOSTS:
  name='crossref_title_'+host.replace(' ','_')+'.json'
  url='https://api.crossref.org/works?'+urlencode({'query.title':host.replace('-',' ')+' planet','rows':8})
  try:
   body=fetch(name,url);items=json.loads(body)['message']['items']
  except Exception as exc:
   errors.append({'host':host,'url':url,'error':str(exc)});print(host,'search incomplete:',str(exc));continue
  matched=0
  for item in items:
   title=' '.join(item.get('title',[]))
   if normalized(host) in normalized(title):
    rows.append({'host':host,'title':title,'doi':item.get('DOI',''),'url':item.get('URL',''),'date':str(item.get('published',{}).get('date-parts','')),'abstract':item.get('abstract',''),'type':item.get('type',''),'score':item.get('score','')});matched+=1
  print(host,matched,'matching title records (finite scope)')
 with (ROOT/'pilot_crossref_title_matches.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=['host','title','doi','url','date','abstract','type','score']);w.writeheader();w.writerows(rows)
 (ROOT/'pilot_crossref_errors.json').write_text(json.dumps(errors,indent=2)+'\n')
