"""Write known transit ephemerides, preserving NASA confirmed values before TOI values."""
import csv,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent
HOSTS=['306996324','79748331','243185500','284441182','235678745','354944123','144401492','261257684','307210830']

def same_period(p,q):
 return abs(float(p)/float(q)-1)<0.001

if __name__=='__main__':
 confirmed=list(csv.DictReader((ROOT/'confirmed_all.csv').open()))
 toi=list(csv.DictReader((ROOT/'toi.csv').open()))
 result=[]
 for tid in HOSTS:
  matched=[r for r in confirmed if r['tic_id']=='TIC '+tid and r['tran_flag']=='1']
  tt=[r for r in toi if r['tid']==tid and r['pl_orbper'] and r['pl_tranmid']]
  for r in matched:
   t=next((t for t in tt if same_period(t['pl_orbper'],r['pl_orbper'])),None)
   duration=r['pl_trandur'] or (t['pl_trandurh'] if t else '')
   result.append({'tic':tid,'name':r['pl_name'],'period':r['pl_orbper'],'epoch_bjd':r['pl_tranmid'],'duration_hours':duration,'source':'NASA pscomppars confirmed transit; '+r['pl_orbper_reflink']})
  for t in tt:
   if not any(same_period(t['pl_orbper'],r['pl_orbper']) for r in matched):
    result.append({'tic':tid,'name':t['toidisplay'],'period':t['pl_orbper'],'epoch_bjd':t['pl_tranmid'],'duration_hours':t['pl_trandurh'],'source':'NASA TOI table, disposition='+t['tfopwg_disp']+'; rowupdate='+t['rowupdate']})
 path=ROOT.parent/'known_planets.csv'
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=['tic','name','period','epoch_bjd','duration_hours','source']);w.writeheader();w.writerows(result)
 path.with_suffix('.csv.manifest.json').write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'row_count':len(result),'parent_files':['catalogs/confirmed_all.csv','catalogs/toi.csv'],'matching_relative_period_tolerance':.001,'epoch_timescale':'Catalog BJD, no offset applied; inspect primary sources for mixed BJD time standards. TESS uses BJD TDB minus 2457000.','sha256':hashlib.sha256(path.read_bytes()).hexdigest()},indent=2)+'\n')
 for r in result:print(r['tic'],r['name'],r['period'],r['epoch_bjd'],r['duration_hours'])
