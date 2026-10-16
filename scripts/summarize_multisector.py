"""Read-only, deduplicated accounting for the completed six-sector search."""
from pathlib import Path
import csv,gzip,hashlib,json,sqlite3
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];CACHE=Path('/workspace/tess-research-data/multisector');OUT=ROOT/'results/multisector'
BATCHES=[('b01','survey.sqlite','error_recovery.sqlite'),('b02_096_105','b02_096_105/survey.sqlite','b02_096_105/error_recovery.sqlite')]
TABLES=[('b01_producer','results/multisector/signals.csv.gz'),('b01_empirical','results/multisector_error_recovery/signals.csv'),('b02_producer','results/multisector_batches/b02_096_105/signals.csv'),('b02_empirical','results/multisector_batches/b02_096_105/error_recovery/signals.csv'),('b01_producer_pooled','results/multisector/combined_signals.csv.gz'),('b01_mixed_pooled','results/multisector_mixed_pool/signals.csv.gz'),('b02_pooled','results/multisector_batches/b02_096_105/pooled/signals.csv.gz'),('cross_batch_pooled','results/multisector_cross_batch/signals.csv.gz')]

def read(db):
 with sqlite3.connect((CACHE/db).as_uri()+'?mode=ro',uri=True)as c:
  rows=c.execute('select tic,sector,status,result_json from observations order by sector,tic').fetchall()
 return {(tic,sector):(status,json.loads(body)if body else {})for tic,sector,status,body in rows}

def main():
 rows=[];allkeys=set();batchsets={};sectors={}
 for batch,main,recovery in BATCHES:
  a,b=read(main),read(recovery);batchsets[batch]=set()
  for (tic,sector),(status,r)in a.items():
   if (tic,sector)in allkeys:raise ValueError('Overlapping target-sector observations between distinct batches')
   allkeys.add((tic,sector));rs,rr=b.get((tic,sector),(None,{}));searched=status=='screened'or rs=='screened';assumption='producer_errors'if status=='screened'else'empirical_adjacent_difference_MAD'if rs=='screened'else'not_searched'
   row=dict(batch=batch,tic=tic,sector=sector,period_searched=searched,producer_status=status,empirical_status=rs,error_assumption=assumption,ra=r.get('ra'),dec=r.get('dec'),tmag=r.get('tmag'),teff=r.get('teff'),radius_solar=r.get('radius'),raw_url=r.get('url'),raw_path=r.get('path'),raw_sha256=r.get('sha256'),raw_bytes=r.get('raw_bytes'))
   rows.append(row)
   if searched:batchsets[batch].add(tic)
   z=sectors.setdefault(sector,dict(selected=0,producer_searched=0,empirical_searched=0,unusable=0));z['selected']+=1;z['producer_searched']+=int(status=='screened');z['empirical_searched']+=int(rs=='screened');z['unusable']+=int(not searched)
 target=OUT/'all_observations.csv'
 with target.open('w',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 old=set()
 with gzip.open(ROOT/'results/broad/targets.csv.gz','rt')as f:
  reader=csv.DictReader(f)
  for row in reader:
   if row.get('status')=='screened':old.add(int(row['tic']))
 searched={r['tic']for r in rows if r['period_searched']};summaries=[]
 for label,name in TABLES:
  path=ROOT/name;opener=gzip.open if name.endswith('.gz')else open
  with opener(path,'rt',newline='')as f:count=sum(1 for _ in csv.DictReader(f))
  summaries.append(dict(label=label,path=name,rows=count,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),counting_unit='correlated_peak_rows_not_discoveries'))
 result=dict(updated_utc=datetime.now(timezone.utc).isoformat(),selected_observations=len(rows),selected_unique_tics=len({r['tic']for r in rows}),period_searched_observations=sum(r['period_searched']for r in rows),period_searched_unique_tics=len(searched),unusable_observations=sum(not r['period_searched']for r in rows),sectors=sectors,batches={k:dict(searched_unique_tics=len(v))for k,v in batchsets.items()},shared_searched_tics=len(set.intersection(*batchsets.values())),old_sector1_searched_tics=len(old),overlap_with_old_sector1=len(searched&old),outside_old_sector1=len(searched-old),search_tables=summaries,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),notes=['All four databases are read-only. Each target-sector observation appears once, with the successful uncertainty assumption explicit.','Follow-up observations, control stars and repeated pooled searches are not counted as additional survey stars.','The search samples part of six sectors; it is not a complete census of any sector.','Stored catalogue parameters are selection estimates, not precise stellar characterization.'])
 assert result['period_searched_observations']==5058 and result['period_searched_unique_tics']==4391
 (OUT/'campaign_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in ['search_tables','notes']},indent=2))

if __name__=='__main__':main()
