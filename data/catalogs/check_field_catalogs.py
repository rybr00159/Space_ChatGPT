"""Match the independent pilot against cached public planet catalogs by TIC and sky position."""
import csv,json,math,hashlib
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent

def sep_arcsec(ra1,dec1,ra2,dec2):
 ra1,dec1,ra2,dec2=map(math.radians,[ra1,dec1,ra2,dec2])
 val=math.sin((dec2-dec1)/2)**2+math.cos(dec1)*math.cos(dec2)*math.sin((ra2-ra1)/2)**2
 return math.degrees(2*math.asin(min(1,math.sqrt(val))))*3600

if __name__=='__main__':
 field=json.loads((ROOT.parent/'field_sample.json').read_text())['selected']
 catalogspec=[('confirmed_all.csv','tic_id','ra','dec','pl_name','pl_orbper'),('toi.csv','tid','ra','dec','toidisplay','pl_orbper'),('ctoi.csv','TIC ID','RA','Dec','CTOI','Period (days)')]
 nearest=[];matches=[]
 for name,tidcol,racol,deccol,namecol,pcol in catalogspec:
  catalog=list(csv.DictReader((ROOT/name).open()))
  for star in field:
   tid=str(star['tic']);distances=[]
   for r in catalog:
    try:d=sep_arcsec(star['ra'],star['dec'],float(r[racol]),float(r[deccol]))
    except (ValueError,TypeError):continue
    row={'tic':tid,'catalog':name,'matched_name':r[namecol],'matched_tic':r[tidcol].replace('TIC ',''),'period_days':r[pcol],'separation_arcsec':d,'identifier_match':r[tidcol].replace('TIC ','')==tid,'within_30_arcsec':d<=30}
    distances.append(row)
    if row['identifier_match'] or row['within_30_arcsec']:matches.append(row)
   if distances:nearest.append(min(distances,key=lambda r:r['separation_arcsec']))
 def save(name,rows):
  keys=['tic','catalog','matched_name','matched_tic','period_days','separation_arcsec','identifier_match','within_30_arcsec']
  with (ROOT/name).open('w') as f:
   w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 save('field_planet_catalog_matches.csv',matches);save('field_nearest_catalog_entries.csv',nearest)
 meta={'checked_utc':datetime.now(timezone.utc).isoformat(),'position_radius_arcsec':30,'field_definition':'data/field_sample.json','row_count_matches':len(matches),'rows_nearest':len(nearest),'catalog_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name,*_ in catalogspec},'limitations':'Only cached catalogs and 30-arcsec coordinate criterion. Nonmatches do not establish novelty or lack of prior analysis.'}
 (ROOT/'field_catalog_check.manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
 print('matches',len(matches))
 for r in nearest:print(r['tic'],r['catalog'],r['matched_name'],round(r['separation_arcsec'],1),'arcsec')
