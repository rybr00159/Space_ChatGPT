"""Reproduce the initial pilot sample from cached NASA Exoplanet Archive CSVs."""
import csv,collections,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent
HOSTS=['TOI-776','TOI-1064','TOI-1468','TOI-1470','TOI-2095','TOI-4342','TOI-1803','TOI-904']

def write_table(name,rows,fields,parents):
 p=ROOT/name
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 m={'generated_utc':datetime.now(timezone.utc).isoformat(),'row_count':len(rows),'columns':fields,'parent_files':parents,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 p.with_suffix('.csv.manifest.json').write_text(json.dumps(m,indent=2)+'\n')

if __name__=='__main__':
 confirmed=list(csv.DictReader((ROOT/'confirmed_transiting.csv').open()))
 toi=list(csv.DictReader((ROOT/'toi.csv').open()))
 planets=[r for h in HOSTS for r in confirmed if r['hostname']==h]
 write_table('pilot_known_planets.csv',planets,list(confirmed[0]),['confirmed_transiting.csv'])
 fields=['hostname','tic_id','ra','dec','st_teff','st_rad','st_mass','sy_tmag','sy_vmag','sy_dist','sy_pnum','known_transiting_planets','known_periods_days','known_toi_labels','known_toi_periods_days','toi_sectors']
 rows=[];matched=[]
 for h in HOSTS:
  rr=[r for r in planets if r['hostname']==h];r=rr[0];tid=r['tic_id'].replace('TIC ','')
  matches=[t for t in toi if t['tid']==tid];matched.extend(matches)
  row={k:r[k] for k in fields[:11]}
  row.update(tic_id=tid,known_transiting_planets=';'.join(t['pl_name'] for t in rr),known_periods_days=';'.join(t['pl_orbper'] for t in rr),known_toi_labels=';'.join(t['toidisplay'] for t in matches),known_toi_periods_days=';'.join(t['pl_orbper'] for t in matches),toi_sectors=';'.join(t['sectors'] for t in matches))
  rows.append(row)
 write_table('pilot_targets.csv',rows,fields,['confirmed_transiting.csv','toi.csv'])
 write_table('pilot_toi_matches.csv',matched,list(toi[0]),['toi.csv'])
 for r in rows:print(r['hostname'],r['tic_id'],r['known_periods_days'],r['known_toi_periods_days'])

 ctoi=list(csv.DictReader((ROOT/'ctoi.csv').open()))
 tids={r['tic_id'] for r in rows}
 matched_ctoi=[r for r in ctoi if r['TIC ID'] in tids]
 write_table('pilot_ctoi_matches.csv',matched_ctoi,list(ctoi[0]),['ctoi.csv','pilot_targets.csv'])
