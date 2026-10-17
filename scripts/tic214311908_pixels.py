"""Native-pixel conditional localization of independently selected deep dips."""
from pathlib import Path
import os,sys,json
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import multisector_followup as f
from pooled_south_vet import transit_pixels
f.OUT=ROOT/'candidates';f.CACHE=Path('/workspace/tess-research-data/tic214311908')

def main():
 p=ROOT/'candidates/tic214311908/pooled_w0.75.json';d=json.loads(p.read_text());peak=d['peaks'][0];ephem=dict(period=peak['period'],epoch_btjd=peak['t0_btjd'],duration_days=peak['duration_days']);records=[]
 for r in d['sources']:
  a,m=f.read_source(r);records.append(dict(sector=r['sector'],source=dict(r,ra=m['metadata']['RA_OBJ'],dec=m['metadata']['DEC_OBJ'])))
 summary=dict(ephemeris=ephem,discovery=records)
 for sector in [66,93]:
  out=f.OUT/'tic214311908'/f'sector{sector}_transit_difference'
  if (out/'localization.json').exists():continue
  transit_pixels(214311908,summary,sector=sector,size=9)
  result=json.loads((out/'localization.json').read_text())
  f.save(out/'driver.json',dict(script_sha256=f.digest(__file__),pooled_input=dict(path=str(p.relative_to(ROOT)),sha256=f.digest(p)),limitations='Deep dips were selected by BLS. Difference images use the same photons and selected template; not independent detection significance or precise unresolved-source assignment.'))
  print('NATIVE_COMPLETED',sector,[(r['max_difference_snr'],r['offset_target_pixels'],r['localization_valid'])for r in result['tests']],flush=True)

if __name__=='__main__':main()
