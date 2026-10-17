"""Fixed pooled-ephemeris raw/local/aperture tests; preserve zero-event coverage."""
from pathlib import Path
import os,sys,json,traceback
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import multisector_followup as f
from pooled_south_vet import incomplete_observation
f.OUT=ROOT/'candidates';OUT=f.OUT/'tic214311908'

def main():
 pool=OUT/'pooled_w0.75.json';d=json.loads(pool.read_text());ep=d['peaks'][0];ephem=dict(period=ep['period'],epoch_btjd=ep['t0_btjd'],duration_days=ep['duration_days']);results=[]
 for record in d['sources']:
  try:r=f.vet_observation(record,dict(period_days=ep['period']),fixed_ephemeris=ephem)
  except Exception as exc:
   r=incomplete_observation(record,ephem,exc);r['traceback']=traceback.format_exc();f.save(OUT/f'sector{record["sector"]}_zero_event_fallback.json',r)
  results.append(r)
 f.save(OUT/'raw_local_vetting.json',dict(tic=214311908,ephemeris=ephem,methods={p:f.digest(ROOT/p)for p in ['scripts/tic214311908_vet.py','src/multisector_followup.py','scripts/pooled_south_vet.py']},pooled_input=dict(path=str(pool.relative_to(ROOT)),sha256=f.digest(pool)),observations=results,limitations='Selected pooled ephemeris; diagnostics and local polynomial baselines are conditional. A predicted event in a gap is not a physical null.'))

if __name__=='__main__':main()
