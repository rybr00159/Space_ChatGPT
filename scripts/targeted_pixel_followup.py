"""Native-pixel localization attempts for remaining deep/sparse new-batch leads."""
from pathlib import Path
import json,os,sys
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import multisector_followup as f
from pooled_south_vet import transit_pixels

def main():
 for tic,sectors in [(183373155,[96,105]),(71274048,[105])]:
  if tic==183373155:
   path=ROOT/'results/multisector_final_pool_followup'/f'tic{tic}.json';d=json.loads(path.read_text());summary=dict(discovery=[dict(sector=r['sector'],source=r['source'])for r in d['observations']],ephemeris=d['ephemeris'])
  else:
   path=f.OUT/f'tic{tic}'/'batch02_followup.json';d=json.loads(path.read_text());summary=d
  for item in summary['discovery']:
   a,meta=f.read_source(item['source']);item['source'].update(ra=meta['metadata']['RA_OBJ'],dec=meta['metadata']['DEC_OBJ'])
  for sector in sectors:
   out=f.OUT/f'tic{tic}'/f'sector{sector}_transit_difference';output=out/'localization.json'
   if output.exists():continue
   result=transit_pixels(tic,summary,sector=sector,size=9)
   f.save(out/'targeted_driver.json',dict(input_path=str(path.relative_to(ROOT)),input_sha256=f.digest(path),script_sha256=f.digest(__file__),warning='Source localization alone does not measure mass or establish a planetary origin; low SNR/bound-width fits are inconclusive.'))
   print('PIXELS_COMPLETE',tic,sector,flush=True)

if __name__=='__main__':main()
