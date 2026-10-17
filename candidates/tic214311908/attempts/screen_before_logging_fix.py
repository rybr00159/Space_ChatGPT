"""Independent multi-sector transit search for user-nominated TIC214311908.

All good-quality finite producer-error photometry is retained. Screening scores
are conditional and uncalibrated for correlated noise or trials. Data-product
and source checks are required before a planetary or novelty interpretation.
"""
from pathlib import Path
import argparse,json,os,sys,hashlib
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from multisector_screen import read_curve,screen_arrays
from multisector_combine import search_prepared
from fast_screen import preprocess_curve,screen_curve,event_stats
import multisector_followup as f
OUT=ROOT/'candidates/tic214311908';CACHE=Path('/workspace/tess-research-data/tic214311908')

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--pool',action='store_true');args=parser.parse_args();records=json.loads((ROOT/'data/tic214311908/observations.json').read_text());sources=[]
 methods={p:f.digest(ROOT/p)for p in ['scripts/tic214311908_screen.py','src/multisector_screen.py','src/multisector_combine.py','src/fast_screen.py']}
 for record in records:
  if record['http_status']!=200:continue
  curve=read_curve(record['path'],214311908,record['sector'])
  if curve['sha256']!=record['sha256']:raise ValueError('Original hash differs')
  source=dict(tic=214311908,sector=record['sector'],path=record['path'],sha256=record['sha256'],url=record['url']);sources.append(source)
  if not args.pool:
   path=OUT/f'sector{record["sector"]}_screen.json'
   if path.exists():
    d=json.loads(path.read_text())
    if d['methods']!=methods:raise ValueError('Methods changed')
    continue
   scans=screen_arrays(curve['t'],curve['f'],curve['e'],curve['q'])
   # Additional shorter-duration/short-period branch, explicitly separate.
   extra=[]
   for width in [.75,1.5]:
    config=dict(min_period=.2,max_period=3.,durations=(.0125,.02,.03,.05,.075),n_peaks=8,n_refine=32,compute_full_stats=True,grid_drift_durations=.8,width_days=width,scan_single_events=False)
    scan=screen_curve(curve['t'],curve['f'],curve['e'],curve['q'],config=config);scan.update(width_days=width,label='short_period_short_duration',config=config);extra.append(scan)
   scans['scans']+=extra;scans.update(tic=214311908,source=source,metadata=curve['metadata'],fits_checksums=curve['fits_checksums'],methods=methods);f.save(path,scans)
   for scan in scans['scans']:
    print('SCREEN',record['sector'],scan['label'],scan['width_days'],[(round(p['period'],7),round(p.get('local_snr',0),2),round(p.get('local_depth',0)*1e6),p.get('n_positive_events'),p.get('guard_failures'))for p in scan.get('peaks',[])[:4]],flush=True)
 if args.pool:
  for width in [.75,1.5]:
   path=OUT/f'pooled_w{width:g}.json'
   if path.exists():continue
   chunks=[]
   for source in sources:
    a=read_curve(source['path'],214311908,source['sector']);chunk=preprocess_curve(a['t'],a['f'],a['e'],a['q'],config=dict(width_days=width));chunk['sector']=np.full(len(chunk['t']),source['sector']);chunks.append(chunk)
   arrays={k:np.concatenate([c[k]for c in chunks])for k in ['t','f','e','sector']};order=np.argsort(arrays['t']);arrays={k:v[order]for k,v in arrays.items()};np.savez_compressed(OUT/f'pooled_w{width:g}_arrays.npz',**arrays)
   result=search_prepared(arrays['t'],arrays['f'],arrays['e']);result.update(tic=214311908,sources=sources,width_days=width,methods=methods)
   for peak in result['peaks']:
    peak['per_sector']={str(s):event_stats(arrays['t'][arrays['sector']==s],arrays['f'][arrays['sector']==s],arrays['e'][arrays['sector']==s],peak['period'],peak['t0_btjd'],peak['duration_days'])for s in np.unique(arrays['sector'])}
   f.save(path,result);print('POOLED',width,[(round(p['period'],8),round(p['local_snr'],2),round(p['local_depth']*1e6),{s:round(z['local_snr'],2)for s,z in p['per_sector'].items()})for p in result['peaks']],flush=True)

if __name__=='__main__':main()
