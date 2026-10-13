"""Create searchable target/signal ledgers, retaining unvetted ranking statistics."""
import csv,json
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1]

def main():
    out=ROOT/'results';observations=json.loads((ROOT/'data'/'observations.json').read_text())
    target_meta={}
    fs=ROOT/'data'/'field_sample.json'
    if fs.exists():
        for r in json.loads(fs.read_text())['selected']:target_meta[r['tic']]={'name':str(r['tic']),'track':'field','tmag':r['Tmag'],'radius':r['rad'],'teff':r['Teff']}
    with (ROOT/'data'/'catalogs'/'pilot_targets.csv').open() as f:
        for r in csv.DictReader(f):target_meta[int(r['tic_id'])]={'name':r['hostname'],'track':'known_host','tmag':r['sy_tmag'],'radius':r['st_rad'],'teff':r['st_teff']}
    target_meta[307210830]={'name':'L 98-59','track':'control','tmag':10.232,'radius':.303,'teff':3415}
    sig=[];scans=defaultdict(int)
    for family in ('scans','joint','robust'):
        for path in sorted((out/family).glob('*.json')):
            d=json.loads(path.read_text());tic=d['tic'];scans[tic]+=1
            for rank,peak in enumerate(d['peaks'],1):
                sig.append({'tic':tic,'target':target_meta.get(tic,{}).get('name',''),'track':target_meta.get(tic,{}).get('track','unknown'),'method':family,
                    'sector_set':';'.join(map(str,d.get('sectors',[d.get('sector')]))),'width_days':d['width'],'rank':rank,
                    'period_days':peak['period'],'t0_btjd':peak['t0_btjd'],'duration_hours':peak['duration_days']*24,
                    'depth_ppm':peak['depth']*1e6,'depth_error_ppm_diagonal':peak['depth_err']*1e6,'snr_diagonal':peak['snr'],
                    'sde_robust':peak['sde_robust'],'supported_events':peak['n_events'],'odd_even_sigma_diagonal':peak['odd_even_sigma'],
                    'known_masked':d.get('known_masked',True),'status':'screening statistic; see vetted reports','evidence_file':str(path.relative_to(ROOT))})
    sig.sort(key=lambda r:r['snr_diagonal'],reverse=True)
    with (out/'signals.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(sig[0]));writer.writeheader();writer.writerows(sig)
    targets=[]
    for tic in sorted({r['tic'] for r in observations}):
        rows=[r for r in observations if r['tic']==tic]
        best=max((r['snr_diagonal'] for r in sig if r['tic']==tic and r['known_masked']),default=None)
        targets.append({'tic':tic,**target_meta.get(tic,{}),'cached_sectors':';'.join(str(r['sector']) for r in sorted(rows,key=lambda r:r['sector'])),
            'cached_products':len(rows),'scan_variants_completed':scans[tic],'max_unvetted_snr':best,'candidate_status':'no promoted candidate; review pending/completed diagnostics'})
    with (out/'targets.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(targets[0]));writer.writeheader();writer.writerows(targets)
    print('Targets',len(targets),'screened peaks',len(sig),'cached products',len(observations))
if __name__=='__main__':main()
