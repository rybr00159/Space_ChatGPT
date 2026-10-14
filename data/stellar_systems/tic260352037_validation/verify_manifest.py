"""Record or verify this research checkpoint, including immutable raw inputs.

Default is read-only verification. --record rebuilds a current manifest after
running the scientific modes; it never edits observations or scientific results.
"""
from pathlib import Path
import argparse,sys,json,importlib.metadata as im
from datetime import datetime,timezone
import numpy as np
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
import stellar_giant_validate as x

def record():
    out=x.OUT;data=x.DATA;curves=x.load();sha=x.digest(Path(x.__file__));inputs=[]
    for c in curves:
        r=c['source'];p=Path(r['path'])
        with fits.open(p) as h:
            timesys=h[1].header.get('TIMESYS');bjd=h[1].header.get('BJDREFI',0)+h[1].header.get('BJDREFF',0)
        if timesys!='TDB' or bjd!=2457000:raise ValueError('Wrong time coordinate')
        if abs(np.median(np.diff(c['t']))/c['dt']-1)>.01:raise ValueError('Exposure/timestamp mismatch')
        inputs.append(dict(tic=x.TIC,sector=c['sector'],era=c['era'],path=str(p),sha256=x.digest(p),bytes=p.stat().st_size,url=r.get('url'),native_exposure_seconds=c['dt']*86400,time_system=timesys,bjd_reference=bjd,n_good=len(c['t'])))
    for folder in out.glob('sector*'):
        if not folder.is_dir():continue
        m=json.loads((folder/'raw_manifest.json').read_text());p=Path(m['path'])
        if x.digest(p)!=m['sha256']:raise ValueError('Cube bytes changed')
        inputs.append(dict(kind='pixel_cube',sector=m['sector'],path=str(p),sha256=m['sha256'],bytes=m['bytes'],url=m['url'],checks=m['checks']))
    controls=json.loads((out/'photometry.json').read_text());pool={z['group']:z for z in controls['depth_summaries'] if z['baseline_degree']==2}
    if len(curves)!=41 or len(pool['later']['primary']['sectors'])!=16:raise ValueError('Unexpected target sample')
    if not all(r['depth']>0 for family in ['primary','secondary'] for r in pool['later'][family]['sectors']):raise ValueError('Holdout result changed')
    for name in ['photometry','morphology','pixels','timing_stability','wide_pair_dynamics']:
        if json.loads((out/f'{name}.json').read_text())['code_sha256']!=sha:raise ValueError('Analysis method identity mismatch')
    outputs=[dict(path=str(p.relative_to(ROOT)),sha256=x.digest(p),bytes=p.stat().st_size) for p in out.rglob('*') if p.is_file()]
    mast=json.loads((data/'mast_observations.json').read_text())
    manifest=dict(checked_utc=datetime.now(timezone.utc).isoformat(),python=sys.version,dependencies={k:im.version(k) for k in ['numpy','scipy','astropy','pandas','matplotlib','requests']},code=[dict(path=str(p),sha256=x.digest(ROOT/p)) for p in [Path('src/stellar_giant_validate.py'),Path('src/broad_vet.py'),Path('src/archive.py')]],inputs=inputs,outputs=outputs,coverage_query_url='https://mast.stsci.edu/tesscut/api/v0.1/sector?ra=94.787911339&dec=-57.078220725',coverage_query_sha256=x.digest(data/'tesscut_coverage.json'),mast_target_query=dict(service='Mast.Caom.Filtered',params={'columns':'*','filters':[{'paramName':'obs_collection','values':['TESS']},{'paramName':'target_name','values':['260352037']}]},response_sha256=x.digest(data/'mast_observations.json'),n_rows=len(mast['data'])),checks=['All 41 original QLP byte hashes match immutable manifests and FITS checksums verified where present.','All 41 curves declare BTJD TDB and native exposure matches retained timestamp interval within 1%.','All 4 raw pixel cubes match manifests and science-coverage checks.','All 16 later sectors have positive primary and secondary depths at the old fixed ephemeris.','All five quantitative analysis files identify the current method source hash.','Scientific analyses completed and authored Python compiles.'],limitations=['Raw-cache storage persistence unverified; immutable product URLs and hashes retained.','Clean-cache full rerun not independently executed.','Conditional statistical uncertainties exclude uncalibrated profile/color/intra-pixel systematics.'])
    x.save('analysis_manifest.json',manifest,data)
    print('RECORDED',len(inputs),'immutable raw products;',len(outputs),'derived artifacts')

def verify():
    manifest=json.loads((x.DATA/'analysis_manifest.json').read_text());n=0
    for group in ['code','inputs','outputs']:
        for r in manifest[group]:
            p=Path(r['path']);p=p if p.is_absolute() else ROOT/p
            if x.digest(p)!=r['sha256']:raise ValueError(f'Changed {group} bytes: {p}')
            n+=1
    print('VERIFIED',n,'code/raw/derived byte hashes')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--record',action='store_true');a=p.parse_args()
    if a.record:record()
    else:verify()
