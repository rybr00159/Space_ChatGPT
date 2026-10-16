"""Pooled searches that explicitly include reconstructed-error observations.

Only repeated targets needing at least one separately recovered observation
enter this branch. Producer-error pooled results remain untouched. No array is
written over a producer-prepared array, and each uncertainty assumption is
part of the result identity.
"""
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,csv,gzip,json,sqlite3

import numpy as np

from fast_screen import preprocess_curve,event_stats
from multisector_screen import read_curve,CONFIG as SCREEN_CONFIG,WIDTHS,seed_reasons
from multisector_combine import search_prepared
from multisector_combine_parallel import atomic_json,atomic_bytes,sha,utc
from multisector_missing_errors import estimate,diagnostics
from multisector_recognition import FreshCatalogIndex

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/workspace/tess-research-data/multisector')
DEST=CACHE/'combined_mixed'
OUT=ROOT/'results/multisector_mixed_pool'
INDEX=None
SOURCES=['src/multisector_mixed_pool.py','src/multisector_combine.py','src/multisector_screen.py',
         'src/fast_screen.py','src/multisector_missing_errors.py']


def inputs():
    grouped=defaultdict(dict)
    for name,assumption in [('survey.sqlite','producer_reported'),('error_recovery.sqlite','empirical_adjacent_difference_MAD')]:
        with sqlite3.connect((CACHE/name).as_uri()+'?mode=ro',uri=True)as conn:
            rows=conn.execute("SELECT result_json FROM observations WHERE status='screened' ORDER BY tic,sector").fetchall()
        for row in rows:
            record=json.loads(row[0]);record['error_assumption']=assumption
            key=record['tic'],record['sector']
            if record['sector']in grouped[record['tic']]:raise ValueError('Duplicate successful source assumption: '+str(key))
            grouped[record['tic']][record['sector']]=record
    return {tic:[rows[s]for s in sorted(rows)]for tic,rows in grouped.items()
            if len(rows)>=2 and any(r['error_assumption']!='producer_reported'for r in rows.values())}


def prepared_sources(tic,rows):
    sources=[]
    for record in rows:
        sector=record['sector'];paths={}
        source=dict(sector=sector,path=record['path'],raw_sha256=record['sha256'],error_assumption=record['error_assumption'])
        if record['error_assumption']=='producer_reported':
            for width in WIDTHS:paths[str(width)]=CACHE/'prepared'/f'tic{tic}_s{sector:04d}_w{width:g}.npz'
        else:
            curve=read_curve(record['path'],tic,sector)
            if curve['sha256']!=record['sha256']:raise ValueError('Original observation differs')
            if not diagnostics(curve)['recovery_eligible']:raise ValueError('Recovery preconditions no longer pass')
            noise=estimate(curve)
            if noise!=record['empirical_error']:raise ValueError('Empirical noise reconstruction differs')
            source.update(empirical_sigma_flux=noise['empirical_sigma_flux'],estimator_sha256=sha(ROOT/'src/multisector_missing_errors.py'))
            for width in WIDTHS:
                path=CACHE/'prepared_mixed'/f'tic{tic}_s{sector:04d}_w{width:g}_empirical.npz'
                if not path.exists():
                    flat=preprocess_curve(curve['t'],curve['f'],np.full(curve['e'].shape,noise['empirical_sigma_flux']),
                                          curve['q'],config=dict(SCREEN_CONFIG,width_days=width))
                    path.parent.mkdir(parents=True,exist_ok=True)
                    np.savez_compressed(path,**{k:flat[k]for k in ['t','f','e','segment']})
                paths[str(width)]=path
        source['prepared']={width:dict(path=str(path),sha256=sha(path))for width,path in paths.items()}
        sources.append(source)
    return sources


def initialize(snapshot):
    global INDEX
    INDEX=FreshCatalogIndex(snapshot=ROOT/snapshot)


def compute(job):
    tic,rows,identity=job;scans=[]
    for width in WIDTHS:
        chunks=[]
        for source in identity['sources']:
            prepared=source['prepared'][str(width)];path=Path(prepared['path'])
            if sha(path)!=prepared['sha256']:raise ValueError('Prepared input differs')
            with np.load(path)as data:chunk={k:data[k]for k in ['t','f','e']}
            chunk['sector']=np.full(len(chunk['t']),source['sector']);chunks.append(chunk)
        joined={k:np.concatenate([c[k]for c in chunks])for k in ['t','f','e','sector']}
        order=np.argsort(joined['t']);joined={k:v[order]for k,v in joined.items()}
        known_keep=INDEX.transit_mask(tic,joined['t'],float(rows[0]['ra']),float(rows[0]['dec']))
        for label,mask in [('unmasked',np.ones(len(joined['t']),bool)),('known_masked_residual',known_keep)]:
            if label!='unmasked'and mask.all():continue
            if mask.sum()<50:continue
            t,f,e=(joined[k][mask]for k in ['t','f','e']);sectors=joined['sector'][mask]
            scan=search_prepared(t,f,e);scan.update(width_days=width,label=label,n_known_masked=int((~mask).sum()))
            for peak in scan['peaks']:
                peak['per_sector']={}
                for sector in sorted(set(sectors)):
                    take=sectors==sector
                    peak['per_sector'][str(sector)]=event_stats(t[take],f[take],e[take],peak['period'],peak['t0_btjd'],peak['duration_days'])
                peak['positive_support_sectors']=sum(s['n_positive_events']>=1 and s['local_snr']>=2.5 for s in peak['per_sector'].values())
            scans.append(scan)
    result=dict(tic=tic,sectors=[r['sector']for r in rows],identity=identity,scans=scans,checked_utc=utc(),
                limitations='At least one uncertainty array reconstructed explicitly; pooled-fitted per-sector statistics are conditional, not independent confirmations or calibrated FAP.')
    path=DEST/f'tic{tic}.json'
    if path.exists():raise FileExistsError('Result appeared during calculation; refusing replacement')
    atomic_json(path,result)
    return tic


def export(completed,methods,snapshot,planned):
    index=FreshCatalogIndex(snapshot=ROOT/snapshot);table=[];sources=[]
    for tic in sorted(completed):
        path=DEST/f'tic{tic}.json';result=json.loads(path.read_text())
        sources.append(dict(tic=tic,path=str(path),sha256=sha(path)))
        for scan in result['scans']:
            for rank,peak in enumerate(scan['peaks'],1):
                match=index.match_signal(tic,peak['period'])
                table.append(dict(tic=tic,sectors=';'.join(map(str,result['sectors'])),width_days=scan['width_days'],label=scan['label'],rank=rank,
                                  period_days=peak['period'],epoch_btjd=peak['t0_btjd'],duration_days=peak['duration_days'],depth_ppm=peak['depth']*1e6,
                                  local_snr=peak['local_snr'],n_events=peak['n_events'],n_positive_events=peak['n_positive_events'],
                                  positive_support_sectors=peak['positive_support_sectors'],seed_retained=not seed_reasons(peak),
                                  known_period_match=match['known_period_match'],known_harmonic_match=match['known_harmonic_match'],
                                  recognition_status=match['recognition_status'],error_assumption='one_or_more_empirical_adjacent_difference_errors',evidence_path=str(path)))
    OUT.mkdir(parents=True,exist_ok=True)
    if table:
        import io
        stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
        atomic_bytes(OUT/'signals.csv.gz',gzip.compress(stream.getvalue().encode(),mtime=0))
    atomic_json(OUT/'summary.json',dict(checked_utc=utc(),methods=methods,catalog_snapshot=snapshot,planned_repeated_tics=planned,
                                      completed_repeated_tics=len(completed),signal_rows=len(table),sources=sources,
                                      both_sector_retained_rows=sum(r['seed_retained']and r['positive_support_sectors']>=2 for r in table),
                                      limitation='Mixed-error branch overlaps the same selected stars; do not add as new unique stars.'))
    print('MIXED POOL',len(completed),'/',planned,flush=True)


def main(workers):
    DEST.mkdir(parents=True,exist_ok=True);grouped=inputs();index=FreshCatalogIndex();snapshot=index.metadata()['fresh_snapshot']
    methods={name:sha(ROOT/name)for name in SOURCES};jobs=[];completed=set()
    for tic,rows in sorted(grouped.items()):
        identity=dict(methods=methods,sources=prepared_sources(tic,rows),catalog_snapshot=snapshot)
        path=DEST/f'tic{tic}.json'
        if path.exists():
            if json.loads(path.read_text())['identity']!=identity:raise ValueError('Accepted mixed identity differs; preserve and fork')
            completed.add(tic)
        else:jobs.append((tic,rows,identity))
    export(completed,methods,snapshot,len(grouped))
    with ProcessPoolExecutor(max_workers=workers,initializer=initialize,initargs=(snapshot,))as pool:
        for future in as_completed([pool.submit(compute,job)for job in jobs]):
            completed.add(future.result());export(completed,methods,snapshot,len(grouped))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);args=parser.parse_args();main(args.workers)
