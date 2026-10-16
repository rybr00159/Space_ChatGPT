"""Inspect current QLP release status and real products before latest-sector search.

Imports the frozen inventory implementation; writes only fresh sector metadata
and this separately named probe. Raw HTTP responses and FITS remain retained.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import tarfile

import numpy as np
import requests
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import multisector_inventory as inventory
import multisector_survey as survey
from multisector_screen import read_curve

OUT = ROOT/'data/multisector/latest_probe'
CACHE = inventory.CACHE/'latest_probe'


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def status_page():
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE/'official_qlp.html'
    manifest_path = path.with_suffix('.manifest.json')
    url = 'https://archive.stsci.edu/hlsp/qlp'
    if path.exists():
        m = json.loads(manifest_path.read_text())
        if sha(path) != m['sha256']:
            raise ValueError('Retained official status page changed')
    else:
        response = requests.get(url, timeout=(15,120))
        path.write_bytes(response.content)
        m = dict(url=url, status=response.status_code, path=str(path),
                 sha256=sha(path), retrieved_utc=utc(), bytes=path.stat().st_size)
        manifest_path.write_text(json.dumps(m, indent=2)+'\n')
        response.raise_for_status()
    return m


def probe_sector(sector):
    rows, manifest = inventory.selected_inventory(sector,8192)
    with Path(manifest['path']).open() as stream:
        official_count = sum(1 for _ in csv.DictReader(stream))
    samples = []
    for row in rows[:2]:
        tic = int(row['#TIC_ID'])
        record = survey.download(dict(tic=tic, sector=sector, ra=row['RA'], dec=row['DEC'],
                                      url=survey.url_for(tic,sector)))
        curve = read_curve(record['path'], tic, sector)
        with fits.open(record['path'], memmap=False) as hd:
            header = hd[1].header
            names = hd[1].data.names
            finite = np.isfinite(curve['t'])
            good_flux = finite & np.isfinite(curve['f']) & (curve['f']>0) & (curve['q']==0)
            good = good_flux & np.isfinite(curve['e']) & (curve['e']>0)
            orbit_ids = sorted(set(int(x) for x in hd[1].data['ORBITID'])) if 'ORBITID' in names else []
            native_step = float(np.median(np.diff(np.sort(curve['t'][finite])))*86400)
            details = dict(tic=tic, sector=sector, path=record['path'], url=record['url'],
                sha256=record['sha256'], raw_bytes=record['raw_bytes'],
                n_rows=len(curve['t']), n_good=int(np.sum(good)),
                n_good_flux=int(np.sum(good_flux)), n_finite_error=int(np.sum(np.isfinite(curve['e']))),
                producer_error_search_usable=bool(np.sum(good)>=1000),
                timesys=header.get('TIMESYS'), timeref=header.get('TIMEREF'),
                bjd_reference=header.get('BJDREFI',0)+header.get('BJDREFF',header.get('BJDREFR',0)),
                timedel=header.get('TIMEDEL'), median_cadence_seconds=native_step,
                min_btjd=float(np.nanmin(curve['t'])), max_btjd=float(np.nanmax(curve['t'])),
                good_flux_min_btjd=float(np.min(curve['t'][good_flux])) if np.any(good_flux) else None,
                good_flux_max_btjd=float(np.max(curve['t'][good_flux])) if np.any(good_flux) else None,
                columns=names, fits_checksums=curve['fits_checksums'], orbit_ids=orbit_ids,
                metadata=curve['metadata'])
            if details['n_good_flux'] < 1000 or not 190<native_step<210:
                raise ValueError('Unexpected latest-sector timing/quality')
            samples.append(details)
    result = dict(sector=sector, official_inventory_rows=official_count, inventory=manifest,
                  selected_preflux_ids=len(rows), raw_probes=samples,
                  note='Two deterministic pre-flux probes establish interface only, not whole-sector data quality or transit sensitivity.')
    (OUT/f's{sector:04d}_probe.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PROBE_READY',sector,official_count,[(s['n_rows'],s['n_good_flux'],s['n_good'],s['orbit_ids']) for s in samples],flush=True)
    return result


def archive_queries(sectors):
    """Preserve exact TIC inputs; static full inventories have retrieval hashes."""
    paths = set(CACHE.glob('*'))
    sources = [Path(__file__), ROOT/'src/multisector_inventory.py',
               ROOT/'src/multisector_screen.py', OUT/'probe_summary.json']
    for sector in sectors:
        manifest = ROOT/f'data/multisector/inventory_s{sector:04d}_manifest.json'
        m = json.loads(manifest.read_text())
        for request in m['raw_metadata_requests']:
            path = Path(request['path'])
            if sha(path) != request['sha256']:
                raise ValueError('Metadata response changed: '+str(path))
            paths.update([path,path.with_suffix('.manifest.json')])
        for kind in ['metadata.csv','eligible.csv','manifest.json']:
            sources.append(manifest.with_name(f'inventory_s{sector:04d}_{kind}'))
        paths.add(Path(m['inventory']['path']).with_suffix('.manifest.json'))
    paths = sorted(p for p in paths if p.is_file())
    archive = OUT/'query_snapshots.tar.xz'
    records = [dict(path=str(p.relative_to(inventory.CACHE)),bytes=p.stat().st_size,sha256=sha(p)) for p in paths]
    with tarfile.open(archive,'w:xz',preset=6) as stream:
        for p in paths:
            stream.add(p,arcname=str(p.relative_to(inventory.CACHE)),recursive=False)
    manifest = dict(created_utc=utc(), archive_path=str(archive.relative_to(ROOT)),
                    sha256=sha(archive), bytes=archive.stat().st_size, files=records,
                    source_files=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in sources],
                    note='Exact public TIC query snapshots/status and failed inventory probes. Complete successful inventories remain cached with static URLs and hashes in the inventory manifests; raw FITS remain separately cached.')
    (OUT/'query_archive.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('ARCHIVED',len(records),'responses/manifests',archive.stat().st_size,'bytes',flush=True)


def checkpoint(mode):
    # Reuse the previously validated exact-byte archive verifier/restorer.
    import multisector_inventory_checkpoint as implementation
    implementation.MANIFEST = OUT/'query_archive.json'
    {'verify':implementation.verify,'restore-queries':implementation.restore_queries}[mode](inventory.CACHE)


def retrieve_inventories(sectors):
    for sector in sectors:
        m = json.loads((ROOT/f'data/multisector/inventory_s{sector:04d}_manifest.json').read_text())['inventory']
        target = Path(m['path'])
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():
            if sha(target)!=m['sha256']:
                raise ValueError('Refusing changed existing inventory: '+str(target))
        else:
            response=requests.get(m['url'],timeout=(15,150))
            body=response.content
            if response.status_code!=200 or hashlib.sha256(body).hexdigest()!=m['sha256']:
                unexpected=target.with_suffix('.unexpected_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
                unexpected.write_bytes(body)
                raise ValueError('Changed/failed live inventory preserved separately: '+str(unexpected))
            target.write_bytes(body)
        meta=target.with_suffix('.manifest.json')
        if not meta.exists():
            meta.write_text(json.dumps(m,indent=2)+'\n')
        print('VERIFIED_INVENTORY',sector,sha(target),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['probe','metadata','all','archive','verify','restore-queries','download-inventories'])
    parser.add_argument('--sectors',type=int,nargs='+',default=[105,96])
    parser.add_argument('--size',type=int,default=16384)
    args=parser.parse_args()
    if any(99<=s<=104 for s in args.sectors):
        raise ValueError('Do not search currently timing-warned QLP sectors 99..104')
    OUT.mkdir(parents=True,exist_ok=True)
    if args.mode in ['probe','all']:
        with ThreadPoolExecutor(max_workers=3) as executor:
            page_future=executor.submit(status_page)
            futures=[executor.submit(probe_sector,s) for s in args.sectors]
            results=[f.result() for f in futures]
            page=page_future.result()
        existing=OUT/'probe_summary.json'
        created=json.loads(existing.read_text())['created_utc'] if existing.exists() else utc()
        summary=dict(created_utc=created, official_status_page=page, sectors=sorted(results,key=lambda x:x['sector']),
                     methods={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'src/multisector_inventory.py',ROOT/'src/multisector_screen.py']})
        (OUT/'probe_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.mode in ['metadata','all']:
        for sector in args.sectors:
            existing=ROOT/f'data/multisector/inventory_s{sector:04d}_manifest.json'
            if existing.exists():
                m=json.loads(existing.read_text())
                if m['selected_ids_before_flux']!=args.size or m['method_sha256']!=sha(Path(inventory.__file__)):
                    raise ValueError('Refusing to overwrite a differently frozen metadata selection')
                for kind in ['metadata.csv','eligible.csv']:
                    if not existing.with_name(f'inventory_s{sector:04d}_{kind}').exists():
                        raise ValueError('Incomplete retained metadata selection')
                print('CACHED_SELECTION_READY',sector,m['eligible_rows'],flush=True)
                continue
            inventory.run(sector,args.size,512,6)
    if args.mode=='archive':
        archive_queries(args.sectors)
    if args.mode in ['verify','restore-queries']:
        checkpoint(args.mode)
    if args.mode=='download-inventories':
        retrieve_inventories(args.sectors)

if __name__=='__main__':
    main()
