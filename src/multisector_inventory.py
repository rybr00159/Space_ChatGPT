"""Pre-flux, reproducible multi-sector QLP inventory and public TIC metadata.

No variability, transit score or downloaded light curve influences selection.
Raw archive responses and exact requests remain in the scientific cache.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import heapq
import json
import math
import time

import requests

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path('/workspace/tess-research-data/multisector')
ENDPOINT = 'https://mast.stsci.edu/api/v0/invoke'
SERVICE = 'Mast.Catalogs.Filtered.Tic'
COLUMNS = 'ID,ra,dec,Tmag,Teff,rad,contratio,lumclass'


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selection_hash(tic):
    return hashlib.sha256(('multisector-v1:' + str(int(tic))).encode()).hexdigest()


def selected_inventory(sector, size):
    directory = CACHE / 'inventories'
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f'qlp_s{sector:04d}_targets.csv'
    manifest_path = path.with_suffix('.manifest.json')
    url = f'https://archive.stsci.edu/hlsps/qlp/target_lists/s{sector:04d}.csv'
    if not path.exists():
        response = requests.get(url, timeout=(15, 150))
        response.raise_for_status()
        path.write_bytes(response.content)
        manifest_path.write_text(json.dumps(dict(sector=sector, url=url,
            status=response.status_code, path=str(path), sha256=sha(path),
            bytes=path.stat().st_size, retrieved_utc=utc()), indent=2) + '\n')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('status') != 200 or sha(path) != manifest['sha256']:
        raise ValueError('Inventory is failed or changed: ' + str(path))
    with path.open() as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['#TIC_ID', 'RA', 'DEC']:
            raise ValueError('Unexpected official inventory schema')
        rows = heapq.nsmallest(size, reader, key=lambda r: selection_hash(r['#TIC_ID']))
    if len({r['#TIC_ID'] for r in rows}) != len(rows):
        raise ValueError('Duplicate selected TIC identifiers')
    return rows, manifest


def fetch_batch(ids):
    request = dict(service=SERVICE, params=dict(columns=COLUMNS,
        filters=[dict(paramName='ID', values=[str(x) for x in ids])]),
        format='json', pagesize=10000, page=1)
    digest = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
    directory = CACHE / 'tic_metadata'; directory.mkdir(parents=True, exist_ok=True)
    path = directory / (digest + '.json')
    manifest_path = path.with_suffix('.manifest.json')
    if path.exists():
        manifest = json.loads(manifest_path.read_text())
        if sha(path) != manifest['sha256'] or manifest.get('request') != request:
            raise ValueError('Cached TIC response mismatch')
        body = json.loads(path.read_text())
    else:
        for attempt in range(3):
            try:
                response = requests.post(ENDPOINT, data={'request': json.dumps(request)}, timeout=(15, 100))
                if response.status_code != 200:
                    failed = directory / f'{digest}.failed{attempt}.txt'
                    failed.write_bytes(response.content)
                    failed.with_suffix('.manifest.json').write_text(json.dumps(dict(url=ENDPOINT,
                        request=request, status=response.status_code, sha256=sha(failed),
                        retrieved_utc=utc()), indent=2) + '\n')
                response.raise_for_status()
                body = response.json()
                if body.get('status') != 'COMPLETE':
                    raise ValueError('Incomplete TIC query: ' + str(body.get('status')))
                path.write_bytes(response.content)
                manifest = dict(endpoint=ENDPOINT, request=request, status=response.status_code,
                    path=str(path), sha256=sha(path), bytes=len(response.content),
                    retrieved_utc=utc(), rows=len(body.get('data', [])))
                manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
                break
            except (requests.RequestException, ValueError):
                if attempt == 2:
                    raise
                time.sleep(attempt + 1)
    rows = body.get('data', [])
    returned = [int(r['ID']) for r in rows]
    if len(set(returned)) != len(returned) or not set(returned) <= set(ids):
        raise ValueError('Wrong/duplicated IDs returned by exact TIC request')
    if len(rows) > len(ids):
        raise ValueError('Unexpected pagination/truncation')
    return rows, manifest


def numeric(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def run(sector, size=8192, batch_size=512, workers=6):
    selected, inventory = selected_inventory(sector, size)
    ids = [int(r['#TIC_ID']) for r in selected]
    batches = [ids[i:i + batch_size] for i in range(0, len(ids), batch_size)]
    data, manifests = {}, []
    with ThreadPoolExecutor(workers) as executor:
        for future in as_completed([executor.submit(fetch_batch, batch) for batch in batches]):
            rows, manifest = future.result(); manifests.append(manifest)
            for row in rows:
                data[int(row['ID'])] = (row, manifest)
            print('TIC_BATCH', sector, len(data), '/', size, flush=True)
    output = ROOT / 'data/multisector'; output.mkdir(parents=True, exist_ok=True)
    normalized = []
    for rank, row in enumerate(selected):
        tic = int(row['#TIC_ID']); raw, manifest = data.get(tic, ({}, {}))
        tmag, teff, radius = (numeric(raw.get(k)) for k in ['Tmag', 'Teff', 'rad'])
        reasons = []
        if not raw: reasons.append('TIC_query_no_row')
        if tmag is None or tmag > 12.5: reasons.append('missing_or_faint_Tmag')
        if teff is None or not 3000 <= teff <= 6500: reasons.append('missing_or_outside_Teff')
        if radius is None or not .15 <= radius <= 2.5: reasons.append('missing_or_outside_radius')
        normalized.append(dict(sector=sector, tic=tic, ra=float(row['RA']), dec=float(row['DEC']),
            tmag=tmag, teff=teff, radius=radius, contamination=numeric(raw.get('contratio')),
            lumclass=raw.get('lumclass'), eligible=not reasons, rejection_reasons=';'.join(reasons),
            sourcehash=manifest.get('sha256'), inventory_sourcehash=inventory['sha256'],
            selection_hash=selection_hash(tic), preflux_order=rank))
    columns = list(normalized[0])
    for name, rows in [('metadata', normalized), ('eligible', [r for r in normalized if r['eligible']])]:
        path = output / f'inventory_s{sector:04d}_{name}.csv'
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, columns, lineterminator='\n'); writer.writeheader(); writer.writerows(rows)
    result = dict(created_utc=utc(), sector=sector, inventory=inventory,
        selection_rule="First N ascending SHA256('multisector-v1:'+decimal_TIC_ID); same salt in all sectors",
        selected_ids_before_flux=len(ids), returned_metadata_rows=len(data),
        eligible_rows=sum(r['eligible'] for r in normalized), missing_metadata_ids=sorted(set(ids)-set(data)),
        eligibility='Finite Tmag<=12.5, 3000<=Teff<=6500K, .15<=radius<=2.5Rsun; catalog estimates only',
        raw_metadata_requests=sorted(manifests, key=lambda m:m['path']), method_sha256=sha(Path(__file__)),
        warning='Inventory/TIC metadata coverage is not an analyzed light-curve count. No flux affected selection.')
    (output / f'inventory_s{sector:04d}_manifest.json').write_text(json.dumps(result, indent=2) + '\n')
    print('SELECTION_READY', sector, 'eligible', result['eligible_rows'], flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sectors', type=int, nargs='+', default=[83, 86, 95, 98])
    parser.add_argument('--size', type=int, default=8192)
    parser.add_argument('--batch-size', type=int, default=512)
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    for sector in args.sectors:
        run(sector, args.size, args.batch_size, args.workers)
