"""Select a deterministic small southern cool-dwarf pilot; never download FITS.

Raw queries live in TESS_DATA_DIR/field_queries. Catalog absence is a selection
filter only: it is not evidence of novelty, and TIC parameters are provisional.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time
from datetime import datetime, timezone
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data')) / 'field_queries'
ENDPOINT = 'https://mast.stsci.edu/api/v0/invoke'


def invoke_cached(name, service, params, timeout=90):
    request = {'service': service, 'params': params, 'format': 'json', 'pagesize': 10000, 'page': 1}
    path = CACHE / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        cached = json.loads(path.read_text())
        if cached.get('request') == request and cached.get('response', {}).get('status') == 'COMPLETE':
            return cached
    payload = urlencode({'request': json.dumps(request)}).encode()
    for poll in range(8):
        with urlopen(Request(ENDPOINT, data=payload), timeout=timeout) as response:
            result = json.load(response)
        record = {'request': request, 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                  'endpoint': ENDPOINT, 'response': result}
        path.write_text(json.dumps(record, indent=2) + '\n')
        if result.get('status') == 'COMPLETE':
            paging = result.get('paging', {})
            if paging.get('pagesFiltered', paging.get('pages', 1)) > 1:
                raise RuntimeError('Query is paginated; narrow selection rather than truncate.')
            return record
        if result.get('status') != 'EXECUTING':
            raise RuntimeError('MAST status: ' + str(result.get('status')) + ' ' + str(result.get('msg', '')))
        time.sleep(2)
    raise RuntimeError('MAST query still EXECUTING; raw response retained.')


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def catalog_exclusions():
    excluded = {}
    snapshots = []
    for name, columns in [('confirmed_all.csv', ('tic_id',)), ('confirmed_transiting.csv', ('tic_id',)),
                          ('toi.csv', ('tid',)), ('ctoi.csv', ('TIC ID', 'TIC', 'TIC ID ', 'tid', 'tic_id'))]:
        path = ROOT / 'data' / 'catalogs' / name
        if not path.exists():
            continue
        body = path.read_bytes()
        snapshots.append({'file': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(body).hexdigest(),
                          'bytes': len(body)})
        with path.open(newline='') as stream:
            for row in csv.DictReader(stream):
                value = next((row.get(c) for c in columns if row.get(c)), None)
                match = re.fullmatch(r'(?:TIC\s*)?(\d+)(?:\.0)?', str(value or '').strip())
                if match:
                    excluded.setdefault(int(match.group(1)), []).append(name)
    return excluded, snapshots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=12)
    parser.add_argument('--max-observation-queries', type=int, default=20)
    parser.add_argument('--ra', type=float, default=100)
    parser.add_argument('--dec', type=float, default=-60)
    parser.add_argument('--radius', type=float, default=3)
    parser.add_argument('--min-sectors', type=int, default=4)
    parser.add_argument('--max-contamination', type=float, default=0.05)
    args = parser.parse_args()
    if args.limit <= 0 or args.max_observation_queries <= 0:
        parser.error('Limits must be positive.')
    params = {'columns': '*', 'ra': args.ra, 'dec': args.dec, 'radius': args.radius, 'filters': [
        {'paramName': 'Tmag', 'values': [{'min': 9, 'max': 11.5}]},
        {'paramName': 'Teff', 'values': [{'min': 3200, 'max': 4800}]},
        {'paramName': 'rad', 'values': [{'min': 0.2, 'max': 0.8}]}]}
    filename = 'tic_' + hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:12] + '.json'
    cache = invoke_cached(filename, 'Mast.Catalogs.Filtered.Tic.Position.Rows', params)
    excluded, snapshots = catalog_exclusions()
    if not any(s['file'].endswith('confirmed_all.csv') for s in snapshots):
        raise RuntimeError('All confirmed-planet catalog is required to exclude RV hosts.')
    if not any(s['file'].endswith('toi.csv') for s in snapshots):
        raise RuntimeError('TOI snapshot is required.')
    ranked = []
    initial_rejected = []
    for row in cache['response']['data']:
        tic = int(row['ID'])
        if tic in excluded:
            initial_rejected.append({'tic': tic, 'reason': 'catalog match', 'catalogs': sorted(set(excluded[tic]))})
            continue
        contam = number(row.get('contratio'))
        if contam is None or contam > args.max_contamination:
            initial_rejected.append({'tic': tic, 'reason': 'missing or high TIC contamination', 'contratio': contam})
            continue
        if str(row.get('disposition', '')).upper() in ('DUPLICATE', 'ARTIFACT', 'SPLIT'):
            initial_rejected.append({'tic': tic, 'reason': 'TIC disposition', 'disposition': row.get('disposition')})
            continue
        vals = {k: number(row.get(k)) for k in ('Tmag', 'Teff', 'rad', 'ra', 'dec')}
        if any(v is None for v in vals.values()) or not (9 <= vals['Tmag'] <= 11.5 and 3200 <= vals['Teff'] <= 4800 and .2 <= vals['rad'] <= .8):
            initial_rejected.append({'tic': tic, 'reason': 'parameter validation failed'})
            continue
        ranked.append({'tic': tic, **vals, 'contratio': contam, 'TIC_disposition': row.get('disposition'),
                       'gaia': row.get('GAIA'), 'objType': row.get('objType')})
    ranked.sort(key=lambda row: (row['Tmag'], row['contratio'], row['tic']))
    selected = []
    checked = []
    for row in ranked[:args.max_observation_queries]:
        if len(selected) >= args.limit:
            break
        tic = row['tic']
        try:
            q = invoke_cached(f'caom_tic{tic}.json', 'Mast.Caom.Filtered', {'columns': '*', 'filters': [
                {'paramName': 'obs_collection', 'values': ['TESS']},
                {'paramName': 'target_name', 'values': [str(tic)]}]})
            products = []
            for obs in q['response']['data']:
                if re.search(r'[-_]s_lc\.fits$', str(obs.get('dataURL', ''))):
                    products.append({'sector': obs.get('sequence_number'), 'obsid': obs.get('obsid'),
                                     'obs_id': obs.get('obs_id'), 'data_uri': obs['dataURL'],
                                     't_exptime': obs.get('t_exptime')})
            products.sort(key=lambda obs: (obs['sector'] or 0, obs['data_uri']))
            sectors = sorted({int(obs['sector']) for obs in products if obs['sector'] is not None})
            result = {**row, 'sectors': sectors, 'sector_count': len(sectors), 'products': products,
                      'caom_query': str((CACHE / f'caom_tic{tic}.json')), 'selection_status': 'selected' if len(sectors) >= args.min_sectors else 'insufficient SPOC sectors'}
            checked.append(result)
            if len(sectors) >= args.min_sectors:
                selected.append(result)
                print(f"Selected TIC {tic}, Tmag={row['Tmag']:.3f}, Teff={row['Teff']:.0f} K, sectors={sectors}", flush=True)
            else:
                print(f"Checked TIC {tic}: only {len(sectors)} SPOC sectors", flush=True)
        except Exception as exc:
            checked.append({**row, 'selection_status': 'query failed', 'error': f'{type(exc).__name__}: {exc}'})
            print(f'Query failed TIC {tic}: {type(exc).__name__}', flush=True)
        # Preserve incremental target IDs even if the process is interrupted.
        output = {'created_utc': datetime.now(timezone.utc).isoformat(), 'strategy': 'bright cool dwarfs; no match in saved confirmed/TOI catalogs',
                  'selection_rule': vars(args), 'ordering': ['Tmag ascending', 'contratio ascending', 'TIC ascending'],
                  'tic_query_file': str(CACHE / filename), 'tic_query': cache['request'],
                  'catalog_snapshots': snapshots, 'eligible_tic_count': len(ranked), 'selected': selected,
                  'checked': checked, 'initial_rejected': initial_rejected,
                  'limitations': ['Short-cadence selection does not establish little prior analysis.', 'Catalog exclusion does not establish novelty.',
                                  'No variability or pixel-level checks have yet been performed.', 'TIC contamination and stellar parameters are provisional.']}
        (ROOT / 'data' / 'field_sample.json').write_text(json.dumps(output, indent=2) + '\n')
    if not checked:
        raise RuntimeError('No eligible targets; adjust field using a documented selection change.')
    lines = ['# Independent cool-dwarf field selection', '', f"Selection completed at {datetime.now(timezone.utc).isoformat()}.", '',
             f"The field is centered at RA {args.ra}°, Dec {args.dec}° with radius {args.radius}°. TIC cuts: Tmag 9–11.5, Teff 3200–4800 K, radius 0.2–0.8 solar radii, finite contamination ≤ {args.max_contamination}. Exclude all matches in current saved confirmed-planet and TOI catalogs, including dispositions other than planet candidate. Sort remaining stars by magnitude, contamination, then TIC ID. Check up to {args.max_observation_queries} stars, accept the first {args.limit} with at least {args.min_sectors} actual SPOC `-s_lc.fits` observation products. No flux or BLS result influenced selection.", '',
             '| TIC | Tmag | Teff (K) | Radius (solar) | TIC contamination | Sector count |', '| --- | --- | --- | --- | --- | --- |']
    for row in selected:
        lines.append(f"| {row['tic']} | {row['Tmag']:.3f} | {row['Teff']:.0f} | {row['rad']:.3f} | {row['contratio']:.4f} | {row['sector_count']} |")
    lines += ['', f"Selected {len(selected)} targets after checking {len(checked)} of {len(ranked)} eligible catalog targets. Product URIs, sectors, exact request, catalog hashes, exclusion reasons, and observation-query statuses are saved in `data/field_sample.json`. Complete raw JSON responses remain in `{CACHE}`. FITS downloads are a separate pipeline step.", '',
              'This is a comparison pilot; these targets are not proven to be unstudied. Before any novelty claim, check CTOIs, TCEs, aliases, and literature. Incomplete catalog coverage remains a limitation. TIC temperatures/radii/contamination are catalog estimates and require verification for survivors.', '',
              f'Reproduce with `python src/select_field_sample.py --limit {args.limit} --max-observation-queries {args.max_observation_queries} --ra {args.ra} --dec {args.dec} --radius {args.radius} --min-sectors {args.min_sectors} --max-contamination {args.max_contamination}`. Existing exact-query cache entries avoid repeated archive traffic. Apply the same catalog snapshots or record updated hashes. A larger/different field must be logged as a strategy change.']
    (ROOT / 'reports' / 'field_selection.md').write_text('\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
