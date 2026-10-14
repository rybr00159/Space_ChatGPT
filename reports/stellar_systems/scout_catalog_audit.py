"""Expanded prior-EB matching for the new stellar scout, retaining snapshots."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
import csv
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/stellar_systems'
RAW = OUT / 'scout_catalogs'
DEEP = [260503860, 115188534, 300652458, 52285789, 115152243, 403150111,
        294206969, 89018908]


def parse(path):
    lines = [line for line in Path(path).read_text().splitlines()
             if line and not line.startswith('#')]
    if len(lines) < 4 or '\t' not in lines[0]:
        return []
    return [{key: value.strip() if value else '' for key, value in row.items()}
            for row in csv.DictReader([lines[0]] + lines[3:], delimiter='\t')]


def fetch(item):
    tic, ra, dec, name, source = item
    url = 'https://vizier.cds.unistra.fr/viz-bin/asu-tsv?' + urlencode({
        '-source': source, '-c': f'{ra} {dec}', '-c.rm': '2',
        '-out.max': 'unlimited', '-out.all': '', '-out.form': 'tsv'})
    path = RAW / f'{name}_tic{tic}.tsv'
    manifest = path.with_suffix('.manifest.json')
    if not path.exists():
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        path.write_bytes(response.content)
        manifest.write_text(json.dumps({'url': url, 'source': source, 'radius_arcmin': 2,
                                        'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                                        'sha256': hashlib.sha256(response.content).hexdigest(),
                                        'bytes': len(response.content), 'http_status': response.status_code}, indent=2) + '\n')
    else:
        if hashlib.sha256(path.read_bytes()).hexdigest() != json.loads(manifest.read_text())['sha256']:
            raise ValueError('Snapshot changed')
    rows = parse(path)
    return {'tic': tic, 'catalog': name, 'rows': rows, 'source': source,
            'raw_file': str(path.relative_to(ROOT)),
            'retrieved_utc': json.loads(manifest.read_text())['retrieved_utc']}


def main():
    targets = pd.read_csv(OUT / 'scout_raw_measurements.csv')
    host_coords = SkyCoord(targets.ra.to_numpy() * u.deg, targets.dec.to_numpy() * u.deg)
    matches = []
    for name in ['shi2022', 'tang2025', 'battley2022']:
        rows = [row for row in parse(RAW / f'{name}.tsv') if row.get('TIC', '').isdigit()]
        if name == 'shi2022':
            coords = SkyCoord([row['RAJ2000'] for row in rows], [row['DEJ2000'] for row in rows],
                              unit=(u.hourangle, u.deg))
            pkey = 'Per'
        else:
            ra_key, dec_key = ('_RA', '_DE') if name == 'battley2022' else ('RAJ2000', 'DEJ2000')
            coords = SkyCoord([float(row[ra_key]) for row in rows] * u.deg,
                              [float(row[dec_key]) for row in rows] * u.deg)
            pkey = 'Period' if name == 'battley2022' else 'Per'
        for _, target in targets.iterrows():
            separations = SkyCoord(target.ra * u.deg, target.dec * u.deg).separation(coords).arcsec
            keep = (separations <= 120) | (np.array([int(row['TIC']) for row in rows]) == int(target.tic))
            for index in np.flatnonzero(keep):
                row = rows[index]
                period = float(row[pkey]) if row[pkey] else None
                matches.append({'tic': int(target.tic), 'catalog': name,
                                'catalog_tic': int(row['TIC']), 'separation_arcsec': float(separations[index]),
                                'exact_identifier_match': int(row['TIC']) == int(target.tic),
                                'period_days': period, 'row': row})
    (OUT / 'scout_expanded_catalog_matches.json').write_text(json.dumps(matches, indent=2) + '\n')
    pd.DataFrame([{key: row[key] for key in row if key != 'row'} for row in matches]).to_csv(
        OUT / 'scout_expanded_catalog_matches.csv', index=False)
    jobs = []
    for tic in DEEP:
        target = targets[targets.tic == tic].iloc[0]
        for name, source in [('vsx', 'B/vsx/vsx'), ('asassn', 'II/366/catv2021'), ('sb9', 'B/sb9/main')]:
            jobs.append((tic, target.ra, target.dec, name, source))
    results = list(ThreadPoolExecutor(max_workers=6).map(fetch, jobs))
    (OUT / 'scout_deep_novelty.json').write_text(json.dumps(results, indent=2) + '\n')
    print('Expanded whole-catalog host/neighbor matches:', len(matches))
    for result in results:
        print(result['tic'], result['catalog'], len(result['rows']), repr(result['rows'])[:1800])


if __name__ == '__main__':
    main()
