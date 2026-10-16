#!/usr/bin/env python3
"""Audit historical coverage and rerun frozen screening on a reported candidate.

This is a retrospective recovery control, not a blind discovery or independent
validation of the paper's false-positive probabilities. Original survey outputs
are never modified. Sector 98 uses the current QLP DET_FLUX schema adapter.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import gzip
import hashlib
import json
import os
import sqlite3
import sys

import numpy as np
import requests
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive import verify_fits_checksums
from fast_screen import screen_curve, preprocess_curve, event_stats, phase
from broad_survey import CONFIG
from broad_prioritize import classify
from broad_crossmatch import CatalogIndex

OUT = ROOT / 'candidates/tic4206066'
TIC = 4206066
PERIOD, EPOCH, DURATION = 3.182785, 3998.736, 1.87 / 24
PAPER = Path('/workspace/tess-research-data/tic4206066/manuscript.pdf')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coverage_audit():
    checks = []
    for name, key in [
        ('data/broad_inventory/qlp_s0001_targets.csv.gz', '#TIC_ID'),
        ('data/broad_inventory/qlp_s0001_selection.csv.gz', 'tic'),
        ('results/broad/targets.csv.gz', 'tic'),
        ('results/broad/signals.csv.gz', 'tic'),
        ('results/targets.csv', 'tic'),
    ]:
        path = ROOT / name
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(path, 'rt') as stream:
            rows = csv.DictReader(stream)
            if key not in rows.fieldnames:
                raise ValueError('Unexpected target-ID column in ' + name)
            count, matches = 0, 0
            for row in rows:
                count += 1
                matches += int(row[key]) == TIC
        checks.append(dict(path=name, sha256=sha(path), rows=count, tic_matches=matches))
    db = Path('/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite')
    if db.exists():
        with sqlite3.connect(db.as_uri() + '?mode=ro', uri=True) as conn:
            rows = conn.execute('SELECT tic,status FROM stars WHERE tic=?', (TIC,)).fetchall()
        checks.append(dict(path=str(db), query='SELECT tic,status FROM stars WHERE tic=4206066', matches=rows))
    return checks


def retrieve(record, allow_download):
    path = Path(record['path'])
    if not path.exists() and allow_download:
        response = requests.get(record['url'], timeout=(20, 120))
        response.raise_for_status()
        body = response.content
        if hashlib.sha256(body).hexdigest() != record['sha256']:
            changed = path.with_suffix('.changed')
            changed.parent.mkdir(parents=True, exist_ok=True)
            changed.write_bytes(body)
            raise ValueError('Archive changed; unexpected bytes retained separately')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
    if sha(path) != record['sha256']:
        raise ValueError('Original observation hash mismatch: ' + str(path))
    return path


def plot(curves):
    fig, axes = plt.subplots(3, 2, figsize=(11, 8), squeeze=False)
    for ax, (sector, data) in zip(axes, curves):
        t, f, e = (data[k] for k in ('t', 'f', 'e'))
        x = phase(t, PERIOD, EPOCH) * 24
        use = abs(x) < 5
        ax[0].scatter(t, (f - 1) * 1e6, s=1, alpha=.5)
        ax[0].set(xlabel='BTJD (TDB)', ylabel='Detrended flux (ppm)', title=f'Sector {sector}')
        ax[1].scatter(x[use], (f[use] - 1) * 1e6, s=2, alpha=.15)
        edges = np.arange(-5, 5.001, .25)
        xs, ys, es = [], [], []
        for left, right in zip(edges[:-1], edges[1:]):
            take = (x >= left) & (x < right)
            if not take.any():
                continue
            weight = 1 / e[take] ** 2
            xs.append((left + right) / 2)
            ys.append(np.average(f[take] - 1, weights=weight) * 1e6)
            es.append(1 / np.sqrt(weight.sum()) * 1e6)
        ax[1].errorbar(xs, ys, yerr=es, fmt='o', ms=3, color='black')
        for edge in [-DURATION * 12, DURATION * 12]:
            ax[1].axvline(edge, color='tab:red', lw=.8)
        ax[1].set(xlabel='Hours from published 3.182785-day ephemeris', ylabel='Flux (ppm)', xlim=(-5, 5), ylim=(-2000, 2000))
    fig.suptitle('Retrospective recovery control; published ephemeris used for folding')
    fig.tight_layout()
    fig.savefig(OUT / 'recovery.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    index = CatalogIndex()
    results, curves, observations = [], [], []
    for manifest in sorted(OUT.glob('*_manifest.json')):
        record = json.loads(manifest.read_text())
        if record.get('http_status') != 200 or 'sector' not in record:
            continue
        path = retrieve(record, args.download)
        with fits.open(path, memmap=False) as hd:
            hd.verify('exception')
            integrity = verify_fits_checksums(path, hd)
            header, data = hd[1].header, hd[1].data
            if header.get('TIMESYS') != 'TDB' or header.get('BJDREFI', 0) + header.get('BJDREFF', 0) != 2457000:
                raise ValueError('Unexpected time system/reference')
            if int(header.get('TICID', hd[0].header.get('TICID'))) != TIC:
                raise ValueError('Wrong star')
            column = 'DET_FLUX' if record['sector'] == 98 else 'PDCSAP_FLUX'
            t, f, e = (np.asarray(data[k], dtype=float) for k in ('TIME', column, column + '_ERR'))
            q = np.asarray(data['QUALITY'])
            metadata = {k: header.get(k, hd[0].header.get(k)) for k in ['RADIUS', 'TEFF', 'TESSMAG', 'PROCVER']}
        obs = dict(**record, fits_checksums=integrity, column=column, n_input=len(t))
        observations.append(obs)
        print('OBSERVATION', record['sector'], column, 'verified', flush=True)
        for label, config in [('original_broad', CONFIG), ('refined_defaults', {})]:
            result = screen_curve(t, f, e, quality=q, config=config)
            # Run the actual old prioritizer against its retained catalogue snapshots.
            # It cannot establish novelty relative to the supplied public manuscript.
            parent = dict(tic=TIC, ra=80.8598409617, dec=-1.3272782749, metadata=metadata)
            selections = []
            for peak in result['peaks']:
                row, failures = classify(parent, result, peak, index)
                selections.append(dict(period_days=peak['period'], retained_by_prioritizer=row is not None,
                                       rejection_reasons=failures, category=row.get('category') if row else None))
            result.update(sector=record['sector'], label=label, flux_column=column, prioritizer=selections)
            results.append(result)
            print('SCAN', record['sector'], label, [(round(p['period'], 6), round(p['local_snr'], 2)) for p in result['peaks']], flush=True)
        prepared = preprocess_curve(t, f, e, quality=q)
        fixed = event_stats(prepared['t'], prepared['f'], prepared['e'], PERIOD, EPOCH, DURATION)
        results.append(dict(sector=record['sector'], label='published_3d_ephemeris', period_days=PERIOD,
                            epoch_btjd=EPOCH, duration_days=DURATION, **fixed))
        curves.append((record['sector'], prepared))
        np.savez_compressed(OUT / f'sector{record["sector"]}_prepared.npz', **{k: prepared[k] for k in ['t', 'f', 'e', 'segment']})
    summary = dict(checked_utc=datetime.now(timezone.utc).isoformat(),
                   method_sha256={name: sha(ROOT / name) for name in ['src/fast_screen.py', 'src/broad_survey.py', 'src/broad_prioritize.py', 'scripts/tic4206066_recovery.py']},
                   manuscript=dict(title='Two Transit-like Signals in TESS Photometry of the Nearby K Dwarf TIC 4206066',
                                   author='Pavel Rabtsevich', draft_date='2026-09-25', path=str(PAPER), sha256=sha(PAPER) if PAPER.exists() else None),
                   historical_coverage=coverage_audit(), observations=observations, runs=results,
                   limitations=['Target and data sectors selected from the supplied paper; not a discovery.',
                                'No signal injection was used in these actual-data scans.',
                                'Blind-period scans do not use the published ephemeris; plotted folds and fixed tests do.',
                                'QLP98 DET_FLUX replaces legacy KSPSAP_FLUX at the file-reading interface only.',
                                'Reported diagonal SNR is not LEO-Vetter MES or a calibrated false-alarm probability.',
                                'No independent centroid, transit fit, companion model or TRICERATOPS validation performed.',
                                'Cached-catalogue prioritizer does not include all later public source reports.'])
    (OUT / 'recovery.json').write_text(json.dumps(summary, indent=2) + '\n')
    plot(curves)


if __name__ == '__main__':
    main()
