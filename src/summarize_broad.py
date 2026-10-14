"""Audit and visualize a stable export of the broad photometric survey.

Run after broad_survey.py finishes exporting. Counts refer to unique TICs and
actual searches, never to the archive inventory or the number of BLS peaks.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/broad'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    paths = [OUT / n for n in ('summary.json', 'targets.csv.gz', 'signals.csv.gz', 'observations.jsonl.gz')]
    before = {p.name: digest(p) for p in paths}
    summary = json.loads(paths[0].read_text())
    targets = pd.read_csv(paths[1], dtype={'tic': 'int64'})
    signals = pd.read_csv(paths[2], dtype={'tic': 'int64'})
    selected = pd.read_csv(ROOT / 'data/broad_inventory/qlp_s0001_selection.csv.gz', dtype={'tic': 'int64'})
    if targets.tic.duplicated().any() or selected.tic.duplicated().any():
        raise ValueError('Duplicate TICs: cannot report unique-star counts')
    if not set(targets.tic) <= set(selected.tic):
        raise ValueError('Export contains targets outside the frozen selection')
    statuses = {str(k): int(v) for k, v in targets.status.value_counts().items()}
    if statuses != summary['status_counts']:
        raise ValueError('Inconsistent export; wait for survey export to finish')
    searched_ids = set(targets.loc[targets.status == 'screened', 'tic'])
    if not set(signals.tic) <= searched_ids:
        raise ValueError('Signals assigned to a star without a completed search')
    obs_ids, bytes_total, checksum_statuses = set(), 0, Counter()
    with gzip.open(paths[3], 'rt') as f:
        for line in f:
            row = json.loads(line)
            if row['tic'] in obs_ids:
                raise ValueError('Duplicate observation for single-sector survey')
            obs_ids.add(row['tic'])
            bytes_total += row['raw_bytes']
            checksum_statuses.update(r['status'] for r in row['fits_checksums'])
    retained_ids = set(targets.loc[targets.status.isin(['screened', 'insufficient_data']), 'tic'])
    if obs_ids != retained_ids:
        raise ValueError('Observation manifest and retained-target export disagree')
    if before != {p.name: digest(p) for p in paths}:
        raise ValueError('Survey export changed during audit; rerun after export completes')
    guard_counts = Counter()
    for value in signals.guard_failures.dropna():
        guard_counts.update(str(value).split(';'))
    summary.update(
        audited_utc=datetime.now(timezone.utc).isoformat(), export_sha256=before,
        unique_exported_stars=len(targets), remaining_selected_stars=len(selected) - len(targets),
        retained_raw_products=len(obs_ids), retained_raw_bytes=bytes_total,
        raw_hdu_checksum_status_counts=dict(checksum_statuses),
        stored_screening_peaks=len(signals),
        unmasked_peaks=int((signals.label == 'unmasked_broad').sum()),
        residual_peaks=int((signals.label == 'known_masked_residual').sum()),
        stars_with_stored_residual_peaks=int(signals.loc[signals.label == 'known_masked_residual', 'tic'].nunique()),
        initial_guard_passing_peaks=int(signals.screening_pass.eq(True).sum()),
        guard_failure_counts=dict(guard_counts),
        counts_note='Broad sample only; pilot counts are not added. TIC overlap with the pilot is possible.',
        significance_note='BLS scores and guards do not establish a global false-alarm probability.',
    )
    (OUT / 'audit.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')

    searched = targets.loc[targets.status == 'screened']
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), layout='constrained')
    h = axes[0].hist2d(searched.ra, searched.dec, bins=(72, 36), cmap='viridis')
    axes[0].set(xlabel='Right ascension (degrees)', ylabel='Declination (degrees)', title='Actually searched sector-1 stars')
    fig.colorbar(h[3], ax=axes[0], label='Stars per angular bin (unequal sky areas)')
    valid_mag = searched.tmag[np.isfinite(searched.tmag) & searched.tmag.between(0, 25)]
    axes[1].hist(valid_mag, bins=np.arange(5, 16.25, .25), color='#3c718c')
    axes[1].set(xlabel='QLP header TESS magnitude', ylabel='Searched stars', title='Photometric sample')
    snr = searched.max_unvetted_snr.replace([np.inf, -np.inf], np.nan).dropna()
    axes[2].hist(np.minimum(snr, 30), bins=np.arange(0, 31.01, 1), color='#987844')
    axes[2].set(xlabel='Largest stored diagonal BLS S/N (capped at 30)', ylabel='Searched stars', title='Unvalidated screening scores', yscale='log')
    fig.suptitle(f'{len(searched):,} actual period searches; no significance or discovery inference')
    fig.savefig(OUT / 'survey_overview.png', dpi=150)
    plt.close(fig)
    print(json.dumps({k: summary[k] for k in ('status_counts', 'remaining_selected_stars', 'retained_raw_products', 'retained_raw_bytes', 'stored_screening_peaks', 'stars_with_stored_residual_peaks')}, indent=2))


if __name__ == '__main__':
    main()
