"""Separate sensitivity analysis for QLP products with unusable error columns.

Read the main multi-sector database strictly read-only. Preserve producer flux,
errors and QUALITY, and explicitly substitute an empirical constant uncertainty
only in a separately labeled analysis. This does not repair the original data
or establish white noise, calibrated significance, or candidate reliability.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import time

import numpy as np

from multisector_screen import read_curve, screen_arrays, seed_reasons, CONFIG, WIDTHS

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path('/workspace/tess-research-data/multisector')
MAIN_DB = CACHE / 'survey.sqlite'
RECOVERY_DB = CACHE / 'error_recovery.sqlite'
OUT = ROOT / 'results/multisector_error_recovery'
ASSUMPTIONS = dict(min_good_flux=400, max_positive_producer_errors_exclusive=20,
    max_adjacent_gap_cadences=2.5, min_difference_pairs=100,
    estimator='1.4826 * median(abs(diff(flux) - median(diff(flux)))) / sqrt(2)',
    weighting='One positive constant in original corrected-flux units per target-sector observation.',
    limitations='Adjacent differences assume locally independent comparable scatter; correlated noise can bias low and rapid variability high. Frozen detrending still applies observed-scatter floors. Empirical errors are not producer errors or a survey false-alarm calibration.')


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signature():
    return {name:sha(ROOT / name) for name in
            ['src/fast_screen.py', 'src/multisector_screen.py', 'src/multisector_missing_errors.py']}


def diagnostics(curve):
    t, f, e, q = (curve[k] for k in ['t', 'f', 'e', 'q'])
    good = np.isfinite(t) & np.isfinite(f) & (f > 0) & (q == 0)
    valid_error = np.isfinite(e) & (e > 0)
    return dict(n_input=len(t), n_quality_zero=int(np.sum(q == 0)),
        n_q0_finite_positive_flux=int(good.sum()),
        n_q0_good_flux_positive_producer_errors=int(np.sum(good & valid_error)),
        n_q0_good_flux_nonfinite_producer_errors=int(np.sum(good & ~np.isfinite(e))),
        n_q0_good_flux_zero_producer_errors=int(np.sum(good & (e == 0))),
        n_q0_good_flux_negative_producer_errors=int(np.sum(good & (e < 0))),
        recovery_eligible=bool(good.sum() >= ASSUMPTIONS['min_good_flux'] and
            np.sum(good & valid_error) < ASSUMPTIONS['max_positive_producer_errors_exclusive']))


def estimate(curve):
    t, f, q = (curve[k] for k in ['t', 'f', 'q'])
    good = np.isfinite(t) & np.isfinite(f) & (f > 0) & (q == 0)
    order = np.argsort(t[good]); ts, fs = t[good][order], f[good][order]
    cadence = curve['metadata'].get('TIMEDEL')
    if not cadence or not np.isfinite(float(cadence)) or float(cadence) <= 0:
        finite_t = np.sort(t[np.isfinite(t)]); dt = np.diff(finite_t)
        cadence = float(np.median(dt[dt > 0]))
    cadence = float(cadence)
    dt = np.diff(ts)
    adjacent = (dt > 0) & (dt <= ASSUMPTIONS['max_adjacent_gap_cadences'] * cadence)
    differences = np.diff(fs)[adjacent]
    if len(differences) < ASSUMPTIONS['min_difference_pairs']:
        raise ValueError('Too few short-gap adjacent differences for empirical error estimate')
    center = float(np.median(differences))
    sigma_diff = float(1.4826 * np.median(np.abs(differences-center)))
    sigma = sigma_diff / np.sqrt(2.)
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('Empirical error estimate is nonpositive or undefined')
    return dict(cadence_days=cadence, n_difference_pairs=len(differences),
        max_adjacent_gap_days=ASSUMPTIONS['max_adjacent_gap_cadences'] * cadence,
        difference_median_flux=center, difference_mad_sigma_flux=sigma_diff,
        empirical_sigma_flux=float(sigma), median_good_flux=float(np.median(fs)),
        empirical_sigma_ppm=float(sigma / np.median(fs) * 1e6),
        error_assumption='empirical_constant_from_adjacent_flux_differences',
        producer_error_column_retained_unchanged=True, quality_column_retained_unchanged=True)


def connect():
    conn = sqlite3.connect(RECOVERY_DB)
    conn.execute('PRAGMA journal_mode=WAL')
    conn.execute('CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT)')
    conn.execute('CREATE TABLE IF NOT EXISTS observations(tic INTEGER,sector INTEGER,status TEXT,result_json TEXT,error TEXT,updated_utc TEXT,PRIMARY KEY(tic,sector))')
    methods = json.dumps(dict(source_sha256=signature(), assumptions=ASSUMPTIONS,
        frozen_screen_config=CONFIG, frozen_screen_widths=WIDTHS), sort_keys=True)
    old = conn.execute("SELECT value FROM metadata WHERE key='methods'").fetchone()
    if old and old[0] != methods:
        raise ValueError('Recovery method changed; preserve and use a separately named database')
    conn.execute('INSERT OR IGNORE INTO metadata VALUES (?,?)', ('methods', methods));conn.commit()
    return conn


def save(conn, tic, sector, status, record, error=None):
    conn.execute('INSERT OR REPLACE INTO observations VALUES(?,?,?,?,?,?)',
        (tic, sector, status, json.dumps(record, allow_nan=False), str(error) if error else None, utc()))
    conn.commit()


def export(conn):
    OUT.mkdir(parents=True, exist_ok=True)
    targets, signals = [], []
    with gzip.open(OUT / 'all_scans.jsonl.gz.tmp', 'wt') as archive:
        for tic, sector, status, raw, error in conn.execute('SELECT tic,sector,status,result_json,error FROM observations ORDER BY sector,tic'):
            record = json.loads(raw); diag = record.get('producer_error_diagnostics', {})
            noise = record.get('empirical_error', {})
            row = dict(tic=tic, sector=sector, recovery_status=status, error=error,
                raw_sha256=record.get('sha256'), raw_path=record.get('path'),
                original_record_sha256=record.get('original_record_sha256'), **diag, **noise)
            scans = record.get('scans', [])
            if scans:
                row.update(n_clean=max(s['n_clean'] for s in scans),
                           n_searched=max(s['n_searched'] for s in scans))
            targets.append(row)
            archive.write(json.dumps(record, allow_nan=False) + '\n')
            for scan in scans:
                for rank, peak in enumerate(scan.get('peaks', []), 1):
                    reasons = seed_reasons(peak)
                    signals.append(dict(tic=tic, sector=sector, width_days=scan['width_days'],
                        rank=rank, period_days=peak['period'], epoch_btjd=peak['t0_btjd'],
                        duration_days=peak['duration_days'], depth_ppm=peak['depth']*1e6,
                        local_depth_ppm=(peak.get('local_depth') or 0)*1e6,
                        local_depth_error_ppm=(peak.get('local_depth_err') or 0)*1e6,
                        snr_diagonal=peak['snr'], local_snr=peak.get('local_snr'),
                        n_events=peak['n_events'], n_positive_events=peak['n_positive_events'],
                        odd_even_sigma=peak.get('odd_even_sigma'),
                        event_snr_max_fraction=peak.get('event_snr_max_fraction'),
                        seed_retained=not reasons, rejection_reasons=';'.join(reasons),
                        raw_path=record['path'], raw_sha256=record['sha256'],
                        error_assumption='empirical_adjacent_difference_MAD',
                        empirical_sigma_flux=noise['empirical_sigma_flux']))
    (OUT / 'all_scans.jsonl.gz.tmp').replace(OUT / 'all_scans.jsonl.gz')
    for name, rows in [('targets.csv', targets), ('signals.csv', signals),
                       ('candidate_seeds.csv', [r for r in signals if r['seed_retained']])]:
        columns = list(dict.fromkeys(k for row in rows for k in row)) or ['tic','sector']
        with (OUT / name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, columns, lineterminator='\n');writer.writeheader();writer.writerows(rows)
    summary = dict(updated_utc=utc(), main_db_read_only=str(MAIN_DB),
        recovery_db=str(RECOVERY_DB), methods=signature(), assumptions=ASSUMPTIONS,
        counts=dict(Counter(r['recovery_status'] for r in targets)),
        sector_counts={str(s):dict(Counter(r['recovery_status'] for r in targets if r['sector']==s))
                       for s in sorted({r['sector'] for r in targets})},
        additionally_period_searched_observations=sum(r['recovery_status']=='screened' for r in targets),
        unique_additionally_period_searched_tics=len({r['tic'] for r in targets if r['recovery_status']=='screened'}),
        seed_rows=sum(r['seed_retained'] for r in signals),
        seed_tics=len({r['tic'] for r in signals if r['seed_retained']}),
        interpretation='Alternative uncertainty-assumption search, separate from frozen producer-error main analysis. No candidate validation or discovery claim.')
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('RECOVERY_CHECKPOINT', summary['counts'], 'seed_tics', summary['seed_tics'], flush=True)


def run(limit=0, sectors=None, report_every=20):
    # Copy a consistent main-DB snapshot and close it immediately; never write or
    # pin a long-lived read transaction while the primary survey is running.
    main = sqlite3.connect(MAIN_DB.resolve().as_uri() + '?mode=ro', uri=True)
    main.execute('PRAGMA query_only=ON'); main.execute('BEGIN')
    source_methods = dict(main.execute('SELECT key,value FROM metadata'))
    rows = main.execute("SELECT tic,sector,result_json FROM observations WHERE status='insufficient_data' AND result_json IS NOT NULL ORDER BY sector,tic").fetchall()
    main.rollback();main.close()
    conn = connect(); done = set(conn.execute('SELECT tic,sector FROM observations'))
    todo = [(tic, sector, raw) for tic, sector, raw in rows if (tic,sector) not in done and (not sectors or sector in sectors)]
    if limit: todo = todo[:limit]
    print('RECOVERY_PLAN', len(todo), 'new insufficient_data records', flush=True)
    for count, (tic, sector, raw) in enumerate(todo, 1):
        original = json.loads(raw)
        base = {k:original[k] for k in ['tic','sector','ra','dec','path','url','sha256','raw_bytes']}
        base.update(original_status='insufficient_data', original_record_sha256=hashlib.sha256(raw.encode()).hexdigest(),
                    original_main_methods=source_methods, source_sha256=signature(), checked_utc=utc())
        try:
            curve = read_curve(original['path'], tic, sector)
            if curve['sha256'] != original['sha256']:
                raise ValueError('Original raw product changed')
            base['producer_error_diagnostics'] = diagnostics(curve)
            base.update(metadata=curve['metadata'], fits_checksums=curve['fits_checksums'], flux_column=curve['column'])
            if not base['producer_error_diagnostics']['recovery_eligible']:
                save(conn, tic, sector, 'ineligible', base); continue
            noise = estimate(curve); base['empirical_error'] = noise
            empirical_errors = np.full(curve['e'].shape, noise['empirical_sigma_flux'], dtype=float)
            result = screen_arrays(curve['t'], curve['f'], empirical_errors, curve['q'])
            result.update(base)
            status = 'screened' if any(s['status']=='screened' for s in result['scans']) else 'insufficient_data_after_recovery'
            result['status'] = status
            save(conn, tic, sector, status, result)
        except Exception as exc:
            save(conn, tic, sector, 'error', base, exc)
            print('RECOVERY_ERROR', tic, sector, str(exc), flush=True)
        if count % report_every == 0:
            export(conn)
    export(conn);conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--sectors', nargs='+', type=int)
    parser.add_argument('--report-every', type=int, default=20)
    args = parser.parse_args()
    run(args.limit, args.sectors, args.report_every)
