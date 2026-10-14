"""Checkpointed EB-oriented reinterpretation of retained QLP sector-1 flux.

This is a new selection experiment, not a modification of the frozen planetary
search. Catalog absence and conditional eclipse statistics never establish
novelty or a stellar companion. Raw files and the authoritative survey are read
only. Small-aperture/large-aperture behavior flags blends, not source identities.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

from broad_crossmatch import CatalogIndex
from broad_vet import local_quadratic, mad, phase, variability_diagnostics

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data'))
OUT = ROOT / 'results/stellar_systems'
METHOD = 'retained-QLP-EB-scout-v1'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, obj):
    def clean(value):
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [clean(item) for item in value]
        if isinstance(value, (float, np.floating)):
            return float(value) if np.isfinite(value) else None
        if isinstance(value, np.generic):
            return value.item()
        return value
    temp = Path(path).with_suffix('.json.tmp')
    temp.write_text(json.dumps(clean(obj), indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def select():
    signals = pd.read_csv(ROOT / 'results/broad/signals.csv.gz')
    targets = pd.read_csv(ROOT / 'results/broad/targets.csv.gz')
    # Retain broad eclipses that the planetary morphology/duty guard rejected,
    # but require repeatable positive events and forbid single-event domination.
    use = ((signals.label == 'unmasked_broad') & (signals.depth_ppm >= 3000)
           & (signals.depth_ppm < 500000) & (signals.snr_diagonal >= 15)
           & (signals.local_snr >= 10) & (signals.n_events >= 3)
           & (signals.n_positive_events / signals.n_events >= .7)
           & (signals.duration_hours / 24 / signals.period_days < .3)
           & ~signals.guard_failures.fillna('').str.contains('single_event_dominates'))
    selected = signals[use].sort_values('snr_diagonal', ascending=False).drop_duplicates('tic')
    columns = ['tic', 'ra', 'dec', 'tmag', 'radius', 'is_eb_host',
               'is_confirmed_host', 'is_toi_host', 'is_ctoi_host', 'is_published_host']
    selected = selected.merge(targets[columns], on='tic', validate='one_to_one')
    selected = selected[~selected.is_eb_host & ~selected.is_confirmed_host].copy()
    # Moderate depths are particularly useful for blends/missed second stars;
    # deep eclipses remain as a separate stratum and positive controls.
    selected['depth_stratum'] = np.where(selected.depth_ppm < 50000, 'moderate', 'deep')
    selected['pilot_priority'] = (np.log10(selected.snr_diagonal.clip(lower=1))
                                  + .2 * np.log10(selected.n_events)
                                  - .5 * (selected.is_toi_host | selected.is_ctoi_host
                                          | selected.is_published_host).astype(float))
    return selected.sort_values('pilot_priority', ascending=False).reset_index(drop=True)


def read(record):
    path = Path(record['path'])
    if digest(path) != record['sha256']:
        raise ValueError('Immutable raw observation checksum mismatch')
    with fits.open(path, memmap=False) as hd:
        hd.verify('exception')
        h, data = hd[1].header, hd[1].data
        if h.get('TICID', hd[0].header.get('TICID')) != record['tic'] or h.get('SECTOR', hd[0].header.get('SECTOR')) != 1:
            raise ValueError('Target/sector mismatch')
        if h.get('TIMESYS') != 'TDB' or h.get('BJDREFI') != 2457000:
            raise ValueError('Unexpected time coordinate')
        arrays = {key: np.asarray(data[key], float) for key in
                  ['TIME', 'QUALITY', 'SAP_FLUX', 'KSPSAP_FLUX', 'KSPSAP_FLUX_ERR',
                   'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']}
    return arrays


def normalized(t, flux, err, quality):
    good = (np.isfinite(t) & np.isfinite(flux) & np.isfinite(err)
            & (quality == 0) & (flux > 0) & (err > 0))
    t, flux, err = t[good], flux[good], err[good]
    order = np.argsort(t)
    t, flux, err = t[order], flux[order], err[order]
    norm = np.median(flux)
    return t, flux / norm, err / norm


def refine(t, flux, err, period):
    # A broad running median removes only slow changes; all local diagnostics
    # below instead fit each eclipse's native SAP baseline independently.
    cadence = np.median(np.diff(t))
    width = max(3., 3 * period)
    size = max(5, int(width / cadence) // 2 * 2 + 1)
    trend = median_filter(flux, size=size, mode='nearest')
    ff = flux / trend
    ee = np.maximum(err / trend, max(mad(np.diff(ff)) / np.sqrt(2), np.median(err)))
    periods = np.linspace(period * .99, period * 1.01, 151)
    durations = np.array([.04, .07, .10, .16, .25, .35])
    durations = durations[durations < .28 * periods.min()]
    bls = BoxLeastSquares(t, ff, ee)
    result = bls.power(periods, durations, objective='snr', oversample=10)
    j = int(np.argmax(result.power))
    return {'period_days': float(result.period[j]), 'epoch_btjd': float(result.transit_time[j]),
            'duration_days': float(result.duration[j]), 'flattened_depth': float(result.depth[j]),
            'diagonal_bls_snr': float(result.power[j]),
            'period_grid_boundary': j in (0, len(periods) - 1)}


def measure(t, flux, err, period, epoch, duration):
    try:
        stats, _ = local_quadratic(t, flux, err, period, epoch, duration)
        return stats
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
        return {'n_events': 0, 'depth': None, 'depth_err': None, 'snr': 0,
                'events': [], 'error': repr(exc)}


def analyze(record, selection, index):
    arrays = read(record)
    t, f, e = normalized(arrays['TIME'], arrays['SAP_FLUX'], arrays['KSPSAP_FLUX_ERR'], arrays['QUALITY'])
    tk, fk, ek = normalized(arrays['TIME'], arrays['KSPSAP_FLUX'], arrays['KSPSAP_FLUX_ERR'], arrays['QUALITY'])
    seed = float(selection['period_days'])
    refined = refine(tk, fk, ek, seed)
    p, t0, duration = (refined[key] for key in ['period_days', 'epoch_btjd', 'duration_days'])
    primary = measure(t, f, e, p, t0, duration)
    secondary = measure(t, f, e, p, t0 + p / 2, duration)
    # Two parity groups refer to consecutive events at the photometric period.
    # A formal difference is only conditional on this ephemeris/shape choice.
    parity = [measure(t, f, e, 2 * p, t0 + k * p, duration) for k in (0, 1)]
    contrast = None
    if all(x.get('depth') is not None and x['n_events'] >= 2 for x in parity):
        diff = parity[0]['depth'] - parity[1]['depth']
        error = np.hypot(parity[0]['depth_err'], parity[1]['depth_err'])
        contrast = {'depth_difference': float(diff), 'depth_difference_err': float(error),
                    'conditional_sigma': float(abs(diff) / error)}
    apertures = {}
    for key in ['KSPSAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']:
        ta, fa, ea = normalized(arrays['TIME'], arrays[key], arrays['KSPSAP_FLUX_ERR'], arrays['QUALITY'])
        stats = measure(ta, fa, ea, p, t0, duration)
        apertures[key] = {k: stats.get(k) for k in ['n_events', 'depth', 'depth_err', 'snr']}
    catalog = index.match_signal(record['tic'], p, record['ra'], record['dec'], tolerance=.01)
    # Rank only robust repeated raw events. Primary eclipses alone are not proof
    # of a star: a secondary/parity contrast increases follow-up priority.
    sec_snr = secondary.get('snr', 0.) if secondary.get('n_events', 0) >= 3 else 0.
    parity_snr = (contrast or {}).get('conditional_sigma', 0.)
    robust = (primary.get('n_events', 0) >= 3 and primary.get('snr', 0) >= 10
              and primary.get('depth', 0.) is not None and primary.get('depth', 0.) >= .003
              and primary.get('n_positive_events', 0) / max(1, primary.get('n_events', 0)) >= .7)
    score = (min(primary.get('snr', 0), 100.) / 10 + min(max(sec_snr, parity_snr), 50.)) if robust else 0.
    return {'tic': record['tic'], 'checked_utc': utc(), 'raw_sha256': record['sha256'],
            'raw_path': record['path'], 'retrieval_url': record['url'],
            'original_selection': selection, 'refined': refined, 'sap_primary': primary,
            'sap_half_phase_secondary': secondary, 'sap_parity': parity,
            'parity_contrast': contrast, 'apertures': apertures,
            'raw_variability': variability_diagnostics(t, f, e, p, t0, duration),
            'catalog': catalog, 'robust_repeat_scout_pass': bool(robust), 'followup_score': float(score),
            'status': 'unvetted_stellar_signal' if robust else 'raw_remeasurement_failed',
            'limitations': 'Single sector; conditional errors, not global significance. No pixel localization, current live EB/VSX membership, independent-sector persistence, physical stellar model or exhaustive novelty check.'}


def summary(records):
    rows = []
    for r in records:
        primary, secondary = r['sap_primary'], r['sap_half_phase_secondary']
        parity = r.get('parity_contrast') or {}
        rows.append({'tic': r['tic'], 'ra': r['original_selection']['ra'], 'dec': r['original_selection']['dec'],
                     'tmag': r['original_selection']['tmag'], **r['refined'],
                     'raw_depth_ppm': (primary.get('depth') or 0) * 1e6,
                     'raw_depth_error_ppm': (primary.get('depth_err') or 0) * 1e6,
                     'raw_primary_snr': primary.get('snr'), 'raw_primary_events': primary['n_events'],
                     'half_phase_depth_ppm': (secondary.get('depth') or 0) * 1e6,
                     'half_phase_error_ppm': (secondary.get('depth_err') or 0) * 1e6,
                     'half_phase_snr': secondary.get('snr'), 'half_phase_events': secondary['n_events'],
                     'parity_difference_ppm': (parity.get('depth_difference') or 0) * 1e6,
                     'parity_conditional_sigma': parity.get('conditional_sigma'),
                     'is_eb_host': r['catalog']['is_eb_host'], 'is_toi_host': r['catalog']['is_toi_host'],
                     'is_ctoi_host': r['catalog']['is_ctoi_host'],
                     'nearby_catalog_object': r['catalog']['nearby_catalog_object'],
                     'robust_repeat_scout_pass': r['robust_repeat_scout_pass'],
                     'followup_score': r['followup_score'], 'status': r['status']})
    return pd.DataFrame(rows).sort_values('followup_score', ascending=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int, default=192)
    parser.add_argument('--deep-limit', type=int, default=64)
    parser.add_argument('--db', type=Path, default=CACHE / 'broad/qlp_s0001/survey.sqlite')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    selection = select()
    selection.to_csv(OUT / 'scout_selection.csv', index=False)
    # Fixed stratum assignment precedes any new native-flux result inspection.
    moderate = selection[selection.depth_stratum == 'moderate'].head(max(0, args.limit - args.deep_limit))
    deep = selection[selection.depth_stratum == 'deep'].head(args.deep_limit)
    pilot = pd.concat([moderate, deep]).sort_values('pilot_priority', ascending=False)
    pilot.to_csv(OUT / 'scout_pilot.csv', index=False)
    index = CatalogIndex(position_radius_arcsec=120)
    method_hash = digest(__file__)
    metadata = {'created_utc': utc(), 'method': METHOD, 'source_sha256': method_hash,
                'selection_stars': len(selection), 'planned_pilot_stars': len(pilot),
                'survey_stars': 160000, 'survey_period_searched': 159505,
                'inputs': {name: digest(ROOT / 'results/broad' / name)
                           for name in ['signals.csv.gz', 'targets.csv.gz']},
                'catalogs': index.metadata(), 'paid_calls': 0, 'discovery_claim': False}
    atomic_json(OUT / 'scout_manifest.json', metadata)
    db = sqlite3.connect('file:' + str(args.db) + '?mode=ro', uri=True)
    records = []
    for i, selected in enumerate(pilot.to_dict('records'), 1):
        tic = int(selected['tic'])
        path = OUT / f'scout_tic{tic}.json'
        if path.exists():
            record = json.loads(path.read_text())
            if record.get('scout_source_sha256') != method_hash:
                raise ValueError('Existing scout result uses another source version; preserve it')
        else:
            original = json.loads(db.execute('SELECT result_json FROM stars WHERE tic=?', (tic,)).fetchone()[0])
            record = analyze(original, selected, index)
            record['scout_source_sha256'] = method_hash
            atomic_json(path, record)
        records.append(record)
        if i % 10 == 0 or i == len(pilot):
            table = summary(records)
            table.to_csv(OUT / 'scout_raw_measurements.csv', index=False)
            metadata.update(completed_pilot_stars=i, raw_repeat_scout_pass=int(table.robust_repeat_scout_pass.sum()),
                            last_checkpoint_utc=utc())
            atomic_json(OUT / 'scout_manifest.json', metadata)
            print(f'{i}/{len(pilot)} native curves; {table.robust_repeat_scout_pass.sum()} repeated raw leads', flush=True)
    print(summary(records).head(30).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
