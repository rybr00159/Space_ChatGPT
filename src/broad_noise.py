"""Resumable, one-process real-QLP noise and paired injection assessment.

Does not change the active survey or its frozen screening implementation.
Empirical surrogate lead rates are not planetary FAPs or discovery claims.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sqlite3
import platform
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import astropy
import scipy
from astropy.io import fits
from scipy.ndimage import median_filter

from broad_crossmatch import CatalogIndex
from fast_screen import FastScreenConfig, phase, preprocess_curve, screen_curve

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data')) / 'broad/noise'
DB = CACHE.parent / 'qlp_s0001/survey.sqlite'
CONFIG = {'n_refine': 0, 'n_peaks': 2, 'compute_full_stats': False}
SEED = 202610137
EXCLUDE = ('is_confirmed_host', 'is_toi_host', 'is_ctoi_host', 'is_eb_host',
           'is_diamante_host', 'is_published_host')
METHODS = ('block_permutation_1', 'block_permutation_2', 'block_sign_1', 'block_sign_2')


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def frozen_sample(count, paired_stars):
    path = CACHE / 'sample.json'
    source_hash = sha(ROOT / 'src/fast_screen.py')
    if path.exists():
        sample = json.loads(path.read_text())
        if sample['method_sha256'] != source_hash or sample['count_requested'] != count or sample['injection_stars'] != paired_stars:
            raise ValueError('Preserve frozen experiment; source or selection changed')
        return sample
    db = sqlite3.connect(DB.resolve().as_uri() + '?mode=ro', uri=True)
    db.execute('PRAGMA query_only=ON')
    db.execute('BEGIN')
    counts = dict(db.execute('SELECT status,count(*) FROM stars GROUP BY status'))
    index = CatalogIndex()
    selected = []
    rejected = Counter()
    with gzip.open(ROOT / 'data/broad_inventory/qlp_s0001_selection.csv.gz', 'rt') as stream:
        for entry in csv.DictReader(stream):
            tic = int(entry['tic'])
            row = db.execute("SELECT result_json FROM stars WHERE tic=? AND status='screened'", (tic,)).fetchone()
            if not row:
                rejected['not_successfully_screened_at_snapshot'] += 1
                continue
            record = json.loads(row[0])
            flags = index.match_target(tic, float(entry['ra']), float(entry['dec']))
            exclusions = [flag for flag in EXCLUDE if flags.get(flag)]
            if exclusions:
                rejected['catalog_host_excluded'] += 1
                continue
            # Absolutely no measured peak, scatter, magnitude, or variability cut.
            selected.append({key: record[key] for key in ('tic', 'order', 'path', 'url', 'sha256', 'metadata')})
            if len(selected) == count:
                break
    db.rollback()
    db.close()
    if len(selected) != count:
        raise RuntimeError('Too few finished stars in the snapshot; wait for meaningful survey progress')
    rng = np.random.default_rng(SEED)
    for rank, entry in enumerate(selected):
        entry['sample_rank'] = rank
        if rank < paired_stars:
            entry['injection_period_days'] = float(np.exp(rng.uniform(np.log(.5), np.log(12.))))
            entry['injection_epoch_phase_fraction'] = float(rng.uniform())
    sample = {'created_utc': utc(), 'seed': SEED, 'source_sqlite': str(DB),
              'survey_status_counts_at_selection': counts, 'count_requested': count,
              'injection_stars': paired_stars, 'method_sha256': source_hash, 'config': CONFIG,
              'selection_sha256': sha(ROOT / 'data/broad_inventory/qlp_s0001_selection.csv.gz'),
              'selection_rule': 'First successfully screened targets in frozen pre-flux randomized survey order, excluding catalog host flags only; no flux/SNR/variability/magnitude selection.',
              'excluded_flags': list(EXCLUDE), 'selection_rejections': dict(rejected),
              'catalog_snapshot_hashes': {str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT / 'data/broad_catalogs').glob('*')) if p.is_file()},
              'surrogate_methods': list(METHODS), 'surrogate_block_days': .5,
              'injection_duration_days': .1, 'injection_depths_ppm': [1000, 3000],
              'injection_exposure_days': 1 / 48, 'targets': selected}
    save(path, sample)
    return sample


def load_raw(entry):
    if sha(entry['path']) != entry['sha256']:
        raise ValueError('Raw product SHA-256 mismatch; preserve and investigate')
    with fits.open(entry['path'], memmap=False) as hd:
        h = hd[1].header
        if (h.get('TIMESYS', hd[0].header.get('TIMESYS')) != 'TDB'
                or abs(h.get('BJDREFI', 0) + h.get('BJDREFF', 0) - 2457000) > 1e-8
                or abs(h.get('TIMEDEL', 0) - 1 / 48) > 1e-8):
            raise ValueError('Unexpected time reference or exposure; do not assume the injection timing')
        data = hd[1].data
        return [np.array(data[column], dtype=float if column != 'QUALITY' else int)
                for column in ('TIME', 'KSPSAP_FLUX', 'KSPSAP_FLUX_ERR', 'QUALITY')]


def surrogate(t, flux, error, quality, rng, method):
    """Shuffle paired residual/errors or randomize signs in ~0.5-d blocks.

    Original time stamps, exposure, quality flags, gaps and local trends remain.
    Correlations inside a block survive, but boundaries/long-time correlations
    do not; real unrecognized transits/variables can remain in this mixture.
    """
    out_flux, out_error = flux.copy(), error.copy()
    valid = np.isfinite(t) & np.isfinite(flux) & np.isfinite(error) & (flux > 0) & (error > 0) & (quality == 0)
    ids = np.flatnonzero(valid)
    ids = ids[np.argsort(t[ids])]
    dt = np.diff(t[ids])
    cadence = float(np.median(dt[dt > 0]))
    starts = np.r_[0, np.flatnonzero(dt >= max(1 / 48, 2.5 * cadence)) + 1, len(ids)]
    block_points = max(2, round(.5 / cadence))
    for a, b in zip(starts[:-1], starts[1:]):
        jj = ids[a:b]
        if len(jj) < 3:
            continue
        width = max(3, int(1 / cadence) // 2 * 2 + 1)
        trend = median_filter(flux[jj], size=width, mode='reflect')
        residual = flux[jj] / trend - 1
        normalized_error = error[jj] / trend
        blocks = [np.arange(k, min(k + block_points, len(jj))) for k in range(0, len(jj), block_points)]
        if method.startswith('block_permutation'):
            permutation = np.concatenate([blocks[k] for k in rng.permutation(len(blocks))])
            residual = residual[permutation]
            normalized_error = normalized_error[permutation]
        else:
            for block in blocks:
                residual[block] *= rng.choice([-1, 1])
        out_flux[jj] = (1 + residual) * trend
        out_error[jj] = normalized_error * trend
    return out_flux, out_error, block_points


def strong_local(peak):
    return bool(peak.get('screening_pass') and (peak.get('local_snr') or 0) >= 8
                and peak.get('n_events', 0) >= 3 and peak.get('n_positive_events', 0) >= 3
                and peak.get('n_positive_events', 0) >= .7 * peak.get('n_events', 0)
                and peak.get('event_snr_max_fraction') is not None
                and peak['event_snr_max_fraction'] <= .65 and peak.get('depth', 0) > 0)


def compact(scan):
    peaks = scan.get('peaks', [])
    return {'status': scan['status'], 'n_searched': scan['n_searched'],
            'period_grid_count': scan.get('search_metadata', {}).get('period_grid_count'),
            'max_diagonal_snr': max([p['snr'] for p in peaks], default=0.),
            'max_local_snr': max([p.get('local_snr') or 0 for p in peaks], default=0.),
            'any_initial_guard_pass_snr_ge_8': any(p['screening_pass'] and p['snr'] >= 8 for p in peaks),
            'any_strong_local_guard_pass': any(strong_local(p) for p in peaks),
            'peaks': peaks}


def inject(t, flux, period, epoch, depth):
    duration, exposure = .1, 1 / 48
    distance = phase(t, period, epoch)
    overlap = np.maximum(0, np.minimum(distance + exposure / 2, duration / 2)
                         - np.maximum(distance - exposure / 2, -duration / 2)) / exposure
    return flux * (1 - depth * overlap)


def recovery(scan, period, epoch):
    exact, alias = [], []
    for peak in scan.get('peaks', []):
        ratio = peak['period'] / period
        direct = abs(ratio - 1) < .01
        harmonic = any(abs(ratio / r - 1) < .01 for r in (.5, 2., 1 / 3, 3.))
        phase_match = abs(float(phase(np.array([peak['t0_btjd']]), period, epoch)[0])) <= .5 * (.1 + peak['duration_days']) + 1 / 48
        if direct and phase_match:
            exact.append(peak)
        if (direct or harmonic) and phase_match:
            alias.append(peak)
    return {'period_only_recovered': any(abs(p['period'] / period - 1) < .01 for p in scan.get('peaks', [])),
            'period_and_epoch_recovered': bool(exact), 'exact_initial_guards_recovered': any(p['screening_pass'] for p in exact),
            'exact_strong_local_recovered': any(strong_local(p) for p in exact),
            'period_epoch_or_simple_alias_recovered': bool(alias),
            'alias_strong_local_recovered': any(strong_local(p) for p in alias)}


def quantiles(values):
    clean = [float(v) for v in values if v is not None and np.isfinite(v)]
    return dict(zip(('minimum', 'p10', 'median', 'p90', 'maximum'), np.quantile(clean, [0, .1, .5, .9, 1]).tolist())) if clean else None


def summarize(sample):
    baseline, nulls, injections = [], [], []
    for entry in sample['targets']:
        path = CACHE / f"tic{entry['tic']}.json"
        if not path.exists():
            continue
        result = json.loads(path.read_text())
        baseline.append(result['baseline'])
        nulls.extend(result['nulls'])
        injections.extend(result['injections'])
    null_counts = Counter(n['method'] for n in nulls)
    summary = {k: v for k, v in sample.items() if k != 'targets'}
    summary.update({'updated_utc': utc(), 'targets': sample['targets'], 'completed_baselines': len(baseline),
                    'completed_surrogates': len(nulls), 'surrogates_by_method': dict(null_counts),
                    'completed_injections': len(injections), 'paid_calls': 0,
                    'sample_tmag': quantiles([b['tmag'] for b in baseline]),
                    'sample_residual_mad_ppm': quantiles([b['residual_mad_ppm'] for b in baseline]),
                    'sample_lag1_correlation': quantiles([b['lag1_correlation'] for b in baseline]),
                    'null_max_diagonal_snr': quantiles([n['max_diagonal_snr'] for n in nulls]),
                    'null_max_local_snr': quantiles([n['max_local_snr'] for n in nulls]),
                    'empirical_best_possible_plus_one_resolution': 1 / (len(nulls) + 1),
                    'planned_survey_160000_bonferroni_0_05_threshold': .05 / 160000,
                    'global_false_alarm_probability': None,
                    'null_rates': {key: {'count': sum(bool(n[key]) for n in nulls), 'trials': len(nulls),
                                       'fraction': sum(bool(n[key]) for n in nulls) / len(nulls) if nulls else None}
                                   for key in ('any_initial_guard_pass_snr_ge_8', 'any_strong_local_guard_pass')},
                    'baseline_rates': {key: {'count': sum(bool(b[key]) for b in baseline), 'trials': len(baseline)}
                                       for key in ('any_initial_guard_pass_snr_ge_8', 'any_strong_local_guard_pass')},
                    'injection_results': injections,
                    'limitations': ['Known-host exclusions are catalog based, not flux based; unknown astrophysical signals can enter nulls.',
                                    'Successful-screen snapshot and known-host exclusions change applicability to all survey targets.',
                                    'Block/sign surrogates preserve short correlations approximately but change long correlations and boundaries.',
                                    'Four trials within each of 128 stars are dependent; descriptive trial fractions are not calibrated FAPs.',
                                    'Only top two coarse peaks are searched for recovery; paired box injections use a fixed 2.4h duration.',
                                    'No whole-survey significance can be established at this empirical resolution.']})
    summary['software_versions'] = {'python': platform.python_version(), 'numpy': np.__version__,
                                    'astropy': astropy.__version__, 'scipy': scipy.__version__}
    summary['assessment_source_sha256'] = sha(Path(__file__))
    summary['catalog_provenance_verified_at_summary'] = CatalogIndex().provenance
    summary['catalog_provenance_verification_utc'] = utc()
    summary['null_maxima_examples'] = sorted(
        [{'tic': n['tic'], 'method': n['method'], 'max_diagonal_snr': n['max_diagonal_snr'],
          'max_local_snr': n['max_local_snr'], 'local_guard_pass': n['any_strong_local_guard_pass']}
         for n in nulls], key=lambda n: -n['max_diagonal_snr'])[:10]
    summary['null_rates_by_star'] = {key: {'stars_with_any_of_four': sum(
        any(n[key] for n in nulls if n['tic'] == b['tic']) for b in baseline),
        'stars': len(baseline)} for key in ('any_initial_guard_pass_snr_ge_8', 'any_strong_local_guard_pass')}
    summary['limitations'].append('Sign reversal can convert positive flares/outliers normally removed by preprocessing into artificial dips; extreme surrogate maxima must not define detection thresholds.')
    rates = {}
    for depth in (1000, 3000):
        rows = [r for r in injections if r['depth_ppm'] == depth]
        keys = ('period_only_recovered', 'period_and_epoch_recovered', 'exact_initial_guards_recovered',
                'exact_strong_local_recovered', 'period_epoch_or_simple_alias_recovered', 'alias_strong_local_recovered')
        rates[str(depth)] = {'trials': len(rows), **{k: sum(bool(r[k]) for r in rows) for k in keys}}
    summary['injection_recovery_counts_by_depth_ppm'] = rates
    strata = {'period_days': [(0., 3.), (3., 6.), (6., 12.01)],
              'residual_mad_ppm': [(0., 1000.), (1000., 3000.), (3000., float('inf'))],
              'tmag': [(0., 11.), (11., 12.), (12., 20.)],
              'observed_injected_cycles_before_edge_cuts': [(0., 3.), (3., 6.), (6., float('inf'))]}
    stratified = {}
    for name, bins in strata.items():
        stratified[name] = {}
        for lower, upper in bins:
            key = f'[{lower},{upper})'
            stratified[name][key] = {}
            for depth in (1000, 3000):
                rows = [r for r in injections if r['depth_ppm'] == depth and r.get(name) is not None
                        and lower <= r[name] < upper]
                stratified[name][key][str(depth)] = {'trials': len(rows),
                    'period_epoch_recovered': sum(r['period_and_epoch_recovered'] for r in rows),
                    'exact_local_recurrence_recovered': sum(r['exact_strong_local_recovered'] for r in rows)}
    summary['injection_recovery_stratified'] = stratified
    for method in METHODS:
        rows = [r for r in nulls if r['method'] == method]
        summary.setdefault('null_rates_by_method', {})[method] = {'trials': len(rows),
            'strong_local_count': sum(bool(r['any_strong_local_guard_pass']) for r in rows),
            'initial_guard_snr8_count': sum(bool(r['any_initial_guard_pass_snr_ge_8']) for r in rows)}
    recovery_table = '| Depth ppm | Injection pairs | Period + epoch | With initial guards | With local recurrence guards |\n|---:|---:|---:|---:|---:|\n'
    for depth, row in rates.items():
        recovery_table += f"| {depth} | {row['trials']} | {row['period_and_epoch_recovered']} | {row['exact_initial_guards_recovered']} | {row['exact_strong_local_recovered']} |\n"
    low = stratified['residual_mad_ppm']['[0.0,1000.0)']['1000']
    high1 = stratified['residual_mad_ppm']['[3000.0,inf)']['1000']
    high3 = stratified['residual_mad_ppm']['[3000.0,inf)']['3000']
    long1 = stratified['period_days']['[6.0,12.01)']['1000']
    long3 = stratified['period_days']['[6.0,12.01)']['3000']
    # All per-target surrogate records remain in the independent cache; metrics
    # and raw product hashes/retrieval URLs are stored in the repository report.
    save(ROOT / 'results/broad/noise_assessment.json', summary)
    text = f'''# Real-observation screening calibration

Updated {summary['updated_utc']}. Frozen broad-screen source SHA-256:
`{sample['method_sha256']}`. Screening config is exactly `{CONFIG}`.
No active survey source, database, or selection was modified.

The first {sample['count_requested']} successfully screened, catalog-unflagged
stars in the pre-flux randomized survey order were selected before inspecting
noise outcomes; no S/N, variability, magnitude, or radius selection was applied.
Excluded host flags: {', '.join(EXCLUDE)}. Snapshot selection exclusions and
catalog hashes are in `results/broad/noise_assessment.json`. Unknown variable,
eclipsing, and planetary signals may remain. Selection conditional on successful
screening and catalog exclusions limits extrapolation to all 160,000 stars.

Completed {len(baseline)} original real-light-curve screens, {len(nulls)}
surrogate screens, and {len(injections)} raw-flux injections. TESS magnitude
range/quantiles: `{summary['sample_tmag']}`; detrended residual MAD in ppm:
`{summary['sample_residual_mad_ppm']}`; adjacent-point correlation:
`{summary['sample_lag1_correlation']}`.

Each target has two within-gap-segment block permutations and two random block
sign reversals. Residuals/errors are transformed relative to the original 1-day
median trend and reconstructed in raw flux units; the frozen pipeline then
detrends and searches them afresh. Blocks are approximately half a day (12
points at 30-minute cadence). Permutations move paired residuals/errors;
sign flips retain errors. Original timestamps, quality flags, gaps, and nominal
exposures are retained. Blocks preserve internal short correlations but create
boundaries and alter longer variability, so these are approximate surrogate
noise/astrophysical-mixture tests, not exchangeable null realizations.

Initial-guard plus diagonal S/N >= 8 lead rate: `{summary['null_rates']['any_initial_guard_pass_snr_ge_8']}`.
Stronger local recurrence guard rate: `{summary['null_rates']['any_strong_local_guard_pass']}`.
The stronger guard requires local diagonal S/N >= 8, >=3 sampled and positive
events, >=70% positive events, <=65% event-power dominance, positive depth,
and the original guards. It omits catalog and stellar-physics exclusions.
The BLS-S/N and local-S/N sets are not nested: the local guard does not impose
a separate BLS S/N threshold, consistent with the survey prioritizer.
Rates count at least one of the top two peaks in each surrogate, rather than
counting every trial period as independent. They are empirical lead fractions,
not candidate false-positive probabilities. Repeated trials per star are
dependent. Null maxima are `{summary['null_max_diagonal_snr']}`.
Star-level counts with any of four transformed trials passing are
`{summary['null_rates_by_star']}`. Sign reversal can turn a positive flare/outlier
normally removed by preprocessing into an artificial dip: the largest maxima
(TIC 214299961) fail the recurrence guards, and must not define a threshold.
The transformed distributions are not exchangeable with the original data;
their rates quantify screening vulnerability to structured/asymmetric noise.

The first {sample['injection_stars']} selected stars each receive matched
1000 and 3000 ppm box injections with 2.4-hour duration, independent frozen
random phase, and log-uniform period 0.5--12 days. The box is integrated across
the 30-minute exposure before multiplying the actual raw KSPSAP flux. No gaps
are filled; existing noise, contamination, variability, and exposure loss remain.
Error arrays are unchanged. Injection parameters are fixed in the sample file
before trial outcomes. Recovery requires period within 1% and epoch agreement
within half the sum of injected/recovered durations plus one exposure; separate
period-only, simple-alias, initial-guard, and stronger-local metrics are retained.

Recovery counts by depth (denominators are explicit):

{recovery_table}

Recovery is also stratified in the JSON by fixed descriptive period boundaries
3/6 days, residual-MAD boundaries 1000/3000 ppm, TESS magnitude boundaries
11/12, and observed-event-count boundaries 3/6. Small strata are descriptive;
they are not a general occurrence-rate or completeness calibration. The
stratification was added after the first aggregate results, for descriptive
heterogeneity assessment; it was not used to select targets or injection pairs.
The
exposure duration is verified from the QLP TIMEDEL header as 1/48 day.

In this fixed sample, for residual MAD below 1000 ppm, 1000-ppm period/epoch
recovery is {low['period_epoch_recovered']}/{low['trials']}
({low['exact_local_recurrence_recovered']}/{low['trials']} also meet local recurrence).
Above 3000-ppm residual MAD, recovery is {high1['period_epoch_recovered']}/{high1['trials']}
at 1000 ppm and {high3['period_epoch_recovered']}/{high3['trials']} at 3000 ppm.
For periods longer than 6 days, exact period/epoch recovery is
{long1['period_epoch_recovered']}/{long1['trials']} at 1000 ppm and
{long3['period_epoch_recovered']}/{long3['trials']} at 3000 ppm;
{long1['exact_local_recurrence_recovered']} and {long3['exact_local_recurrence_recovered']}, respectively,
pass the local recurrence cuts. These small descriptive groups show that coarse
screening is strongly incomplete for weak, noisy, and long-period signals.
Only two top peaks are retained, so stronger existing variability can displace
injections. A negative survey screen cannot rule out planets in these regimes;
multi-sector data, multiple detrending scales, and searches retaining more
independent peaks are scientifically justified next steps.

Even an ideal 512-trial empirical null could resolve a plus-one tail probability
only to 1/513 ~= 0.00195; the 160,000-star Bonferroni 5% target is 3.125e-7.
This experiment cannot establish survey-global significance, and trial/star
dependence and null mismatch make that bound more optimistic than justified.
Candidate-specific raw-flux, independent-sector, pixel, contamination, and
catalog/literature checks remain necessary. No discovery or global FAP is claimed.

Run/resume with `/workspace/openrouter-check-venv/bin/python src/broad_noise.py`.
Independent atomic checkpoints: `{CACHE}`. Raw files are reused and SHA-256
checked against the survey records; exact target URLs/hashes and input catalog
hashes are retained in the repository JSON. Nothing is downloaded or deleted.
Paid external requests: zero.
'''
    (ROOT / 'reports/broad_noise.md').write_text(text)
    return summary


def run(args):
    sample = frozen_sample(args.targets, args.injection_stars)
    for entry in sample['targets']:
        path = CACHE / f"tic{entry['tic']}.json"
        if path.exists():
            continue
        t, f, e, q = load_raw(entry)
        baseline = compact(screen_curve(t, f, e, q, config=CONFIG))
        pp = preprocess_curve(t, f, e, q, config=FastScreenConfig(**CONFIG))
        residual = pp['f'] - 1
        adjacent = (np.diff(pp['t']) < 1.5 * (pp['cadence_days'] or 1 / 48))
        correlation = np.corrcoef(residual[:-1][adjacent], residual[1:][adjacent])[0, 1] if adjacent.sum() > 3 else np.nan
        baseline.update({'tic': entry['tic'], 'tmag': entry['metadata'].get('TESSMAG'),
                         'residual_mad_ppm': float(1.4826 * np.median(np.abs(residual - np.median(residual))) * 1e6),
                         'lag1_correlation': float(correlation) if np.isfinite(correlation) else None})
        nulls = []
        for i, method in enumerate(METHODS):
            rng = np.random.default_rng(np.random.SeedSequence([SEED, entry['tic'], i]))
            nf, ne, block_points = surrogate(t, f, e, q, rng, method)
            scan = compact(screen_curve(t, nf, ne, q, config=CONFIG))
            nulls.append({'tic': entry['tic'], 'method': method, 'block_points': block_points, **scan})
        injections = []
        if entry['sample_rank'] < args.injection_stars:
            good = np.isfinite(t) & np.isfinite(f) & (q == 0)
            period = entry['injection_period_days']
            epoch = float(t[good].min() + period * entry['injection_epoch_phase_fraction'])
            for depth_ppm in (1000, 3000):
                injected = inject(t, f, period, epoch, depth_ppm / 1e6)
                scan = screen_curve(t, injected, e, q, config=CONFIG)
                # Positive exposure-integrated overlap before quality/edge cuts.
                physical_in = good & (np.abs(phase(t, period, epoch)) < .1 / 2 + 1 / 96)
                cycles = np.unique(np.round((t[physical_in] - epoch) / period).astype(int))
                injections.append({'tic': entry['tic'], 'depth_ppm': depth_ppm, 'period_days': period,
                    'epoch_btjd': epoch, 'duration_days': .1, 'observed_injected_cycles_before_edge_cuts': len(cycles),
                    'tmag': baseline['tmag'], 'residual_mad_ppm': baseline['residual_mad_ppm'],
                    **recovery(scan, period, epoch), 'screen': compact(scan)})
        save(path, {'tic': entry['tic'], 'raw_sha256': entry['sha256'], 'method_sha256': sample['method_sha256'],
                    'checked_utc': utc(), 'baseline': baseline, 'nulls': nulls, 'injections': injections})
        if (entry['sample_rank'] + 1) % 16 == 0:
            summary = summarize(sample)
            print('CHECKPOINT', summary['completed_baselines'], summary['completed_surrogates'], summary['completed_injections'], flush=True)
    summary = summarize(sample)
    print('FINAL', json.dumps({key: summary[key] for key in ('completed_baselines', 'completed_surrogates',
          'completed_injections', 'null_rates', 'injection_recovery_counts_by_depth_ppm')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--targets', type=int, default=128)
    parser.add_argument('--injection-stars', type=int, default=64)
    run(parser.parse_args())
