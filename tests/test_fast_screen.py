"""Meaningful coverage tests for the fast broad survey search."""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import time
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from fast_screen import FastScreenConfig, known_mask, period_grid, phase, preprocess_curve, screen_curve, screen_fits


class FastScreenTests(unittest.TestCase):
    def test_log_grid_bounds_nearest_epoch_drift(self):
        cfg = FastScreenConfig()
        baseline = 27.
        grid = period_grid(baseline, cfg)
        self.assertLess(len(grid), 5000)
        for period in np.geomspace(.3001, 12.999, 1000):
            nearest = grid[np.argmin(abs(grid - period))]
            drift = baseline * abs(nearest / period - 1)
            self.assertLessEqual(drift, min(cfg.durations) * 1.01)
        # Interior points have the half-duration nearest-neighbor bound.
        periods = (grid[:-1] + grid[1:]) / 2
        drift = baseline * np.abs(grid[:-1] / periods - 1)
        self.assertLessEqual(drift.max(), min(cfg.durations) * .501)

    def test_gap_safe_detrending_at_thirty_and_sixty_minutes(self):
        a = np.arange(0, 2., 2 / 1440)
        b = np.arange(a[-1] + 31 / 1440, 4., 2 / 1440)
        c = np.arange(b[-1] + 61 / 1440, 6., 2 / 1440)
        t = np.r_[a, b, c]
        f = np.r_[np.full(len(a), 100.), np.full(len(b), 500.), np.full(len(c), 900.)]
        out = preprocess_curve(t, f, f * .001)
        self.assertEqual(out['n_segments'], 3)
        np.testing.assert_allclose(out['f'], 1.)
        self.assertTrue(np.all(out['e'] > 0))

    def test_quality_and_known_ephemerides_are_removed(self):
        t = np.arange(0, 4., 2 / 1440)[::-1]
        known = [dict(period=1., epoch_bjd=2457000.5, duration_hours=1.)]
        quality = np.zeros(len(t), dtype=int)
        quality[t > 3.5] = 1
        out = preprocess_curve(t, np.ones(len(t)), np.full(len(t), .001), quality=quality, known=known)
        self.assertGreater(out['n_known_masked'], 0)
        self.assertLess(out['t'].max(), 3.5)
        self.assertTrue(np.all(known_mask(out['t'], known)))

    def test_recovers_short_duration_and_long_period_injections(self):
        rng = np.random.default_rng(7654)
        t = np.arange(0, 27., 2 / 1440)
        for period, duration in ((.73123, .025), (3.4312, .05), (12.3491, .10)):
            f = 1 + rng.normal(0, .0004, len(t))
            f[np.abs(phase(t, period, 1.173)) < duration / 2] -= .0015
            out = screen_curve(t, f, np.full(len(t), .0004))
            matches = [v for v in out['peaks'] if abs(v['period'] / period - 1) < .002]
            self.assertTrue(matches, (period, out['peaks']))
            self.assertTrue(matches[0]['screening_pass'])
            self.assertGreater(matches[0]['snr'], 10)

    def test_single_event_is_preserved_and_cannot_pass_periodic_guards(self):
        rng = np.random.default_rng(411)
        t = np.arange(0, 27., 2 / 1440)
        f = 1 + rng.normal(0, .0004, len(t))
        f[np.abs(t - 12.123) < .05] -= .002
        out = screen_curve(t, f, np.full(len(t), .0004))
        self.assertTrue(any(abs(v['time'] - 12.123) < .05 for v in out['single_events']))
        self.assertFalse(any(v['screening_pass'] and v['snr'] > 10 for v in out['peaks']))

    def test_thirty_minute_cadence_does_not_become_all_gaps(self):
        rng = np.random.default_rng(66)
        t = np.arange(0, 27., 30 / 1440)
        # Barycentric timing perturbations must not split every cadence.
        t += np.sin(t) * 1e-6
        f = 1 + rng.normal(0, .0004, len(t))
        f[np.abs(phase(t, 4.3721, 1.173)) < .1 / 2] -= .002
        result = screen_curve(t, f, np.full(len(t), .0004), config=dict(n_refine=0, n_peaks=2))
        self.assertEqual(result['n_segments'], 1)
        self.assertGreater(result['n_searched'], 1200)
        self.assertEqual(result['search_metadata']['minimum_duration_searched_days'], .05)
        self.assertTrue(any(abs(v['period'] / 4.3721 - 1) < .005 for v in result['peaks']))


def benchmark():
    """Reproduce the recorded controls; cached public observations are optional.

    Run: python tests/test_fast_screen.py --benchmark
    Exact reference FITS products are retained in data/observations.json and the
    broad archive manifest. This command never downloads or modifies raw data.
    """
    root = Path(__file__).resolve().parents[1]
    broad = dict(n_refine=0, n_peaks=2, compute_full_stats=False)
    rng = np.random.default_rng(20261013)
    rows = []
    nulls = []
    for trial in range(36):
        t = np.arange(0, 27., 30 / 1440)
        # Simulate the mid-sector downlink and scattered bad cadences.
        good = ~((t > 12.6) & (t < 13.8)) & (rng.random(len(t)) > .035)
        t = t[good]
        p = float(np.exp(rng.uniform(np.log(.4), np.log(12.5))))
        d = (.05, .1, .16)[trial % 3]
        epoch = float(rng.uniform(.3, min(p, 2.)))
        inside = np.abs(phase(t, p, epoch)) < d / 2
        requested = (8., 12., 20.)[(trial // 3) % 3]
        sigma = .0008
        depth = float(requested * sigma / np.sqrt(inside.sum()))
        f = 1 + rng.normal(0, sigma, len(t)) + .0006 * np.sin(t / 2.7)
        f[inside] -= depth
        out = screen_curve(t, f, np.full(len(t), sigma), config=broad)
        recovered = any(abs(v['period'] / p - 1) < .005 for v in out['peaks'])
        harmonic = any(min(abs(v['period'] / p - k) for k in (.5, 1., 2.)) < .005 for v in out['peaks'])
        rows.append(dict(trial=trial, period=p, duration=d, epoch=epoch, depth=depth, requested_white_snr=requested,
                         exact_recovered=recovered, exact_or_simple_alias=harmonic,
                         runtime_seconds=out['runtime_seconds'], peaks=out['peaks']))
    for trial in range(64):
        t = np.arange(0, 27., 30 / 1440)
        f = 1 + rng.normal(0, .0008, len(t))
        out = screen_curve(t, f, np.full(len(t), .0008), config=broad)
        nulls.append(dict(trial=trial, best_snr=out['peaks'][0]['snr'], runtime_seconds=out['runtime_seconds']))
    controls = []
    observations = json.loads((root / 'data/observations.json').read_text())
    for sector in (2, 5, 8):
        rec = next(v for v in observations if v['tic'] == 307210830 and v['sector'] == sector)
        if not Path(rec['local_path']).exists():
            continue
        out = screen_fits(rec['local_path'])
        controls.append(dict(tic=rec['tic'], sector=sector, data_uri=rec['data_uri'],
                             expected_period=3.6906764,
                             exact_recovered=any(abs(v['period'] / 3.6906764 - 1) < .005 for v in out['peaks']),
                             result=out))
    qlp = []
    for tic in (8195168, 8195216):
        path = Path(f'/workspace/tess-research-data/broad_probe/products/qlp_tic{tic}_s1.fits')
        if path.exists():
            out = screen_fits(path, column='KSPSAP_FLUX', config=broad)
            qlp.append(dict(tic=tic, sector=1, path=str(path), result=out))
    payload = dict(checked_utc=datetime.now(timezone.utc).isoformat(), seed=20261013,
                   broad_config=broad, injection_cases=rows, white_noise_nulls=nulls,
                   known_controls=controls, qlp_controls=qlp)
    outpath = root / 'results/fast_screen_validation.json'
    outpath.write_text(json.dumps(payload, indent=2, allow_nan=False))
    runtimes = [v['runtime_seconds'] for v in rows]
    lines = [
        '# Fast broad survey validation', '',
        f"Recorded {payload['checked_utc']}; deterministic RNG seed 20261013.", '',
        'Reproduce with `python tests/test_fast_screen.py --benchmark`; run integrity tests with '
        '`python -m unittest discover -s tests -p test_fast_screen.py`.', '',
        '## Search and cadence handling', '',
        'This searches actual light curves across 0.3–13 days with a logarithmic period grid. '
        'The nearest trial period drifts by at most half the shortest trial duration over the '
        'observed baseline, including both grid endpoints. Four nominal durations (0.025, 0.05, '
        '0.10, 0.16 days) are filtered to those >=1.5 effective cadences. For sector-1 QLP '
        '30-minute cadence this removes 0.025-day boxes; 0.05/0.10/0.16-day boxes remain.', '',
        'Quality-zero data are split at gaps >=max(30 minutes, 2.5 median cadences), binned '
        'without crossing gaps, detrended with a one-day robust median, and conservatively '
        'assigned errors at least as large as observed residual MAD. Edge exclusions of 0.06 '
        'days lose some observing completeness. Known ephemerides are masked before binning. '
        'Two sampled local events, two positive local depths, <=0.15 duty cycle, and no event '
        'contributing >85% of summed positive event S/N squared are required by the guard flag.', '',
        'Broad mode uses `n_refine=0, n_peaks=2, compute_full_stats=False`; it still evaluates '
        'the full documented period grid and retains an independent local single-event scan. '
        'Detailed followup uses the default top-20 local grid refinement and five detailed peaks.', '',
        '## Recovery and runtime', '',
        f"Thirty-minute synthetic injections: {sum(v['exact_recovered'] for v in rows)}/36 "
        f"periods recovered within 0.5% in the top two peaks; {sum(v['exact_or_simple_alias'] for v in rows)}/36 "
        'including simple half/double aliases. Injections span randomized 0.4–12.5-day periods, '
        '1.2/2.4/3.84-hour durations, requested ideal white-noise S/N 8/12/20, weak smooth '
        'variability, a 1.2-day downlink gap, and 3.5% missing cadences.', '',
        f"Synthetic broad-screen runtime: median {np.median(runtimes):.4f} s, "
        f"maximum {max(runtimes):.4f} s per star on this machine (one worker; excludes downloads).", '',
    ]
    for snr in (8., 12., 20.):
        subset = [v for v in rows if v['requested_white_snr'] == snr]
        lines.append(f"- Requested S/N {snr:g}: {sum(v['exact_recovered'] for v in subset)}/{len(subset)} exact-period recoveries.")
    lines += ['', f"L 98-59 c known positive control: {sum(v['exact_recovered'] for v in controls)}/{len(controls)} "
              'available two-minute sectors recovered near 3.6906764 days. This verifies sensitivity '
              'to a known planet; it is not a new candidate.']
    for v in qlp:
        r = v['result']
        lines += ['', f"QLP TIC {v['tic']}, sector 1: {r['n_input']} raw, {r['n_clean']} quality-filtered, "
                  f"{r['n_searched']} searched cadences, {r['search_metadata']['period_grid_count']} periods, "
                  f"{r['runtime_seconds']:.4f} s. Best score {r['peaks'][0]['snr']:.2f}; this benchmark "
                  'does not assign a planetary interpretation.']
    maxima = np.array([v['best_snr'] for v in nulls])
    lines += ['', '## Statistical limitations', '',
              f"64 white-noise nulls produced median/max best BLS S/N {np.median(maxima):.2f}/{maxima.max():.2f}. "
              'These are runtime/behavior controls, not a sufficient survey-wide false-alarm calibration. '
              'A many-star search requires empirical null controls preserving correlated/instrumental noise, '
              'trial accounting, and independent confirmation of recurrence and source. No reported BLS '
              'score is a false-alarm probability.', '',
              'This deliberately broad stage can miss shallow, shorter-than-resolution, single/double '
              'transit, boundary, very long-period, or strong-variability signals. A recovered harmonic '
              'must be resolved by event inspection. Failures at requested S/N 8 do not exclude planets. '
              'The injection test is small and synthetic; real-data completeness is not established.', '',
              'Every measurement and case is in `results/fast_screen_validation.json`.']
    (root / 'reports/fast_screen_validation.md').write_text('\n'.join(lines) + '\n')
    print(outpath)
    print(f"exact recovery {sum(v['exact_recovered'] for v in rows)}/36; median runtime {np.median(runtimes):.4f}s")


if __name__ == '__main__':
    if '--benchmark' in sys.argv:
        benchmark()
    else:
        unittest.main()
