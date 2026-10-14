"""Local-depth positive/negative controls on three real cadence patterns.

These controls test local estimators; they are not completeness or global FAP.
No original observations are modified and no paid API is used.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from broad_vet import integrated_trapezoid, local_quadratic, phase

OUT = ROOT / 'results/stellar_systems'
SEED = 20261014


def measure(t, f, e, p, epoch, duration):
    groups = [local_quadratic(t, f, e, 2 * p, epoch + k * p, duration)[0]
              for k in (0, 1)]
    if any(group['depth'] is None or group['n_events'] < 2 for group in groups):
        return {'status': 'insufficient_event_groups'}
    contrast = groups[0]['depth'] - groups[1]['depth']
    error = np.hypot(groups[0]['depth_err'], groups[1]['depth_err'])
    return {'status': 'measured', 'events': [group['n_events'] for group in groups],
            'depths': [group['depth'] for group in groups],
            'contrast': float(contrast), 'contrast_error': float(error),
            'conditional_sigma': float(abs(contrast) / error)}


def main():
    rng = np.random.default_rng(SEED)
    results = []
    for tic in [294272014, 38825979, 281581467]:
        original = json.loads((OUT / f'scout_tic{tic}.json').read_text())
        path = Path(original['raw_path'])
        if hashlib.sha256(path.read_bytes()).hexdigest() != original['raw_sha256']:
            raise ValueError('Raw observation changed')
        with fits.open(path, memmap=False) as hd:
            data = hd[1].data
            t = np.asarray(data['TIME'], float)
            use = np.isfinite(t) & (data['QUALITY'] == 0)
            t = np.sort(t[use])
            exposure = float(hd[1].header['TIMEDEL'])
        p, epoch, duration = [original['refined'][name] for name in
                              ['period_days', 'epoch_btjd', 'duration_days']]
        # Small deterministic slow variability gives a realistic local
        # interpolation challenge without pretending to emulate all raw noise.
        baseline = 1 + .001 * np.sin(2 * np.pi * (t - t[0]) / 7.3)
        patterns = [integrated_trapezoid(phase(t, 2 * p, epoch + k * p),
                                         duration, duration * .15, exposure)
                    for k in (0, 1)]
        e = np.full(len(t), .0005)
        trials = []
        for label, depths in [('equal', [.03, .03]), ('unequal', [.03, .02])]:
            clean = baseline * (1 - sum(depth * shape for depth, shape in zip(depths, patterns)))
            no_noise = measure(t, clean, e, p, epoch, duration)
            draws = [measure(t, clean + rng.normal(0, e), e, p, epoch, duration)
                     for _ in range(20)]
            trials.append({'label': label, 'injected_central_depths': depths,
                           'no_noise_measurement': no_noise, 'noise_draws': draws})
        equal = [draw['conditional_sigma'] for draw in trials[0]['noise_draws']
                 if draw['status'] == 'measured']
        unequal = [draw['conditional_sigma'] for draw in trials[1]['noise_draws']
                   if draw['status'] == 'measured']
        results.append({'tic': tic, 'raw_sha256': original['raw_sha256'],
                        'period_days': p, 'duration_days': duration, 'exposure_days': exposure,
                        'n_cadences': len(t), 'trials': trials,
                        'equal_depth_max_sigma': max(equal),
                        'unequal_depth_min_sigma': min(unequal)})
    payload = {'seed': SEED, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'controls': results,
               'limitations': 'Synthetic independent Gaussian noise on real timestamps; only local estimator checks, not realistic search-null calibration, completeness or global false-alarm probabilities.'}
    (OUT / 'scout_controls.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    print(json.dumps([{key: result[key] for key in ['tic', 'equal_depth_max_sigma',
                                                    'unequal_depth_min_sigma']}
                       for result in results], indent=2))


if __name__ == '__main__':
    main()
