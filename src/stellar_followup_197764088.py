"""Independent multi-sector raw eclipse measurements for TIC197764088.

The target is an aperture label until the pixel source test identifies the
variable star. Timing fits explicitly preserve independent sector ephemerides.
No frozen broad-source or raw observation is changed.
"""
from pathlib import Path
import hashlib
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from broad_vet import local_quadratic, fit_shape, fit_event_timings, mad, phase
from archive import ROOT, verify_fits_checksums
from stellar_scout import atomic_json

TIC = 197764088
OUT = ROOT / 'candidates/stellar_systems/tic197764088'
SEED = 2.682059291551926


def read(record, column):
    path = Path(record['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
        raise ValueError('Raw file changed')
    with fits.open(path, memmap=False) as hd:
        verify_fits_checksums(path, hd)
        h, data = hd[1].header, hd[1].data
        if h.get('TIMESYS') != 'TDB' or h.get('BJDREFI') != 2457000:
            raise ValueError('Unexpected timing')
        t = np.asarray(data['TIME'], float)
        actual_column = column if column in data.names else column.replace('KSPSAP_', 'DET_')
        error_column = 'KSPSAP_FLUX_ERR' if 'KSPSAP_FLUX_ERR' in data.names else 'DET_FLUX_ERR'
        f = np.asarray(data[actual_column], float)
        e = np.asarray(data[error_column], float)
        use = (np.isfinite(t) & np.isfinite(f) & (data['QUALITY'] == 0) & (f > 0))
        t, f, e = t[use], f[use], e[use]
        dt = float(h['TIMEDEL'])
    order = np.argsort(t)
    t, f, e = t[order], f[order] / np.median(f), e[order] / np.median(f)
    # Recent QLP products sometimes contain all-zero formal DET_FLUX_ERR.
    # Retain their observations and substitute a documented empirical error.
    if np.any(~np.isfinite(e) | (e <= 0)):
        empirical = max(mad(np.diff(f)) / np.sqrt(2), 1e-8)
        e = np.where(np.isfinite(e) & (e > 0), e, empirical)
    return t, f, e, dt


def refine(t, f, e, dt):
    trend = median_filter(f, size=max(5, int(3 / dt) // 2 * 2 + 1), mode='nearest')
    ff = f / trend
    noise = max(mad(np.diff(ff)) / np.sqrt(2), np.median(e / trend))
    ee = np.maximum(e / trend, noise)
    grid = np.linspace(SEED * .99, SEED * 1.01, 1001)
    power = BoxLeastSquares(t, ff, ee).power(grid, [.04, .06, .08, .10], oversample=15, objective='snr')
    j = int(np.argmax(power.power))
    return {'period_days': float(power.period[j]), 'epoch_btjd': float(power.transit_time[j]),
            'duration_days': float(power.duration[j]), 'depth': float(power.depth[j]),
            'diagonal_snr': float(power.depth_snr[j]), 'grid_boundary': j in (0, len(grid) - 1)}


def analyze(record):
    sector = record['sector']
    t, f, e, dt = read(record, 'KSPSAP_FLUX')
    refined = refine(t, f, e, dt)
    p, t0, duration = [refined[key] for key in ['period_days', 'epoch_btjd', 'duration_days']]
    arrays, columns = {}, {}
    for column in ['SAP_FLUX', 'KSPSAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']:
        tc, fc, ec, dt = read(record, column)
        primary, local = local_quadratic(tc, fc, ec, p, t0, duration)
        secondary, slocal = local_quadratic(tc, fc, ec, p, t0 + p / 2, duration)
        shape = fit_shape(local, p, t0, duration, dt)
        timing = None
        if local is not None and shape is not None:
            local['exposure'] = np.full(len(local['t']), dt)
            timing = fit_event_timings([local], [sector], p, t0, shape)
        groups = [local_quadratic(tc, fc, ec, 2 * p, t0 + k * p, duration)[0]
                  for k in (0, 1)]
        contrast = None
        if all(group['n_events'] >= 2 for group in groups):
            difference = groups[0]['depth'] - groups[1]['depth']
            error = np.hypot(groups[0]['depth_err'], groups[1]['depth_err'])
            contrast = {'difference': float(difference), 'error': float(error),
                        'conditional_sigma': float(abs(difference) / error)}
        columns[column] = {'primary': primary, 'half_phase_secondary': secondary,
                           'primary_shape': shape, 'timing': timing,
                           'parity_groups': groups, 'parity_contrast': contrast}
        if column == 'SAP_FLUX':
            arrays = {'time_btjd': tc, 'sap_flux': fc, 'sap_error': ec,
                      'local_primary_t': np.array([]) if local is None else local['t'],
                      'local_primary_f': np.array([]) if local is None else local['f'],
                      'local_secondary_t': np.array([]) if slocal is None else slocal['t'],
                      'local_secondary_f': np.array([]) if slocal is None else slocal['f']}
    np.savez_compressed(OUT / f'sector{sector}_photometry.npz', **arrays)
    result = {'tic': TIC, 'sector': sector, 'exposure_days': dt, 'refined': refined,
              'columns': columns, 'raw_observation': record,
              'column_aliases': 'QLP sectors67 onward use DET_FLUX / DET_FLUX_ERR / DET_FLUX_SML / DET_FLUX_LAG for old KSPSAP labels; SAP_FLUX is unchanged.',
              'missing_error_handling': 'Missing/nonpositive formal DET_FLUX_ERR replaced with robust first-difference scatter/sqrt(2), before independent local baseline/scatter floors. Sectors95/105 formal errors are all NaN.'}
    atomic_json(OUT / f'sector{sector}_measurements.json', result)
    return result


def linear_fit(events, pseed, epoch):
    frame = pd.DataFrame(events)
    if len(frame) < 4:
        return None
    if not np.isfinite(pseed) or pseed <= 0 or not np.isfinite(epoch):
        raise ValueError('Invalid global timing seed')
    cycles = np.rint((frame.time_btjd.to_numpy() - epoch) / pseed).astype(int)
    times, errors = frame.time_btjd.to_numpy(), frame.time_error_days.to_numpy()
    ref = float(cycles.mean())
    a = np.column_stack([np.ones(len(frame)), cycles - ref])
    aw = a / errors[:, None]
    beta = np.linalg.lstsq(aw, times / errors, rcond=None)[0]
    residual = times - a @ beta
    chi = float(np.sum((residual / errors) ** 2))
    cov = np.linalg.pinv(aw.T @ aw) * max(1, chi / (len(frame) - 2))
    frame['global_cycle'] = cycles
    frame['oc_days'] = residual
    frame.to_csv(OUT / 'primary_event_timings.csv', index=False)
    transform = np.array([1, -ref])
    return {'n_events': len(frame), 'period_days': float(beta[1]),
            'period_error_days_scatter_scaled': float(np.sqrt(cov[1, 1])),
            'epoch_btjd': float(beta[0] - beta[1] * ref),
            'epoch_error_days_scatter_scaled': float(np.sqrt(transform @ cov @ transform)),
            'rms_oc_minutes': float(np.sqrt(np.mean(residual ** 2)) * 1440),
            'chi2': chi, 'dof': len(frame) - 2, 'cycle_seed_period_days': pseed,
            'limitations': 'Conditional integer-cycle assignment; covariance scatter inflated, source/blend/PRF systematics not included.'}


def main():
    observations = json.loads((OUT / 'observations.json').read_text())
    records = [analyze(record) for record in observations if record['status'] == 'retained']
    rows, events = [], []
    for record in records:
        raw = record['columns']['SAP_FLUX']
        timing = raw['timing'] or {}
        events.extend(timing.get('events', []))
        shape = raw['primary_shape'] or {}
        secondary = raw['half_phase_secondary']
        rows.append({'sector': record['sector'], 'period_days_independent': record['refined']['period_days'],
                     'primary_depth_ppm': (raw['primary']['depth'] or 0) * 1e6,
                     'primary_error_ppm': (raw['primary'].get('depth_err') or 0) * 1e6,
                     'primary_n_events': raw['primary']['n_events'],
                     'secondary_depth_ppm': (secondary['depth'] or 0) * 1e6,
                     'secondary_error_ppm': (secondary.get('depth_err') or 0) * 1e6,
                     'secondary_n_events': secondary['n_events'],
                     'shape_central_depth': shape.get('depth'),
                     'shape_duration_days': shape.get('duration_days'),
                     'shape_ingress_fraction': shape.get('ingress_fraction'),
                     'timing_period_days_local': (timing.get('fitted_ephemeris') or {}).get('period')})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / 'sector_depths.csv', index=False)
    pseed = float(np.median([row['timing_period_days_local'] for row in rows
                             if row['timing_period_days_local'] is not None]))
    epoch = records[0]['refined']['epoch_btjd']
    global_fit = linear_fit(events, pseed, epoch)
    result = {'tic': TIC, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'n_retained_sectors': len(records), 'sectors': [r['sector'] for r in records],
              'primary_ephemeris': global_fit, 'sector_measurements': rows,
              'assessment': 'Aperture contains strong two-eclipse periodic variability. Pixel localization and independent prior checks required before host or novelty assignment.',
              'timing_exclusions': 'QLP v01 sectors102–104 omitted because published corrections are needed; these omissions do not assert lack of detector coverage.'}
    atomic_json(OUT / 'multisector_summary.json', result)
    fig, axes = plt.subplots(2, 1, figsize=(11, 8))
    for record in records:
        curves = np.load(OUT / f"sector{record['sector']}_photometry.npz")
        p, t0 = record['refined']['period_days'], record['refined']['epoch_btjd']
        axes[0].plot(phase(curves['local_primary_t'], p, t0) * 24, curves['local_primary_f'],
                     '.', ms=2, alpha=.3, label=f"S{record['sector']}")
        axes[1].plot(phase(curves['local_secondary_t'], p, t0 + p / 2) * 24,
                     curves['local_secondary_f'], '.', ms=2, alpha=.3)
    for ax, label in zip(axes, ['Primary', 'Half-phase eclipse']):
        ax.set_xlim(-3, 3)
        ax.set_ylabel('Local SAP relative flux')
        ax.set_title(label)
    axes[0].legend(ncol=7, markerscale=3)
    axes[1].set_xlabel('Hours from independently fitted sector eclipse center')
    fig.suptitle('TIC197764088 aperture: source assignment requires pixels')
    fig.tight_layout()
    fig.savefig(OUT / 'multisector_eclipses.png', dpi=160)
    plt.close(fig)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
