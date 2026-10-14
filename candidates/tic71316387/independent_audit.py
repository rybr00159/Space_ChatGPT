"""Independent event-by-event odd/even audit of TIC 71316387.

Run from the repository root with the research Python environment. Public FITS
sources and their immutable SHA256 values are listed in observations.json.
No broad-survey routines or previously fitted eclipse measurements are used.
"""
from pathlib import Path
import hashlib
import json
import os
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from scipy.stats import t as student_t
os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
P = .6518369465110034
EPOCH = 2067.327192476
WIDTH = 1.46256649 / 24


def integrated_trapezoid(x, width, ingress_fraction, exposure):
    """Exact exposure integral; peak unconvolved depth equals one."""
    ingress = width * ingress_fraction
    positive = lambda z: np.maximum(z, 0.) ** 2
    def integral(z):
        return (positive(z + width / 2) - positive(z + width / 2 - ingress)
                - positive(z - width / 2 + ingress) + positive(z - width / 2)) / (2 * ingress)
    return np.clip((integral(x + exposure / 2) - integral(x - exposure / 2)) / exposure, 0, 1)


def prepare(column):
    events, sources = [], []
    for row in json.loads((HERE / 'observations.json').read_text()):
        if row['status'] != 'retained':
            continue
        path = Path(row['path'])
        if not path.exists():
            cache = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data'))
            path = cache / 'tic71316387' / 'observations' / path.name
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != row['sha256']:
            raise ValueError(f"SHA256 mismatch: {path}")
        with fits.open(path, memmap=False) as hd:
            data, header = hd[1].data, hd[1].header
            assert header['TIMESYS'] == 'TDB' and header['BJDREFI'] == 2457000
            processed = 'KSPSAP_FLUX' if 'KSPSAP_FLUX' in data.names else 'DET_FLUX'
            errcol = processed + '_ERR'
            col = column if column == 'SAP_FLUX' else processed
            time = np.asarray(data['TIME'], float)
            flux = np.asarray(data[col], float)
            error = np.asarray(data[errcol], float)
            error_fallback = None
            if not np.any(np.isfinite(error) & (error > 0)):
                noise_use = ((data['QUALITY'] == 0) & np.isfinite(time)
                             & np.isfinite(flux) & (flux > 0))
                diff = np.diff(flux[noise_use])
                estimate = float(1.4826 * np.median(np.abs(diff - np.median(diff))) / np.sqrt(2))
                error = np.full(len(flux), estimate)
                error_fallback = dict(method='All supplied flux errors are invalid; substituted robust first-difference scatter divided by sqrt(2), then retained local out-of-event floors.',
                                      estimated_flux_error=estimate)
            use = ((data['QUALITY'] == 0) & np.isfinite(time) & np.isfinite(flux)
                   & np.isfinite(error) & (flux > 0) & (error > 0))
            time, flux, error = time[use], flux[use], error[use]
            order = np.argsort(time)
            time, flux, error = time[order], flux[order], error[order]
            median = np.median(flux)
            flux, error = flux / median, error / median
            exposure = float(header['TIMEDEL'])
        sources.append(dict(sector=row['sector'], path=str(path), sha256=actual,
                            flux_column=col, error_column=errcol, exposure_days=exposure,
                            n_quality_zero=len(time), error_fallback=error_fallback))
        cycles = np.rint((time - EPOCH) / P).astype(int)
        for cycle in np.unique(cycles):
            x = time - EPOCH - cycle * P
            near = np.abs(x) < .20
            tx, fx, ex = x[near], flux[near], error[near]
            outside = np.abs(tx) > 1.20 * WIDTH
            inside = np.abs(tx) < WIDTH / 2 + exposure / 2
            if (len(tx) < 12 or inside.sum() < 2
                    or (outside & (tx < 0)).sum() < 3
                    or (outside & (tx > 0)).sum() < 3):
                continue
            # No polynomial spanning a telemetry/data-quality gap near eclipse.
            if np.max(np.diff(tx)) > max(3.05 * exposure, .065):
                continue
            baseline = np.column_stack([np.ones(len(tx)), tx, tx * tx])
            beta = np.linalg.lstsq(baseline[outside], fx[outside], rcond=None)[0]
            residual = fx[outside] - baseline[outside] @ beta
            floor = max(1.4826 * np.median(np.abs(residual - np.median(residual))),
                        np.median(ex))
            ex = np.maximum(ex, floor)
            events.append(dict(sector=row['sector'], cycle=int(cycle),
                               x=tx, f=fx, e=ex, exposure=exposure,
                               noise_floor=floor, baseline=baseline))
    return events, sources


def fit_event(event, shape):
    # Separate shapes for alternating eclipses; shared shift and recurrence P.
    shift, dp, w0, w1, r0, r1 = shape
    parity = event['cycle'] % 2
    width, ingress = ([w0, w1][parity], [r0, r1][parity])
    x = event['x'] - shift - event['cycle'] * dp
    template = integrated_trapezoid(x, width, ingress, event['exposure'])
    design = np.column_stack([event['baseline'], -template])
    a, y = design / event['e'][:, None], event['f'] / event['e']
    beta = np.linalg.lstsq(a, y, rcond=None)[0]
    residual = y - a @ beta
    covariance = np.linalg.pinv(a.T @ a) * max(1., np.sum(residual ** 2) / (len(y) - 4))
    return dict(depth=float(beta[-1]), error=float(np.sqrt(covariance[-1, -1])),
                residual=residual, baseline=design[:, :3] @ beta[:3],
                model=design @ beta, relative_time=x, template=template,
                chi2=float(np.sum(residual ** 2)))


def statistics(records, seed=71316387):
    by_sector = {}
    for sector in sorted({r['sector'] for r in records}):
        use = [r for r in records if r['sector'] == sector]
        parity, depth, error = (np.asarray([r[key] for r in use]) for key in ['parity', 'depth', 'error'])
        means = [float(np.mean(depth[parity == k])) for k in range(2)]
        sem = [max(float(np.std(depth[parity == k], ddof=1) / np.sqrt(np.sum(parity == k))),
                   float(np.sqrt(np.sum(error[parity == k] ** 2)) / np.sum(parity == k))) for k in range(2)]
        by_sector[str(sector)] = dict(n_events=len(use), parity_counts=[int(np.sum(parity == k)) for k in range(2)],
                                     mean_even_odd=means, error_even_odd=sem,
                                     odd_minus_even=means[1] - means[0],
                                     difference_error=float(np.hypot(*sem)),
                                     conditional_snr=float(abs(means[1] - means[0]) / np.hypot(*sem)))
    # Comparing neighboring cycles reduces sensitivity to slow dilution/trends.
    lookup = {(r['sector'], r['cycle']): r for r in records}
    pairs = []
    for r in records:
        if r['parity'] != 0:
            continue
        other = lookup.get((r['sector'], r['cycle'] + 1))
        if other is not None:
            pairs.append(dict(sector=r['sector'], even_cycle=r['cycle'],
                              odd_minus_even=other['depth'] - r['depth'],
                              error=float(np.hypot(r['error'], other['error']))))
    difference = np.array([r['odd_minus_even'] for r in pairs])
    error = np.array([r['error'] for r in pairs])
    mean = float(np.mean(difference))
    sem = max(float(np.std(difference, ddof=1) / np.sqrt(len(difference))),
              float(np.sqrt(np.sum(error ** 2)) / len(error)))
    weights = 1 / error ** 2
    weighted_mean = float(np.sum(weights * difference) / np.sum(weights))
    weighted_chi2 = float(np.sum(weights * (difference - weighted_mean) ** 2))
    weighted_diagonal_error = float(np.sqrt(max(1., weighted_chi2 / (len(pairs) - 1)) / np.sum(weights)))
    sectors = np.asarray([r['sector'] for r in pairs])
    scores = [float(np.sum(weights[sectors == s] * (difference[sectors == s] - weighted_mean)))
              for s in np.unique(sectors)]
    clustered_error = float(np.sqrt(np.sum(np.asarray(scores) ** 2) / np.sum(weights) ** 2
                                    * len(scores) / (len(scores) - 1)))
    max_error_index = int(np.argmax(error))
    without_noisiest = np.delete(difference, max_error_index)
    rng = np.random.default_rng(seed)
    samples = np.mean(rng.choice(difference, size=(20000, len(difference)), replace=True), axis=1)
    # Also use contiguous four-pair blocks, kept inside observing sectors.
    blocks = []
    for s in sorted({r['sector'] for r in pairs}):
        d = np.asarray([r['odd_minus_even'] for r in pairs if r['sector'] == s])
        blocks.extend([d[start:start + 4] for start in range(0, len(d), 4)])
    block_mean = []
    for _ in range(20000):
        chosen = [blocks[k] for k in rng.integers(0, len(blocks), len(blocks))]
        block_mean.append(np.mean(np.concatenate(chosen)))
    return dict(per_sector=by_sector, adjacent_pairs=pairs,
                weighted_paired_contrast=dict(odd_minus_even=weighted_mean,
                    residual_inflated_diagonal_error=weighted_diagonal_error,
                    diagonal_conditional_gaussian_snr=abs(weighted_mean) / weighted_diagonal_error,
                    sector_cluster_error=clustered_error,
                    conservative_cluster_conditional_snr=abs(weighted_mean) / max(clustered_error, weighted_diagonal_error),
                    n_sector_clusters=len(scores), chi2=weighted_chi2, dof=len(pairs) - 1,
                    note='Inverse-variance weights downweight one poorly constrained sector28 eclipse. Cluster uncertainty permits an independent mean offset in each sector; only five clusters means its tail calibration is unreliable. Error-model choices materially change the quoted sigma.'),
                sensitivity_omitting_one_largest_error_pair=dict(
                    excluded_pair=pairs[max_error_index],
                    mean_odd_minus_even=float(np.mean(without_noisiest)),
                    scatter_standard_error=float(np.std(without_noisiest, ddof=1) / np.sqrt(len(without_noisiest))),
                    note='Diagnostic only: retained in the principal analyses. Selection uses the measurement error rather than its eclipse depth.'),
                paired_contrast=dict(n_pairs=len(pairs), odd_minus_even=mean, error=sem,
                                     conditional_gaussian_snr=abs(mean) / sem,
                                     descriptive_student_t_pvalue=float(2 * student_t.sf(abs(mean) / sem, len(difference) - 1)),
                                     pair_bootstrap_95_percent=np.quantile(samples, [.025, .975]).tolist(),
                                     four_pair_block_bootstrap_95_percent=np.quantile(block_mean, [.025, .975]).tolist(),
                                     pair_bootstrap_standard_deviation=float(np.std(samples, ddof=1)),
                                     block_bootstrap_standard_deviation=float(np.std(block_mean, ddof=1)),
                                     limitations='Conditional fixed-template test. Shape/ephemeris were fit on these data; pairs may share systematics. Bootstrap confidence intervals do not establish Gaussian tail probabilities or survey false-alarm rates.'))


def analyze(column):
    events, sources = prepare(column)
    initial = [0., 0., WIDTH, WIDTH, .35, .35]
    def residual(params):
        return np.concatenate([fit_event(event, params)['residual'] for event in events])
    fit = least_squares(residual, initial,
                        bounds=([-.025, -3e-6, .025, .025, .02, .02],
                                [.025, 3e-6, .12, .12, .5, .5]),
                        x_scale=[.003, 1e-6, .06, .06, .3, .3], max_nfev=150)
    shape = fit.x
    records, plotdata = [], []
    for event in events:
        result = fit_event(event, shape)
        records.append(dict(sector=event['sector'], cycle=event['cycle'], parity=event['cycle'] % 2,
                            time_btjd=float(EPOCH + shape[0] + event['cycle'] * (P + shape[1])),
                            depth=result['depth'], error=result['error'], chi2=result['chi2'],
                            n_points=len(event['x']), empirical_noise_floor=event['noise_floor']))
        plotdata.append(dict(sector=event['sector'], parity=event['cycle'] % 2,
                             x=result['relative_time'], y=event['f'] - result['baseline'] + 1))
    return dict(sources=sources, shape_fit=dict(success=bool(fit.success),
                message=fit.message, epoch_btjd=float(EPOCH + shape[0]),
                recurrence_period_days=float(P + shape[1]),
                even_odd_durations_days=shape[2:4].tolist(),
                even_odd_ingress_fractions=shape[4:6].tolist(),
                chi2=float(np.sum(fit.fun ** 2))),
                event_depths=records, depth_statistics=statistics(records)), plotdata


def numerical_controls():
    """Check finite-exposure integrals independently by dense quadrature."""
    x = np.linspace(-.10, .10, 500)
    errors = {}
    for exposure in [1 / 48, 1 / 144, 1 / 432]:
        analytical = integrated_trapezoid(x, .061, .3, exposure)
        offset = np.linspace(-exposure / 2, exposure / 2, 4001)
        instantaneous = np.clip(np.minimum(
            (x[:, None] + offset + .061 / 2) / (.061 * .3),
            (.061 / 2 - x[:, None] - offset) / (.061 * .3)), 0, 1)
        numerical = np.trapezoid(instantaneous, offset, axis=1) / exposure
        errors[str(exposure)] = float(np.max(np.abs(analytical - numerical)))
        assert errors[str(exposure)] < 1e-7
    events, _ = prepare('SAP_FLUX')
    params = [.001, -1e-6, .061, .061, .3, .3]
    depths, parities = [], []
    for event in events:
        xx = event['x'] - params[0] - event['cycle'] * params[1]
        event['f'] = (1 + .001 * event['x'] + .004 * event['x'] ** 2
                      - .003 * integrated_trapezoid(xx, .061, .3, event['exposure']))
        depths.append(fit_event(event, params)['depth'])
        parities.append(event['cycle'] % 2)
    depths, parities = np.asarray(depths), np.asarray(parities)
    max_error = float(np.max(np.abs(depths - .003)))
    contrast = float(np.mean(depths[parities == 1]) - np.mean(depths[parities == 0]))
    assert max_error < 1e-10 and abs(contrast) < 1e-10
    return dict(integration_max_error_by_exposure=errors,
                noiseless_actual_cadence_equal_depth=dict(n_events=len(events),
                    true_depth=.003, max_recovered_depth_error=max_error,
                    odd_minus_even=contrast,
                    note='Numerical control of cadence integration and baseline/depth projection at a known fixed shape. Not an injection-recovery completeness measurement or noise false-alarm test.'))


def main():
    output = dict(tic=71316387, checked_utc=datetime.now(timezone.utc).isoformat(),
                  method='Independent raw-flux local quadratic baseline plus finite-exposure trapezoids. Separate even/odd shape; unconstrained per-event depths; empirical out-of-event error floor and fit-residual inflation.',
                  source_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  catalog_ephemeris=dict(period_days=P, epoch_btjd=EPOCH, duration_days=WIDTH),
                  caveats=['Amplitudes describe unconvolved trapezoid peaks, not physical eclipse depths or radius ratios.',
                           'Four-pair blocks only approximate within-sector noise correlations.',
                           'Alternating eclipses support a double-period binary but do not locate the eclipsing source or confirm stellar masses.',
                           'A planet at the doubled period could show primary/secondary eclipses, but a several-thousand-ppm optical occultation of a Sun-like host would require an exceptionally luminous object.'])
    output['numerical_controls'] = numerical_controls()
    plots = {}
    for col in ['SAP_FLUX', 'processed']:
        output[col], plots[col] = analyze(col)
        print(col, output[col]['shape_fit'], output[col]['depth_statistics']['paired_contrast'])
        print('weighted', output[col]['depth_statistics']['weighted_paired_contrast'])
        print('sectors', output[col]['depth_statistics']['per_sector'])
    (HERE / 'independent_audit.json').write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    fig, axes = plt.subplots(3, 2, figsize=(12, 11))
    colors = {1:'tab:blue', 28:'tab:orange', 68:'tab:green', 95:'tab:red', 105:'tab:purple'}
    for colindex, col in enumerate(['SAP_FLUX', 'processed']):
        result = output[col]
        for parity in [0, 1]:
            ax = axes[parity, colindex]
            for sector in colors:
                rows = [r for r in plots[col] if r['sector'] == sector and r['parity'] == parity]
                if rows:
                    ax.scatter(np.concatenate([r['x'] for r in rows]) * 24,
                               (1 - np.concatenate([r['y'] for r in rows])) * 1e6,
                               s=3, alpha=.25, color=colors[sector], label=f'Sector {sector}')
            ax.set_xlim(-2.5,2.5); ax.set_ylim(9000,-4000)
            ax.set_xlabel('Time from eclipse center (hours)'); ax.set_ylabel('Local depth (ppm)')
            ax.set_title(f'{col}, {"even/shallow" if parity == 0 else "odd/deep"} events')
        rows = result['event_depths']
        for parity, marker in [(0,'o'),(1,'s')]:
            rr = [r for r in rows if r['parity'] == parity]
            axes[2,colindex].errorbar([r['time_btjd'] for r in rr], [r['depth'] * 1e6 for r in rr],
                                    [r['error'] * 1e6 for r in rr], fmt=marker, ms=3, alpha=.7, label=f'Parity {parity}')
        axes[2,colindex].set_xlabel('BTJD (TDB)'); axes[2,colindex].set_ylabel('Fitted peak depth (ppm)')
        axes[2,colindex].legend()
    axes[0,0].legend(fontsize=8)
    fig.suptitle('TIC 71316387 independent finite-cadence event-depth audit')
    fig.tight_layout(); fig.savefig(HERE / 'independent_audit.png', dpi=160)


if __name__ == '__main__':
    main()
