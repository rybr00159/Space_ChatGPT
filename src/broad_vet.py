"""Independent raw-flux and aperture vetting for broad QLP survey leads.

Reads the ongoing survey SQLite database read-only. Never alters the frozen
survey method or raw data. Every signal remains a lead until source/pixel,
multi-sector, catalog/literature, and survey-wide null tests are complete.
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

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.optimize import least_squares

from archive import verify_fits_checksums
from fast_screen import FastScreenConfig, event_stats, phase, preprocess_curve, screen_curve

os.environ.setdefault('MPLCONFIGDIR', '/workspace/tess-research-data/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATABASE = Path(os.environ.get('TESS_DATA_DIR', '/workspace/tess-research-data')) / 'broad/qlp_s0001/survey.sqlite'
OUT = ROOT / 'candidates/broad_s1'
METHOD_VERSION = 'raw-aperture-all-orbital-phase-raw-variability-v3'


def method_identity():
    return dict(version=METHOD_VERSION, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                frozen_screen_sha256=hashlib.sha256((ROOT / 'src/fast_screen.py').read_bytes()).hexdigest(),
                mandatory_tests=['raw_local_quadratic', 'three_median_windows', 'four_apertures',
                                 'odd_even', 'half_phase_secondary', 'all_phase_secondary', 'qlp_centroids',
                                 'robust_raw_orbital_and_half_period_variability'])


def atomic_json(path, payload):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(payload, indent=2, allow_nan=False))
    temp.replace(path)


def mad(x):
    return float(1.4826 * np.median(abs(x - np.median(x)))) if len(x) else 0.


def read_source(record):
    path = Path(record['path'])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != record['sha256']:
        raise ValueError('Survey source hash mismatch')
    with fits.open(path, memmap=False) as hd:
        checksums = verify_fits_checksums(path, hd)
        h, d = hd[1].header, hd[1].data
        if h.get('TIMESYS') != 'TDB' or h.get('BJDREFI', 0) + h.get('BJDREFF', 0) != 2457000:
            raise ValueError('Unrecognized time coordinate')
        arrays = {k: np.asarray(d[k], dtype=float) for k in
                  ('TIME', 'QUALITY', 'KSPSAP_FLUX', 'KSPSAP_FLUX_ERR', 'SAP_FLUX',
                   'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG', 'SAP_X', 'SAP_Y', 'SAP_BKG')}
    return arrays, checksums


def local_quadratic(t, flux, error, period, epoch, duration, degree=2):
    """Fit local out-of-event polynomials independently around every recurrence.

    Baselines use both sides outside 1.15 transit widths and do not bridge >2.5
    cadence gaps. Parameter errors have an out-of-event residual MAD floor.
    Polynomial extrapolation error is included in cadence errors approximately;
    event scatter provides an additional empirical mean-depth error.
    """
    good = np.isfinite(t) & np.isfinite(flux) & np.isfinite(error) & (error > 0) & (flux > 0)
    t, flux, error = t[good], flux[good], error[good]
    order = np.argsort(t)
    t, flux, error = t[order], flux[order], error[order]
    flux, error = flux / np.median(flux), error / np.median(flux)
    cadence = float(np.median(np.diff(t)))
    cycles = np.round((t - epoch) / period).astype(int)
    residuals, times, errors, event_ids = [], [], [], []
    events = []
    for cycle in np.unique(cycles):
        center = epoch + cycle * period
        half_window = min(.35 * period, max(3.5 * duration, .25))
        use = abs(t - center) < half_window
        tx, fx, ex = t[use] - center, flux[use], error[use]
        if len(tx) < 12:
            continue
        inside = abs(tx) < .5 * duration
        outside = abs(tx) > 1.15 * duration
        if inside.sum() < 2 or (outside & (tx < 0)).sum() < 3 or (outside & (tx > 0)).sum() < 3:
            continue
        # Require a contiguous baseline/event island; do not fit through a
        # downlink or rejected-quality gap across the predicted transit.
        if np.max(np.diff(tx)) > max(.06, 3 * cadence):
            continue
        design = np.column_stack([tx ** k for k in range(degree + 1)])
        use_oot = outside.copy()
        for _ in range(2):
            a = design[use_oot] / ex[use_oot, None]
            b = fx[use_oot] / ex[use_oot]
            beta, _, _, _ = np.linalg.lstsq(a, b, rcond=None)
            trend = design @ beta
            scatter = mad(fx[use_oot] - trend[use_oot])
            use_oot = outside & (abs(fx - trend) < max(5 * scatter, 5 * np.median(ex)))
        if use_oot.sum() < degree + 5:
            continue
        normalized = fx / trend
        floor = max(mad(normalized[use_oot]), float(np.median(ex / trend)))
        en = np.maximum(ex / trend, floor)
        covariance = np.linalg.pinv(a.T @ a) * max(1., np.sum(((fx[use_oot] - trend[use_oot]) / ex[use_oot]) ** 2) / max(1, use_oot.sum() - degree - 1))
        baseline_var = np.einsum('ij,jk,ik->i', design, covariance, design) / trend ** 2
        en = np.sqrt(en ** 2 + np.maximum(baseline_var, 0))
        wi = 1 / en[inside] ** 2
        depth = float(np.sum(wi * (1 - normalized[inside])) / wi.sum())
        err = float(1 / np.sqrt(wi.sum()))
        events.append(dict(cycle=int(cycle), time=float(center), n=int(inside.sum()), depth=depth,
                           depth_err=err, snr=depth / err, n_baseline=int(use_oot.sum())))
        residuals.append(normalized)
        times.append(tx + center)
        errors.append(en)
        event_ids.append(np.full(len(tx), cycle))
    if not events:
        return dict(n_events=0, depth=None, depth_err=None, snr=0., events=[], odd_even_sigma=None), None
    d, e = np.array([v['depth'] for v in events]), np.array([v['depth_err'] for v in events])
    w = 1 / e ** 2
    depth = float(np.sum(d * w) / w.sum())
    diagerr = float(1 / np.sqrt(w.sum()))
    # This floor prevents dozens of inconsistent events from manufacturing a
    # precision unsupported by their observed dispersion.
    scattererr = float(np.std(d, ddof=1) / np.sqrt(len(d))) if len(d) > 1 else diagerr
    err = max(diagerr, scattererr)
    parity = np.array([v['cycle'] % 2 for v in events])
    odd_even = None
    if (parity == 0).sum() >= 2 and (parity == 1).sum() >= 2:
        means = [np.sum(w[parity == k] * d[parity == k]) / w[parity == k].sum() for k in (0, 1)]
        pe = [max(1 / np.sqrt(w[parity == k].sum()), np.std(d[parity == k], ddof=1) / np.sqrt((parity == k).sum())) for k in (0, 1)]
        odd_even = float(abs(means[0] - means[1]) / np.hypot(*pe))
    stats = dict(n_events=len(events), n_positive_events=int((d > 0).sum()), depth=depth, depth_err=err,
                 diagonal_depth_err=diagerr, event_scatter_depth_err=scattererr, snr=depth / err,
                 events=events, odd_even_sigma=odd_even)
    return stats, dict(t=np.concatenate(times), f=np.concatenate(residuals), e=np.concatenate(errors), cycle=np.concatenate(event_ids))


def integrated_trapezoid(x, duration, ingress, exposure):
    """Exact box-exposure integral of a unit-height symmetric trapezoid."""
    a, r = duration / 2, max(ingress, 1e-7)
    def primitive(z):
        positive = lambda v: np.maximum(v, 0.) ** 2
        return (positive(z + a) - positive(z + a - r) - positive(z - a + r) + positive(z - a)) / (2 * r)
    return np.clip((primitive(x + exposure / 2) - primitive(x - exposure / 2)) / exposure, 0, 1)


def fit_shape(data, period, epoch, duration, exposure):
    if data is None or len(data['t']) < 20:
        return None
    x, y, e = phase(data['t'], period, epoch), data['f'], data['e']
    use = abs(x) < min(.35 * period, max(2 * duration, .25))
    x, y, e = x[use], y[use], e[use]
    if np.ndim(exposure):
        exposure = np.asarray(exposure)[use]
    depth = max(.0001, float(1 - np.mean(y[abs(x) < duration / 2])))
    def residual(par):
        shift, dur, ratio, dep, offset = par
        return (y - (1 + offset - dep * integrated_trapezoid(x - shift, dur, dur * ratio, exposure))) / e
    upper_d = min(.3 * period, max(duration * 2.5, .15))
    initial = [0., duration, .2, min(depth, .19), 0.]
    fit = least_squares(residual, initial,
                        bounds=([-.7 * duration, max(float(np.min(exposure)) * .6, .01), .02, 0, -.025],
                                [.7 * duration, upper_d, .5, .2, .025]))
    shift, dur, ratio, dep, offset = fit.x
    chi = float(np.sum(fit.fun ** 2))
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * max(1., chi / max(1, len(x) - len(fit.x)))
    uncertainty = np.sqrt(np.maximum(np.diag(covariance), 0))
    flat_chi = float(np.sum(((y - np.average(y, weights=1 / e ** 2)) / e) ** 2))
    # Wide ingress/high duration or a depth driven by the bounds must be treated
    # as shape warnings, never automatic stellar/planet classifications.
    return dict(epoch_btjd=float(epoch + shift), duration_days=float(dur), ingress_fraction=float(ratio),
                depth=float(dep), baseline_offset=float(offset), depth_err=float(uncertainty[3]),
                duration_err=float(uncertainty[1]), epoch_err=float(uncertainty[0]),
                chi2=chi, n_points=len(x), delta_chi2_flat=flat_chi - chi,
                warning='Fixed period; diagonal/covariance errors approximate; no limb darkening or source validation.')


def variability_diagnostics(t, f, e, period, epoch, duration, additional_epoch=None, additional_duration=None):
    """Robust weighted orbital harmonics outside generous eclipse windows.

    Reports five coefficients for c0+s1*sin(theta)+c1*cos(theta)+s2*sin(2theta)
    +c2*cos(2theta), theta=2*pi*(BTJD-epoch)/period. Covariance errors have an
    observed residual floor and reduced-chi-square inflation. These remain
    conditional diagonal estimates, not global false-alarm probabilities.
    """
    good = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0) & (f > 0)
    t, f, e = t[good], f[good], e[good]
    f, e = f / np.median(f), e / np.median(f)
    oot = (abs(phase(t, period, epoch)) > duration) & (abs(phase(t, period, epoch + period / 2)) > duration)
    if additional_epoch is not None:
        oot &= abs(phase(t, period, additional_epoch)) > (additional_duration or duration)
    x = (t - epoch) / period * 2 * np.pi
    design = np.column_stack([np.ones(len(t)), np.sin(x), np.cos(x), np.sin(2 * x), np.cos(2 * x)])
    if oot.sum() < 20:
        return dict(status='insufficient_out_of_eclipse_data', n_fit=int(oot.sum()), raw_mad=mad(f), residual_mad=None,
                    orbital_semiamplitude=None, half_period_semiamplitude=None,
                    orbital_semiamplitude_err=None, half_period_semiamplitude_err=None,
                    orbital_snr=None, half_period_snr=None)
    use = oot.copy()
    ee = e.copy()
    for _ in range(3):
        a, b = design[use] / ee[use, None], f[use] / ee[use]
        beta = np.linalg.lstsq(a, b, rcond=None)[0]
        residual = f - design @ beta
        scatter = mad(residual[use])
        ee = np.maximum(e, max(scatter, float(np.median(e[use]))))
        next_use = oot & (abs(residual - np.median(residual[use])) < 5 * ee)
        if next_use.sum() < 20:
            break
        use = next_use
    # Refit with the final robust selection and final conservative error floor.
    a, b = design[use] / ee[use, None], f[use] / ee[use]
    beta = np.linalg.lstsq(a, b, rcond=None)[0]
    residual = f - design @ beta
    chi2 = float(np.sum((residual[use] / ee[use]) ** 2))
    covariance = np.linalg.pinv(a.T @ a) * max(1., chi2 / max(1, use.sum() - 5))
    amplitudes, errors = [], []
    for indices in ((1, 2), (3, 4)):
        amplitude = float(np.hypot(*beta[list(indices)]))
        # Largest covariance-axis error is conservative under an unknown phase.
        error = float(np.sqrt(max(0., np.linalg.eigvalsh(covariance[np.ix_(indices, indices)]).max())))
        amplitudes.append(amplitude)
        errors.append(error)
    return dict(status='measured', raw_mad=mad(f), residual_mad=mad(residual[use]),
                n_out_of_eclipse=int(oot.sum()), n_fit=int(use.sum()), n_robust_rejected=int(oot.sum() - use.sum()),
                orbital_semiamplitude=amplitudes[0], half_period_semiamplitude=amplitudes[1],
                orbital_semiamplitude_err=errors[0], half_period_semiamplitude_err=errors[1],
                orbital_snr=amplitudes[0] / errors[0] if errors[0] else None,
                half_period_snr=amplitudes[1] / errors[1] if errors[1] else None,
                coefficient_order=['intercept', 'sin_orbital', 'cos_orbital', 'sin_half_period', 'cos_half_period'],
                coefficients=[float(v) for v in beta], coefficient_covariance=covariance.tolist(),
                residual_chi2=chi2, residual_dof=int(use.sum() - 5),
                eclipse_windows=dict(primary_half_width_days=duration, half_phase_half_width_days=duration,
                                     additional_epoch_btjd=additional_epoch, additional_half_width_days=additional_duration),
                model='c0+s1*sin(theta)+c1*cos(theta)+s2*sin(2*theta)+c2*cos(2*theta);theta=2*pi*(BTJD-epoch)/P',
                limitations='Conditional diagonal errors; robust scatter/clipping does not calibrate correlated noise or the search trials.')


def refine_ephemeris(t, f, e, peak):
    """Refine within three conservative coarse-grid cells, preserving provenance."""
    data = preprocess_curve(t, f, e, config=FastScreenConfig(width_days=2.))
    baseline = float(np.ptp(data['t']))
    cadence = float(np.median(np.diff(data['t'])))
    period = peak['period']
    shortest = max(.05, 1.5 * cadence)
    half_range = max(.001 * period, 3 * shortest / baseline * period)
    periods = np.linspace(max(.3, period - half_range), min(13., period + half_range), 301)
    durations = np.unique(np.r_[.05, .075, .10, .125, .16, .20, peak['duration_days']])
    durations = durations[(durations >= 1.5 * cadence) & (durations < .25 * periods.min())]
    if not len(durations):
        durations = np.array([peak['duration_days']])
    result = BoxLeastSquares(data['t'], data['f'], data['e']).power(periods, durations, oversample=15, objective='snr')
    j = int(np.argmax(result.power))
    return dict(period=float(result.period[j]), t0_btjd=float(result.transit_time[j]),
                duration_days=float(result.duration[j]), depth=float(result.depth[j]), snr=float(result.depth_snr[j]),
                period_interval_days=[float(periods[0]), float(periods[-1])], period_grid_step=float(periods[1] - periods[0]),
                note='Local BLS refinement, not period uncertainty; requires independent sector timing.')


def centroid_diagnostic(t, coordinate, period, epoch, duration):
    """QLP aperture centroids provide a warning test, not pixel localization."""
    good = np.isfinite(t) & np.isfinite(coordinate)
    t, coordinate = t[good], coordinate[good]
    event = []
    for k in np.unique(np.round((t - epoch) / period).astype(int)):
        center = epoch + k * period
        use = abs(t - center) < min(.35 * period, max(3.5 * duration, .25))
        x, y = t[use] - center, coordinate[use]
        inside, outside = abs(x) < duration / 2, abs(x) > duration * 1.15
        if inside.sum() < 2 or (outside & (x < 0)).sum() < 3 or (outside & (x > 0)).sum() < 3:
            continue
        a = np.column_stack([np.ones(len(x)), x, x ** 2])
        beta = np.linalg.lstsq(a[outside], y[outside], rcond=None)[0]
        residual = y - a @ beta
        scatter = mad(residual[outside])
        if scatter <= 0:
            continue
        event.append(dict(cycle=int(k), offset_pixels=float(np.mean(residual[inside])),
                          error_pixels=float(scatter / np.sqrt(inside.sum()))))
    if not event:
        return dict(n_events=0, offset_pixels=None, error_pixels=None, sigma=None)
    d, e = np.array([v['offset_pixels'] for v in event]), np.array([v['error_pixels'] for v in event])
    w = 1 / e ** 2
    mean = float(np.sum(w * d) / w.sum())
    error = max(float(1 / np.sqrt(w.sum())), float(np.std(d, ddof=1) / np.sqrt(len(d))) if len(d) > 1 else 0.)
    return dict(n_events=len(event), offset_pixels=mean, error_pixels=error, sigma=abs(mean) / error,
                note='QLP centroid local-baseline warning only; no pixel source attribution.')


def fit_event_timings(local_data, sector_ids, period, epoch, shape):
    """Conditional per-event shape times and linear ephemeris with scatter floor."""
    if shape is None:
        return None
    duration, ingress = shape['duration_days'], shape['duration_days'] * shape['ingress_fraction']
    events = []
    for data, sector in zip(local_data, sector_ids):
        for cycle in np.unique(data['cycle']):
            use = data['cycle'] == cycle
            center = epoch + cycle * period
            x, y, e, exposure = (data[k][use] for k in ('t', 'f', 'e', 'exposure'))
            x = x - center
            if len(x) < 8 or (abs(x) < duration / 2).sum() < 2:
                continue
            def residual(par):
                shift, depth, offset = par
                return (y - (1 + offset - depth * integrated_trapezoid(x - shift, duration, ingress, exposure))) / e
            fit = least_squares(residual, [0., shape['depth'], 0.],
                                bounds=([-.7 * duration, 0., -.025], [.7 * duration, .2, .025]))
            chi = float(np.sum(fit.fun ** 2))
            covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * max(1., chi / max(1, len(x) - 3))
            error = np.sqrt(np.maximum(np.diag(covariance), 0))
            delta = float(np.sum(((y - np.average(y, weights=1 / e ** 2)) / e) ** 2) - chi)
            if delta < 16 or abs(fit.x[0]) > duration * .65 or error[0] <= 0:
                continue
            events.append(dict(sector=int(sector), cycle=int(cycle), time_btjd=float(center + fit.x[0]),
                               time_error_days=float(error[0]), depth=float(fit.x[1]), depth_err=float(error[1]), delta_chi2=delta))
    if len(events) < 4:
        return dict(n_events=len(events), events=events, fitted_ephemeris=None)
    cycles = np.array([v['cycle'] for v in events], float)
    times = np.array([v['time_btjd'] for v in events])
    errors = np.array([v['time_error_days'] for v in events])
    refcycle = float(np.mean(cycles))
    design = np.column_stack([np.ones(len(times)), cycles - refcycle])
    a, b = design / errors[:, None], times / errors
    beta = np.linalg.lstsq(a, b, rcond=None)[0]
    residuals = times - design @ beta
    chi = float(np.sum((residuals / errors) ** 2))
    covariance = np.linalg.pinv(a.T @ a) * max(1., chi / (len(times) - 2))
    uncertainty = np.sqrt(np.maximum(np.diag(covariance), 0))
    epoch_fit = float(beta[0] - beta[1] * refcycle)
    epoch_error = float(np.sqrt(np.array([1., -refcycle]) @ covariance @ np.array([1., -refcycle])))
    leave_sector_out = []
    sectors = np.array([v['sector'] for v in events])
    for sector in np.unique(sectors):
        use = sectors != sector
        bp = np.linalg.lstsq(a[use], b[use], rcond=None)[0]
        leave_sector_out.append(dict(omitted_sector=int(sector), period=float(bp[1])))
    return dict(n_events=len(events), events=events,
                fitted_ephemeris=dict(period=float(beta[1]), period_error_days=float(uncertainty[1]),
                                      epoch_btjd=epoch_fit, epoch_error_days=epoch_error,
                                      timing_chi2=chi, timing_dof=len(times) - 2,
                                      median_time_error_minutes=float(np.median(errors) * 1440),
                                      note='Conditional fixed shape; diagonal errors inflated by event and timing scatter, not a physical confirmation.'),
                leave_one_sector_out=leave_sector_out)


def pool_sector_depths(sector_tests, column, event='secondary', min_events=2):
    valid = [(v['sector'], v['column_tests'][column][event]) for v in sector_tests
             if v['column_tests'][column][event]['n_events'] >= min_events]
    if not valid:
        return None
    depths = np.array([v['depth'] for _, v in valid])
    errors = np.array([v['depth_err'] for _, v in valid])
    w = 1 / errors ** 2
    mean = float(np.sum(w * depths) / w.sum())
    chi2 = float(np.sum(((depths - mean) / errors) ** 2))
    err = float(1 / np.sqrt(w.sum())) * np.sqrt(max(1., chi2 / max(1, len(depths) - 1)))
    if len(depths) > 1:
        err = max(err, float(np.std(depths, ddof=1) / np.sqrt(len(depths))))
    return dict(depth=mean, depth_err=err, snr=mean / err, consistency_chi2=chi2,
                consistency_dof=len(depths) - 1, n_sectors=len(depths), sectors=[v for v, _ in valid],
                note='Targeted fixed-ephemeris measurement; conditional local-baseline errors with observed scatter floors.')


def vet(record, peak, label='unmasked_broad'):
    tic = record['tic']
    arrays, integrity = read_source(record)
    good = arrays['QUALITY'] == 0
    t, e = arrays['TIME'][good], arrays['KSPSAP_FLUX_ERR'][good]
    refined = refine_ephemeris(t, arrays['KSPSAP_FLUX'][good], e, peak)
    p, epoch, duration = refined['period'], refined['t0_btjd'], refined['duration_days']
    cadence = float(np.median(np.diff(np.sort(t[np.isfinite(t)]))))
    columns = ['KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG']
    tests, locals_, intermediates = {}, {}, {}
    for column in columns:
        f = arrays[column][good]
        loc, data = local_quadratic(t, f, e, p, epoch, duration)
        sec, _ = local_quadratic(t, f, e, p, epoch + p / 2, duration)
        locals_[column] = data
        tests[column] = dict(local_quadratic=loc, secondary_half_phase=sec,
                             variability=variability_diagnostics(t, f, e, p, epoch, duration))
        tests[column]['alternative_widths'] = []
        for width in (.5, 1., 2.):
            cfg = FastScreenConfig(width_days=width)
            flat = preprocess_curve(t, f, e, config=cfg)
            stats = event_stats(flat['t'], flat['f'], flat['e'], p, epoch, duration)
            tests[column]['alternative_widths'].append(dict(width_days=width, stats=stats))
            intermediates[f'{column}_w{width:g}_t'] = flat['t']
            intermediates[f'{column}_w{width:g}_f'] = flat['f']
            intermediates[f'{column}_w{width:g}_e'] = flat['e']
    shape = fit_shape(locals_['KSPSAP_FLUX'], p, epoch, duration, cadence)
    centroids = {axis: centroid_diagnostic(t, arrays[axis][good], p, epoch, duration) for axis in ('SAP_X', 'SAP_Y')}
    ft, ff, fe = (intermediates[f'KSPSAP_FLUX_w2_{k}'] for k in ('t', 'f', 'e'))
    primary_width = max(duration, shape['duration_days'] if shape else duration)
    outside_primary = abs(phase(ft, p, epoch)) > 1.25 * primary_width
    secondary_durations = np.array([.04, .06, .08, .10, .12, .16])
    secondary_durations = secondary_durations[(secondary_durations >= 1.5 * cadence) & (secondary_durations < .25 * p)]
    if not len(secondary_durations):
        secondary_durations = np.array([duration])
    additional = BoxLeastSquares(ft[outside_primary], ff[outside_primary], fe[outside_primary]).power([p], secondary_durations,
                                                                                                  oversample=15, objective='snr')
    additional_epoch, additional_duration = float(additional.transit_time[0]), float(additional.duration[0])
    additional_phase = float((additional_epoch - epoch) % p / p)
    for column in columns:
        stats, _ = local_quadratic(t, arrays[column][good], e, p, additional_epoch, additional_duration)
        tests[column]['additional_phase_secondary'] = stats
        tests[column]['variability'] = variability_diagnostics(t, arrays[column][good], e, p, epoch, primary_width,
                                                               additional_epoch=additional_epoch,
                                                               additional_duration=additional_duration)
    additional_report = dict(phase=additional_phase, epoch_btjd=additional_epoch, duration_days=additional_duration,
                             diagonal_snr=float(additional.depth_snr[0]),
                             note='Additional orbital phase selected after masking primary; score requires trials-aware confirmation.')
    primary = tests['KSPSAP_FLUX']['local_quadratic']
    sap = tests['SAP_FLUX']['local_quadratic']
    secondary = tests['KSPSAP_FLUX']['secondary_half_phase']
    flags = []
    if primary['snr'] < 7 or sap['snr'] < 7:
        flags.append('weak_independent_raw_local_detection')
    if primary['n_events'] < 3 or primary.get('n_positive_events', 0) < 3:
        flags.append('fewer_than_three_independent_positive_events')
    if secondary['snr'] > 5 and secondary.get('depth', 0) > .15 * (primary.get('depth') or 1):
        flags.append('significant_half_phase_secondary')
    additional_local = tests['KSPSAP_FLUX']['additional_phase_secondary']
    additional_raw = tests['SAP_FLUX']['additional_phase_secondary']
    if (additional_local['snr'] > 7 and additional_raw['snr'] > 7 and additional_local['n_events'] >= 2
            and (additional_local.get('depth') or 0) > .15 * (primary.get('depth') or 1)):
        flags.append('significant_additional_phase_secondary')
    if (primary.get('odd_even_sigma') or 0) > 4:
        flags.append('odd_even_depth_difference')
    if shape and shape['ingress_fraction'] > .45:
        flags.append('v_shaped_or_unresolved_ingress')
    if shape and shape['duration_days'] / p > .15:
        flags.append('fitted_duty_cycle_exceeds_point15')
    if any(v['sigma'] is not None and v['sigma'] > 5 and abs(v['offset_pixels']) > .001 for v in centroids.values()):
        flags.append('qlp_centroid_shift_warning')
    variability = tests['KSPSAP_FLUX']['variability']
    if primary.get('depth') and max(variability['orbital_semiamplitude'] or 0., variability['half_period_semiamplitude'] or 0.) > primary['depth'] * .5:
        flags.append('large_orbital_phase_variability')
    raw_variability = tests['SAP_FLUX']['variability']
    raw_variability_guard = dict(orbital_amplitude_to_primary_depth=None, half_period_amplitude_to_primary_depth=None,
                                amplitude_fraction_threshold=.5, conditional_snr_threshold=5.,
                                interpretation='A warning for stellar variability/ellipsoidal light; no automatic confirmation or false-positive proof.')
    if primary.get('depth') and primary['depth'] > 0:
        raw_variability_guard['orbital_amplitude_to_primary_depth'] = (raw_variability['orbital_semiamplitude'] or 0.) / primary['depth']
        raw_variability_guard['half_period_amplitude_to_primary_depth'] = (raw_variability['half_period_semiamplitude'] or 0.) / primary['depth']
        for kind, snr_key, ratio_key in (('orbital', 'orbital_snr', 'orbital_amplitude_to_primary_depth'),
                                       ('half_period', 'half_period_snr', 'half_period_amplitude_to_primary_depth')):
            if raw_variability_guard[ratio_key] > .5 and (raw_variability.get(snr_key) or 0.) > 5:
                flags.append(f'large_significant_raw_{kind}_variability')
    if primary.get('depth'):
        sml, lag = tests['KSPSAP_FLUX_SML']['local_quadratic'], tests['KSPSAP_FLUX_LAG']['local_quadratic']
        if sml.get('depth') and lag.get('depth'):
            sigma = abs(sml['depth'] - lag['depth']) / np.hypot(sml['depth_err'], lag['depth_err'])
            if sigma > 4 and abs(sml['depth'] - lag['depth']) > primary['depth'] * .3:
                flags.append('aperture_dependent_depth')
    result = dict(tic=tic, sector=1, checked_utc=datetime.now(timezone.utc).isoformat(), survey_label=label,
                  vetting_method=method_identity(),
                  input_peak=peak, refined_ephemeris=refined, catalog_flags=record['catalog_flags'], metadata=record['metadata'],
                  source={k: record[k] for k in ('url', 'path', 'sha256', 'raw_bytes')}, fits_checksums=integrity,
                  cadence_days=cadence, column_tests=tests, integrated_trapezoid=shape, centroid_tests=centroids,
                  additional_phase_search=additional_report,
                  raw_variability_guard=raw_variability_guard,
                  warnings=flags, assessment='unexplained_lead_requires_independent_source_and_novelty_validation' if not flags else 'fails_or_incomplete_initial_vetting',
                  limitations=['No source localization/pixel dilution model.', 'No survey-wide empirical false-alarm calibration.',
                               'Fixed-period trapezoid is not a limb-darkened transit model.', 'Catalog absence does not establish novelty.',
                               'Raw SAP/aperture errors use KSPSAP relative errors with local residual floors.'])
    directory = OUT / f'tic{tic}_p{peak["period"]:.6f}'
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 'vetting.json').write_text(json.dumps(result, indent=2, allow_nan=False))
    np.savez_compressed(directory / 'diagnostic_arrays.npz', **arrays, **intermediates)
    fig, axes = plt.subplots(4, 2, figsize=(12, 13), layout='constrained')
    for row, column in enumerate(columns):
        f = arrays[column][good]
        axes[row, 0].plot(t, f, '.', ms=2, color='.5')
        centers = np.arange(int(np.floor((np.nanmin(t) - epoch) / p)) - 1, int(np.ceil((np.nanmax(t) - epoch) / p)) + 1) * p + epoch
        for center in centers:
            axes[row, 0].axvline(center, color='tab:red', alpha=.3, lw=.6)
        axes[row, 0].set(xlabel='BTJD', ylabel=column)
        flat = preprocess_curve(t, f, e, config=FastScreenConfig(width_days=2.))
        ph = phase(flat['t'], p, epoch) / p
        axes[row, 1].plot(ph, flat['f'], '.', ms=2, color='.6')
        bins = np.linspace(-.5, .5, 101)
        idx = np.digitize(ph, bins)
        xb, yb = [], []
        for k in range(1, len(bins)):
            ii = idx == k
            if ii.sum():
                xb.append(float(np.mean(ph[ii])))
                yb.append(float(np.mean(flat['f'][ii])))
        axes[row, 1].plot(xb, yb, 'o-', ms=3, color='tab:blue')
        axes[row, 1].set(xlabel=f'Phase at P={p:.6f} d', ylabel='2 d median normalized flux')
    fig.suptitle(f'TIC {tic}; {"; ".join(flags) if flags else "Unexplained lead; validation incomplete"}', fontsize=10)
    fig.savefig(directory / 'diagnostics.png', dpi=130)
    plt.close(fig)
    lines = [f'# TIC {tic}, sector 1 lead', '', f'Refined period {p:.6f} d, epoch BTJD {epoch:.6f}, refined box duration {duration * 24:.2f} h.', '',
             f"Independent local primary: S/N {primary['snr']:.2f}, {primary['n_events']} sampled events; "
             f"raw SAP local S/N {sap['snr']:.2f}; half-phase secondary S/N {secondary['snr']:.2f}.", '',
             f"Assessment: **{result['assessment']}**.", '', 'Warnings: ' + (', '.join(flags) if flags else 'None from these limited tests.'), '',
             'See `vetting.json` for event depths, approximate uncertainty estimates, alternative detrending/apertures, and provenance. '
             '`diagnostic_arrays.npz` preserves small independent scientific arrays; raw archive FITS is cached outside Git.', '',
             'No discovery or confirmation claim. Pixel source attribution, additional sectors, comprehensive literature checks, '
             'and empirical correlated-noise/trials tests remain required.']
    (directory / 'README.md').write_text('\n'.join(lines) + '\n')
    print(tic, f'P={p:.6f}', 'local', round(primary['snr'], 2), 'SAP', round(sap['snr'], 2),
          'secondary', round(secondary['snr'], 2), 'flags', flags, flush=True)
    return result


def load_records(tics=None):
    db = sqlite3.connect(f'file:{DATABASE}?mode=ro', uri=True)
    if tics:
        rows = db.execute('SELECT result_json FROM stars WHERE tic IN (' + ','.join('?' for _ in tics) + ')', tics)
    else:
        rows = db.execute("SELECT result_json FROM stars WHERE status='screened'")
    records = [json.loads(raw) for raw, in rows if raw]
    db.close()
    return records


def vet_multisector(tic, selected_period=None, evidence_id=None):
    """Independent sector epochs, then long-baseline local-period BLS and depths."""
    base = load_records([tic])[0]
    if evidence_id:
        index = evidence_id.split('/')
        if int(index[0].split(':')[1]) != tic:
            raise ValueError('Evidence identifier target differs from selected TIC')
        scan_index, peak_index = int(index[2]), int(index[4])
        peak = base['scans'][scan_index]['peaks'][peak_index]
    elif selected_period is not None:
        peak = min((v for scan in base['scans'] for v in scan['peaks']), key=lambda v: abs(v['period'] - selected_period))
    else:
        peak = base['scans'][0]['peaks'][0]
    manifest = ROOT / f'data/broad_followup_queries/tic{tic}_qlp.json'
    rows = [dict(sector=1, path=base['path'], sha256=base['sha256'], url=base['url'], metadata=base['metadata'])]
    rows += [r for r in json.loads(manifest.read_text()) if r['status'] == 'retained']
    curves, independent, sources = [], [], []
    rawsave = {}
    for row in rows:
        d, checks = read_source(row)
        good = d['QUALITY'] == 0
        t, f, e = d['TIME'][good], d['KSPSAP_FLUX'][good], d['KSPSAP_FLUX_ERR'][good]
        flat = preprocess_curve(t, f, e, config=FastScreenConfig(width_days=2.))
        refined = refine_ephemeris(t, f, e, peak)
        primary, _ = local_quadratic(t, f, e, refined['period'], refined['t0_btjd'], refined['duration_days'])
        secondary, _ = local_quadratic(t, f, e, refined['period'], refined['t0_btjd'] + refined['period'] / 2, refined['duration_days'])
        independent.append(dict(sector=row['sector'], ephemeris=refined, primary=primary, secondary=secondary,
                                n_quality_good=int(good.sum()), cadence_days=float(np.median(np.diff(t[np.isfinite(t)])))))
        curves.append((row, d, good, flat))
        sources.append(dict(sector=row['sector'], url=row['url'], sha256=row['sha256'], path=row['path'], checks=checks))
        for col, array in d.items():
            rawsave[f's{row["sector"]}_{col}'] = array
    t, f, e = (np.concatenate([v[3][k] for v in curves]) for k in ('t', 'f', 'e'))
    order = np.argsort(t)
    t, f, e = t[order], f[order], e[order]
    # This is a targeted local period search after sector-1 detection. The joint
    # score is consequently not blind-search significance or a false-alarm rate.
    periods = np.linspace(peak['period'] - .02, peak['period'] + .02, 4001)
    bls = BoxLeastSquares(t, f, e)
    duration_grid = (.04, .06, .08, .10, .12, .16)
    power = bls.power(periods, duration_grid, oversample=15, objective='snr')
    j = int(np.argmax(power.power))
    p, epoch, duration = float(power.period[j]), float(power.transit_time[j]), float(power.duration[j])
    full, fullsec = [], []
    local_all, local_sectors = [], []
    for row, d, good, flat in curves:
        st, se = d['TIME'][good], d['KSPSAP_FLUX_ERR'][good]
        columns = {}
        for col in ('KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG'):
            primary, data = local_quadratic(st, d[col][good], se, p, epoch, duration)
            secondary, _ = local_quadratic(st, d[col][good], se, p, epoch + p / 2, duration)
            columns[col] = dict(primary=primary, secondary=secondary)
            if col == 'KSPSAP_FLUX' and data is not None:
                local_all.append(data)
                local_sectors.append(row['sector'])
                data['exposure'] = np.full(len(data['t']), float(np.median(np.diff(st[np.isfinite(st)]))))
                for k in data:
                    rawsave[f's{row["sector"]}_local_{k}'] = data[k]
        full.append(dict(sector=row['sector'], metadata=row['metadata'], column_tests=columns))
    joined = dict(t=np.concatenate([v['t'] for v in local_all]), f=np.concatenate([v['f'] for v in local_all]),
                  e=np.concatenate([v['e'] for v in local_all]), exposure=np.concatenate([v['exposure'] for v in local_all])) if local_all else None
    shape = fit_shape(joined, p, epoch, duration, joined['exposure']) if joined else None
    timings = fit_event_timings(local_all, local_sectors, p, epoch, shape)
    pooled_secondaries = {col: pool_sector_depths(full, col) for col in
                          ('KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG')}
    # Search any other orbital phase after removing a generous primary window.
    # This addresses eccentric eclipsing binaries whose secondary is not at 0.5.
    primary_width = max(duration, shape['duration_days'] if shape else duration)
    keep_secondary = abs(phase(t, p, epoch)) > 1.25 * primary_width
    other = BoxLeastSquares(t[keep_secondary], f[keep_secondary], e[keep_secondary]).power([p], duration_grid,
                                                                                      oversample=15, objective='snr')
    other_epoch = float(other.transit_time[0])
    other_phase = float((other_epoch - epoch) % p / p)
    other_stats = event_stats(t, f, e, p, other_epoch, float(other.duration[0]))
    secondary_phase_search = dict(epoch_btjd=other_epoch, phase=other_phase, duration_days=float(other.duration[0]),
                                  diagonal_snr=float(other.depth_snr[0]), local_event_stats=other_stats,
                                  note='Targeted additional-phase ranking; no trials/false-alarm claim; raw local verification required.')
    for report, (row, d, good, flat) in zip(full, curves):
        st, se = d['TIME'][good], d['KSPSAP_FLUX_ERR'][good]
        for col in ('KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG'):
            stats, _ = local_quadratic(st, d[col][good], se, p, other_epoch, float(other.duration[0]))
            report['column_tests'][col]['additional_phase_secondary'] = stats
            report['column_tests'][col]['variability'] = variability_diagnostics(
                st, d[col][good], se, p, epoch, primary_width,
                additional_epoch=other_epoch, additional_duration=float(other.duration[0]))
    pooled_other_phase = {col: pool_sector_depths(full, col, event='additional_phase_secondary') for col in
                          ('KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG')}
    secondary_counterevidence = (pooled_other_phase['SAP_FLUX'] is not None and
                                pooled_other_phase['KSPSAP_FLUX_SML'] is not None and
                                pooled_other_phase['SAP_FLUX']['snr'] > 5 and
                                pooled_other_phase['KSPSAP_FLUX_SML']['snr'] > 5)
    planetary_assessment = dict(status='eclipsing_stellar_source_favored' if secondary_counterevidence else 'transit_like_lead_unresolved',
                                reason='Recurring secondary appears in independently normalized raw SAP and alternate aperture flux.' if secondary_counterevidence else
                                       'These diagnostics do not establish planetary origin or novelty.',
                                confirmed=False)
    means = [v['column_tests']['KSPSAP_FLUX']['primary'] for v in full if v['column_tests']['KSPSAP_FLUX']['primary']['n_events'] >= 3]
    if means:
        depths = np.array([v['depth'] for v in means])
        errors = np.array([v['depth_err'] for v in means])
        w = 1 / errors ** 2
        mean = float(np.sum(w * depths) / w.sum())
        scatter_chi2 = float(np.sum(((depths - mean) / errors) ** 2))
        pooled = dict(depth=mean, depth_err=float(1 / np.sqrt(w.sum())), consistency_chi2=scatter_chi2,
                      consistency_dof=len(means) - 1, n_sectors=len(means))
    else:
        pooled = None
    result = dict(tic=tic, checked_utc=datetime.now(timezone.utc).isoformat(), sources=sources,
                  vetting_method=method_identity(),
                  selected_original_peak=peak, selected_period_argument=selected_period, selected_evidence_id=evidence_id,
                  independent_sector_searches=independent,
                  joint_ephemeris=dict(period=p, epoch_btjd=epoch, duration_days=duration,
                                       snr_diagonal=float(power.depth_snr[j]), period_grid_step=float(periods[1] - periods[0]),
                                       baseline_days=float(np.ptp(t)), grid_interval=[float(periods[0]), float(periods[-1])]),
                  joint_fixed_ephemeris_sector_tests=full, pooled_depth=pooled,
                  approximate_pooled_shape=shape,
                  event_timing_fit=timings,
                  pooled_half_phase_secondaries=pooled_secondaries, additional_phase_search=secondary_phase_search,
                  pooled_additional_phase_secondaries=pooled_other_phase,
                  planetary_assessment=planetary_assessment,
                  limitations=['Targeted local period search has no global false-alarm calibration.',
                               'Pooled shape integrates actual per-sector median exposures; parameter errors are conditional on period.',
                               'Nearby stars can generate coherent recurring eclipses; source localization required.',
                               'No novelty claim from catalog absence.'])
    directory = OUT / f'tic{tic}_multisector'
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 'multisector.json').write_text(json.dumps(result, indent=2, allow_nan=False))
    rawsave.update(joint_t=t, joint_f=f, joint_e=e, joint_periods=periods, joint_snr=np.asarray(power.power))
    np.savez_compressed(directory / 'multisector_arrays.npz', **rawsave)
    fig, axes = plt.subplots(len(curves) + 1, 1, figsize=(11, 2.4 * (len(curves) + 1)), layout='constrained')
    axes[0].plot(periods, power.power, lw=.7)
    axes[0].set(xlabel='Period (days)', ylabel='Joint targeted BLS diagonal S/N', title=f'TIC {tic}; joint P={p:.8f} d')
    for ax, (row, d, good, flat) in zip(axes[1:], curves):
        ph = phase(flat['t'], p, epoch)
        ax.plot(ph, flat['f'], '.', ms=2, color='.5')
        ax.set(xlabel='Phase time (days)', ylabel=f'Sector {row["sector"]}', xlim=(-.30, .30))
    fig.savefig(directory / 'multisector.png', dpi=130)
    plt.close(fig)
    if timings and timings.get('fitted_ephemeris'):
        eph = timings['fitted_ephemeris']
        events = timings['events']
        fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
        for sector in sorted({v['sector'] for v in events}):
            ev = [v for v in events if v['sector'] == sector]
            ts = np.array([v['time_btjd'] for v in ev])
            cycles = np.array([v['cycle'] for v in ev])
            errs = np.array([v['time_error_days'] for v in ev])
            oc = ts - eph['epoch_btjd'] - eph['period'] * cycles
            ax.errorbar(ts, oc * 1440, yerr=errs * 1440, fmt='o', ms=3, label=f'S{sector}')
        ax.axhline(0, color='.5', lw=.7)
        ax.set(xlabel='BTJD', ylabel='Observed minus linear ephemeris (minutes)', title=f'TIC {tic}; conditional event timings')
        ax.legend(ncol=4, fontsize=8)
        fig.savefig(directory / 'event_timing.png', dpi=130)
        plt.close(fig)
    lines = [f'# TIC {tic} independent multi-sector check', '',
             f"Assessment: **{planetary_assessment['status']}**. {planetary_assessment['reason']}", '',
             f'Joint local period: {p:.8f} d, BTJD epoch {epoch:.7f}, box duration {duration * 24:.2f} h.', '',
             'Each sector was first independently searched near the sector-1 period, with an unrestricted local transit epoch. '
             'Then the joint fixed ephemeris was tested in raw SAP, processed flux, and alternative apertures.', '',
             '|Sector|Independent P (d)|Local primary S/N|Fixed joint primary S/N|Fixed joint depth (ppm)|Secondary S/N|',
             '|---|---:|---:|---:|---:|---:|']
    for a, b in zip(independent, full):
        pr = b['column_tests']['KSPSAP_FLUX']['primary']
        sec = b['column_tests']['KSPSAP_FLUX']['secondary']
        lines.append(f"|{a['sector']}|{a['ephemeris']['period']:.6f}|{a['primary']['snr']:.2f}|{pr['snr']:.2f}|{(pr['depth'] or 0) * 1e6:.0f}|{sec['snr']:.2f}|")
    lines += ['', 'This establishes or falsifies recurrence in independent observations, not planetary source attribution. '
              'The bright star CD-79 151 at about 86.6 arcsec remains a possible contaminant pending pixel analysis. '
              'All measurements, raw-source hashes, event depths, and period scores are retained.']
    if timings and timings.get('fitted_ephemeris'):
        eph = timings['fitted_ephemeris']
        lines += ['', f"Conditional timing fit from {timings['n_events']} events: P={eph['period']:.9f} "
                  f"±{eph['period_error_days']:.9f} d, T0={eph['epoch_btjd']:.7f} ±{eph['epoch_error_days']:.7f} BTJD. "
                  f"Timing chi²/dof={eph['timing_chi2']:.1f}/{eph['timing_dof']}. Errors are conditional on "
                  'fixed trapezoid shape and inflated for measured scatter; they are not validation probabilities.']
    if pooled_secondaries['KSPSAP_FLUX']:
        sec = pooled_secondaries['KSPSAP_FLUX']
        lines += ['', f"Pooled half-phase secondary: {sec['depth'] * 1e6:.0f} ±{sec['depth_err'] * 1e6:.0f} ppm, "
                  f"conditional S/N {sec['snr']:.2f}, sector-consistency chi²/dof "
                  f"{sec['consistency_chi2']:.1f}/{sec['consistency_dof']}. "
                  'A significant recurring secondary is evidence favoring an eclipsing stellar source; source localization '
                  'and physical modeling are still required.']
    if pooled_other_phase['KSPSAP_FLUX']:
        sec = pooled_other_phase['KSPSAP_FLUX']
        rawsec = pooled_other_phase['SAP_FLUX']
        lines += ['', f"Additional-phase secondary at phase {other_phase:.5f}: raw local processed-flux depth "
                  f"{sec['depth'] * 1e6:.0f} ±{sec['depth_err'] * 1e6:.0f} ppm, S/N {sec['snr']:.2f}; "
                  f"raw SAP S/N {rawsec['snr']:.2f}. This is a targeted phase search, so its score is not a false-alarm probability. "
                  'The catalog-unlisted recurring primary must be judged together with this counterevidence.']
    (directory / 'README.md').write_text('\n'.join(lines) + '\n')
    print('MULTISECTOR', tic, 'P', p, 'epoch', epoch, 'joint_diagonal_snr', float(power.depth_snr[j]), flush=True)
    for a, b in zip(independent, full):
        pr = b['column_tests']['KSPSAP_FLUX']['primary']
        print('sector', a['sector'], 'independent_period', a['ephemeris']['period'], 'independent_snr', a['primary']['snr'],
              'joint_local_snr', pr['snr'], 'joint_depth', pr['depth'], flush=True)
    return result


def validate_method():
    """Check baseline/shape recovery against independently constructed exposure data."""
    rng = np.random.default_rng(502)
    t = np.arange(0, 27., 10 / 1440)
    period, epoch, duration, ingress, depth = 3.713, 1.373, .105, .025, .006
    ph = phase(t, period, epoch)
    # Construct reference exposure averages numerically from the piecewise ramp,
    # independently of the analytic integration routine used in the fitted model.
    offsets = np.linspace(-5 / 1440, 5 / 1440, 1001)
    sample = ph[:, None] + offsets[None, :]
    pulse = np.mean(np.clip(np.minimum((sample + duration / 2) / ingress,
                                     (duration / 2 - sample) / ingress), 0, 1), axis=1)
    f = (1 + .00002 * t + .000002 * t * t) * (1 - depth * pulse) + rng.normal(0, .0004, len(t))
    e = np.full(len(t), .0004)
    stats, data = local_quadratic(t, f, e, period, epoch, duration)
    shape = fit_shape(data, period, epoch, duration, 10 / 1440)
    x = np.linspace(-.12, .12, 49)
    numerical = []
    for value in x:
        s = np.linspace(value - 1 / 96, value + 1 / 96, 50001)
        numerical.append(np.mean(np.clip(np.minimum((s + duration / 2) / ingress,
                                                   (duration / 2 - s) / ingress), 0, 1)))
    integration_error = float(np.max(abs(np.asarray(numerical) - integrated_trapezoid(x, duration, ingress, 1 / 48))))
    checks = dict(exposure_integration_error_below_5e6=integration_error < 5e-6,
                  injected_depth_within_3sigma=abs(shape['depth'] - depth) < 3 * shape['depth_err'],
                  injected_duration_within_3sigma=abs(shape['duration_days'] - duration) < 3 * shape['duration_err'],
                  injected_epoch_within_3sigma=abs(shape['epoch_btjd'] - epoch) < 3 * shape['epoch_err'],
                  seven_independent_events=stats['n_events'] == 7)
    # Verify that robust *raw* fitting recovers ellipsoidal light even with
    # transit/eclipses, positive outliers, and incomplete phase coverage.
    tv = np.arange(0, 27., 30 / 1440)
    pv, ev, dv = 1.626, 1.173, .11
    angle = 2 * np.pi * (tv - ev) / pv
    fv = 1 + .001 * np.sin(angle) + .0025 * np.cos(2 * angle)
    fv -= .005 * (abs(phase(tv, pv, ev)) < dv / 2)
    fv += rng.normal(0, .0008, len(tv))
    fv[::137] += .03
    vv = variability_diagnostics(tv, fv, np.full(len(tv), .0008), pv, ev, dv)
    checks['raw_ellipsoidal_amplitude_recovered'] = abs(vv['half_period_semiamplitude'] - .0025) < 3 * vv['half_period_semiamplitude_err']
    checks['raw_orbital_amplitude_recovered'] = abs(vv['orbital_semiamplitude'] - .001) < 3 * vv['orbital_semiamplitude_err']
    checks['positive_outliers_robustly_rejected'] = vv['n_robust_rejected'] > 0
    checks['injected_ellipsoidal_exceeds_guard_threshold'] = vv['half_period_snr'] > 5 and vv['half_period_semiamplitude'] > .5 * .004
    transit_only = variability_diagnostics(t, f, e, period, epoch, duration)
    checks['transit_only_does_not_raise_raw_guard'] = not any(
        (transit_only[f'{name}_snr'] or 0) > 5 and transit_only[f'{name}_semiamplitude'] > .5 * depth
        for name in ('orbital', 'half_period'))
    payload = dict(checked_utc=datetime.now(timezone.utc).isoformat(), seed=502,
                   injected=dict(period=period, epoch_btjd=epoch, duration_days=duration,
                                 ingress_days=ingress, depth=depth, noise_sigma=.0004,
                                 exposure_minutes=10, baseline='quadratic'),
                   recovered=shape, n_events=stats['n_events'], checks=checks,
                   raw_ellipsoidal_control=vv, transit_only_raw_harmonics=transit_only,
                   maximum_numeric_exposure_error=integration_error,
                   limitations='A synthetic conditional-model check, not real-data completeness or false-positive calibration.')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'method_validation.json').write_text(json.dumps(payload, indent=2))
    lines = ['# Independent broad-lead vetting method', '',
             'Run `python src/broad_vet.py --tic TIC_ID` for retained sector-1 leads. '
             'Use `--period PERIOD` to select a saved harmonic for one TIC. '
             'Use `--multisector --tic TIC_ID` after `src/followup_qlp.py` has saved a follow-up manifest. '
             'Reproduce the numerical control with `python src/broad_vet.py --validate`.', '',
             'Raw FITS source hashes and available archive checksum fields are checked before use. '
             'Local period refinement stays within three broad-search grid cells, preserving the original peak. '
             'Quality-zero flux is examined with 0.5/1/2-day gap-aware medians and independent local quadratic '
             'out-of-transit fits. SAP, processed, small-aperture, and large-aperture flux are compared. '
             'Local depths include residual-noise and event-scatter error floors. Primary and half-phase secondary '
             'events, odd/even differences, broad variability, and QLP centroids provide warnings.', '',
             'Version 3 also fits raw SAP orbital and half-period harmonics outside generous primary, half-phase, '
             'and selected additional-secondary windows. It uses three iterations of 5-sigma clipping, '
             'an observed residual MAD error floor, and reduced-chi-square covariance inflation. All five '
             'sin/cos/intercept coefficients and their covariance are saved. A raw harmonic with conditional '
             'S/N >5 and amplitude >50% of the measured primary depth raises a variability warning, even '
             'if KSPSAP preprocessing removed the modulation. This warns about stellar/ellipsoidal light '
             'and is not an automatic binary confirmation; spot modulation or instrumental trends remain possible.', '',
             'The fitted trapezoid integrates each sector\'s actual median exposure. A fixed-period shape fit is '
             'conditional and approximate, with no limb darkening. Local box-mean depths differ from fitted '
             'full-depth values because ingress and exposure smear the transit. Independent event timing fits '
             'include scatter inflation and leave-one-sector-out stability checks.', '',
             f'Numerical control: maximum analytic-vs-numeric exposure-average error {integration_error:.3g}; '
             f'injected depth {depth:.6f}, fitted {shape["depth"]:.6f} ±{shape["depth_err"]:.6f}; '
             f'injected duration {duration:.6f} d, fitted {shape["duration_days"]:.6f} ±{shape["duration_err"]:.6f} d. '
             f'All checks passed: {all(checks.values())}. The numerical injection used a quadratic baseline, '
             'seven recurrences, 10-minute integrated samples, and white noise, and is not a completeness estimate.', '',
             'Aperture/centroid warnings are source-contamination diagnostics, not automatic proof of a false positive. '
             'An unexplained lead still requires pixel localization, stellar/transit physical models, comprehensive '
             'catalog and literature checks, and empirical survey-wide correlated-noise/trials calibration. '
             'No folder here establishes planetary discovery or confirmation.']
    (OUT / 'METHOD.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps(checks))
    if not all(checks.values()):
        raise RuntimeError('Scientific method validation failed')


def _ranked_summary(result, priority=None, evidence_id=None):
    primary = result['column_tests']['KSPSAP_FLUX']['local_quadratic']
    raw = result['column_tests']['SAP_FLUX']['local_quadratic']
    additional = result['column_tests']['KSPSAP_FLUX']['additional_phase_secondary']
    p = result['refined_ephemeris']['period']
    shape = result.get('integrated_trapezoid')
    duration = shape['duration_days'] if shape else result['refined_ephemeris']['duration_days']
    score = max(0., min(primary['snr'], raw['snr'], 30.)) * max(.05, 1 - duration / p / .2)
    flags = result['warnings']
    if any(v in flags for v in ('significant_half_phase_secondary', 'significant_additional_phase_secondary',
                                'weak_independent_raw_local_detection', 'fewer_than_three_independent_positive_events')):
        score *= .025
    if any(v in flags for v in ('aperture_dependent_depth', 'qlp_centroid_shift_warning')):
        score *= .5
    if any(v in flags for v in ('odd_even_depth_difference', 'v_shaped_or_unresolved_ingress',
                                'fitted_duty_cycle_exceeds_point15', 'large_orbital_phase_variability',
                                'large_significant_raw_orbital_variability', 'large_significant_raw_half_period_variability')):
        score *= .25
    directory = OUT / f'tic{result["tic"]}_p{result["input_peak"]["period"]:.6f}'
    return dict(tic=result['tic'], source_evidence_id=evidence_id, priority_source=priority,
                ranking_score_unvalidated=score, period_days=p, duration_hours=duration * 24,
                local_depth_ppm=(primary['depth'] or 0.) * 1e6, fitted_depth_ppm=shape['depth'] * 1e6 if shape else None,
                processed_local_snr=primary['snr'], raw_sap_local_snr=raw['snr'],
                n_events=primary['n_events'], n_positive_events=primary.get('n_positive_events', 0),
                half_phase_secondary_snr=result['column_tests']['KSPSAP_FLUX']['secondary_half_phase']['snr'],
                other_phase=result['additional_phase_search']['phase'], other_phase_secondary_snr=additional['snr'],
                odd_even_sigma=primary['odd_even_sigma'], ingress_fraction=shape['ingress_fraction'] if shape else None,
                centroid_max_sigma=max((v['sigma'] or 0. for v in result['centroid_tests'].values()), default=0.),
                raw_orbital_semiamplitude_ppm=(result['column_tests']['SAP_FLUX']['variability']['orbital_semiamplitude'] or 0.) * 1e6,
                raw_half_period_semiamplitude_ppm=(result['column_tests']['SAP_FLUX']['variability']['half_period_semiamplitude'] or 0.) * 1e6,
                raw_orbital_snr=result['column_tests']['SAP_FLUX']['variability'].get('orbital_snr'),
                raw_half_period_snr=result['column_tests']['SAP_FLUX']['variability'].get('half_period_snr'),
                assessment=result['assessment'], warnings=';'.join(flags),
                evidence_file=str((directory / 'vetting.json').relative_to(ROOT)),
                method_sha256=result['vetting_method']['source_sha256'])


def export_vetting_summary(results, checkpoint):
    rows = sorted(results, key=lambda v: -v['ranking_score_unvalidated'])
    OUT.mkdir(parents=True, exist_ok=True)
    atomic_json(OUT / 'batch_checkpoint.json', dict(checkpoint, completed=len(rows), results=rows))
    if rows:
        destination = OUT / 'ranked_vetted_leads.csv'
        temp = destination.with_suffix('.csv.tmp')
        with temp.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        temp.replace(destination)


def preserve_previous_batch():
    """Retain exact prior selection/checkpoint/ranking and compact dossier state."""
    checkpoint_path = OUT / 'batch_checkpoint.json'
    if not checkpoint_path.exists():
        return
    raw = checkpoint_path.read_bytes()
    previous = json.loads(raw)
    version = previous.get('method', {}).get('version', 'unversioned')
    digest = hashlib.sha256(raw).hexdigest()
    destination = OUT / 'batch_history' / f'{version}_{digest[:12]}'
    destination.mkdir(parents=True, exist_ok=True)
    for name in ('batch_checkpoint.json', 'batch_queue.csv', 'ranked_vetted_leads.csv'):
        source = OUT / name
        if source.exists() and not (destination / name).exists():
            (destination / name).write_bytes(source.read_bytes())
    dossiers = destination / 'vetting_reports.jsonl.gz'
    if not dossiers.exists():
        with dossiers.open('wb') as stream, gzip.GzipFile(fileobj=stream, mode='wb', mtime=0) as compressed:
            for result in previous.get('results', []):
                source = ROOT / result['evidence_file']
                if source.exists():
                    compressed.write(json.dumps(json.loads(source.read_text()), separators=(',', ':')).encode() + b'\n')
    atomic_json(destination / 'snapshot_manifest.json', dict(previous_checkpoint_sha256=digest,
                                                            previous_method=previous.get('method'),
                                                            completed=previous.get('completed'),
                                                            queue_sha256=previous.get('queue_sha256'),
                                                            raw_data_deleted=False))


def run_queue(queue_path, top):
    """Freeze a review queue, honor its exact harmonic, and checkpoint every star."""
    raw = Path(queue_path).read_bytes()
    # Read source first: it may itself be the batch_queue.csv to be archived.
    preserve_previous_batch()
    rows = list(csv.DictReader(raw.decode().splitlines()))
    eligible = [r for r in rows if r['category'] == 'unrecognized_lead' and float(r['period_days']) >= .5]
    eligible.sort(key=lambda r: -float(r['priority_score_unvalidated']))
    selected, seen = [], set()
    for row in eligible:
        tic = int(row['tic'])
        if tic in seen:
            continue
        seen.add(tic)
        selected.append(row)
        if len(selected) >= top:
            break
    identity = method_identity()
    checkpoint = dict(started_utc=datetime.now(timezone.utc).isoformat(), queue_path=str(queue_path),
                      queue_sha256=hashlib.sha256(raw).hexdigest(), method=identity,
                      eligible_rows=len(eligible), selected_unique_stars=len(selected), requested_top=top,
                      planned_evidence_ids=[r['evidence_id'] for r in selected], errors=[])
    OUT.mkdir(parents=True, exist_ok=True)
    # Preserve this precise batch selection in ordinary Git; raw data stay cached.
    snapshot = OUT / 'batch_queue.csv'
    temp = snapshot.with_suffix('.csv.tmp')
    with temp.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(selected)
    temp.replace(snapshot)
    db = sqlite3.connect(f'file:{DATABASE}?mode=ro', uri=True)
    summaries = []
    for row in selected:
        tic = int(row['tic'])
        raw_record = db.execute('SELECT result_json FROM stars WHERE tic=?', (tic,)).fetchone()
        try:
            record = json.loads(raw_record[0])
            index = row['evidence_id'].split('/')
            scan_index, peak_index = int(index[2]), int(index[4])
            scan, peak = record['scans'][scan_index], record['scans'][scan_index]['peaks'][peak_index]
            if record['sha256'] != row['raw_sha256'] or abs(peak['period'] - float(row['period_days'])) > 1e-9:
                raise ValueError('Queue evidence identifier/hash/ephemeris no longer matches source')
            destination = OUT / f'tic{tic}_p{peak["period"]:.6f}' / 'vetting.json'
            cached = json.loads(destination.read_text()) if destination.exists() else None
            current = (cached is not None and cached.get('vetting_method') == identity and
                       cached.get('source', {}).get('sha256') == record['sha256'] and
                       cached.get('additional_phase_search') is not None and
                       all('additional_phase_secondary' in cached.get('column_tests', {}).get(k, {}) for k in
                           ('KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG')))
            if current:
                result = cached
                print('CACHED', tic, row['period_days'], flush=True)
            else:
                result = vet(record, peak, scan['label'])
            summaries.append(_ranked_summary(result, row['priority_score_unvalidated'], row['evidence_id']))
        except Exception as err:
            checkpoint['errors'].append(dict(tic=tic, evidence_id=row['evidence_id'],
                                             error=type(err).__name__ + ': ' + str(err)))
            print('ERROR', tic, type(err).__name__, str(err), flush=True)
        checkpoint['updated_utc'] = datetime.now(timezone.utc).isoformat()
        export_vetting_summary(summaries, checkpoint)
    db.close()
    print('DETAILED_BATCH_COMPLETE', len(summaries), 'errors', len(checkpoint['errors']), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tic', type=int, nargs='+')
    parser.add_argument('--top', type=int, default=12)
    parser.add_argument('--period', type=float, help='Select nearest saved broad peak, for one --tic target only.')
    parser.add_argument('--multisector', action='store_true', help='Analyze cached follow-up manifest for one target.')
    parser.add_argument('--evidence-id', help='Exact saved evidence identifier for one --multisector target.')
    parser.add_argument('--validate', action='store_true', help='Run independent numeric/synthetic method controls.')
    parser.add_argument('--queue', type=Path, default=ROOT / 'results/broad/shortlist.csv', help='Use exact ranked evidence rows for --top.')
    args = parser.parse_args()
    if args.validate:
        validate_method()
        return
    if args.period is not None and (not args.tic or len(args.tic) != 1):
        parser.error('--period requires exactly one --tic')
    if args.multisector:
        if not args.tic or len(args.tic) != 1:
            parser.error('--multisector requires exactly one --tic')
        vet_multisector(args.tic[0], selected_period=args.period, evidence_id=args.evidence_id)
        return
    if not args.tic and args.queue.exists():
        run_queue(args.queue, args.top)
        return
    leads = []
    for record in load_records(args.tic):
        for scan in record['scans']:
            for peak in scan['peaks']:
                if args.tic:
                    selected = min(scan['peaks'], key=lambda v: abs(v['period'] - args.period)) if args.period else scan['peaks'][0]
                    if peak == selected:
                        leads.append((peak['local_snr'], record, peak, scan['label']))
                    continue
                cat = peak['catalog']
                if not peak['screening_pass'] or peak['n_events'] < 3 or peak['local_snr'] < 10:
                    continue
                if peak['period'] < 1 or peak['depth'] > .015:
                    continue
                if cat['known_period_match'] or cat['known_harmonic_match'] or record['catalog_flags']['is_eb_host']:
                    continue
                radius = record['metadata'].get('RADIUS')
                if radius and radius > 0 and np.sqrt(peak['depth']) * radius * 109.076 > 20:
                    continue
                leads.append((peak['local_snr'], record, peak, scan['label']))
    done = set()
    results = []
    for _, record, peak, label in sorted(leads, key=lambda v: -v[0]):
        if record['tic'] in done:
            continue
        done.add(record['tic'])
        results.append(vet(record, peak, label))
        if not args.tic and len(done) >= args.top:
            break
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'last_batch_summary.json').write_text(json.dumps([dict(tic=r['tic'], period=r['input_peak']['period'],
                                                                assessment=r['assessment'], warnings=r['warnings']) for r in results], indent=2))


if __name__ == '__main__':
    main()
