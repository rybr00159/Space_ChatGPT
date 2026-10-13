"""Broad, inexpensive TESS screening; statistics rank leads, not discoveries.

The period grid is logarithmic, so the nearest grid point has an accumulated
phase drift <= half a shortest trial duration over the observed baseline. Top
coarse peaks are refined, then subjected to event-count/duty-cycle guards.
Independent local box events retain single-transit leads outside periodic BLS.
Use raw/pixel data, alternative detrending, catalogs, and search-null calibration
before interpreting any reported S/N as planetary evidence.
"""
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import time

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter


@dataclass(frozen=True)
class FastScreenConfig:
    min_period: float = .3
    max_period: float = 13.
    durations: tuple = (.025, .05, .10, .16)
    n_peaks: int = 5
    n_refine: int = 20
    compute_full_stats: bool = True
    scan_single_events: bool = True
    grid_drift_durations: float = 1.
    coarse_oversample: int = 5
    refine_oversample: int = 10
    bin_minutes: float = 6.
    width_days: float = 1.
    gap_days: float = 1 / 48
    edge_cut_days: float = .06
    max_duty: float = .15
    single_event_snr: float = 5.


def _config(config):
    return config if isinstance(config, FastScreenConfig) else FastScreenConfig(**(config or {}))


def _mad(x):
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else 0.


def phase(t, period, epoch):
    return (t - epoch + .5 * period) % period - .5 * period


def known_mask(t, known):
    """Mask catalog ephemerides in BTJD; missing epochs never invent a mask."""
    keep = np.ones(len(t), dtype=bool)
    for row in known or []:
        p, epoch, duration = row.get('period'), row.get('epoch_bjd'), row.get('duration_hours')
        try:
            p, epoch, duration = float(p), float(epoch), float(duration or 3.)
        except (ValueError, TypeError):
            continue
        if not np.isfinite(p) or p <= 0 or not np.isfinite(epoch):
            continue
        if not np.isfinite(duration) or duration <= 0:
            duration = 3.
        keep &= np.abs(phase(t, p, epoch - 2457000)) > max(duration / 24 * .85, .05)
    return keep


def preprocess_curve(t, flux, error, quality=None, known=None, config=None):
    """Quality filter, gap-safe weighted bins, and robust median normalization.

    Gaps >=max(30 min, 2.5 cadences) define independent detrending segments. Bins are relative to
    each segment, so no sample is combined across a gap. Cadences near segment
    edges are excluded; this explicitly sacrifices boundary completeness.
    Known-event samples are removed before binning and trend estimation. Noise
    errors never fall below robust observed scatter after detrending.
    """
    cfg = _config(config)
    t, flux, error = (np.asarray(v, dtype=float) for v in (t, flux, error))
    if not (t.shape == flux.shape == error.shape) or t.ndim != 1:
        raise ValueError('time, flux, error must be equal-length one-dimensional arrays')
    good = np.isfinite(t) & np.isfinite(flux) & np.isfinite(error) & (flux > 0) & (error > 0)
    if quality is not None:
        q = np.asarray(quality)
        if q.shape != t.shape:
            raise ValueError('quality must have the same shape as time')
        good &= q == 0
    t, flux, error = t[good], flux[good], error[good]
    order = np.argsort(t)
    t, flux, error = t[order], flux[order], error[order]
    n_clean = len(t)
    if n_clean < 20:
        return dict(t=np.array([]), f=np.array([]), e=np.array([]), segment=np.array([], dtype=int),
                    n_clean=n_clean, n_known_masked=0, n_segments=0, cadence_days=None, effective_gap_days=None)
    # Segment on original quality-filtered cadence, never on the catalog mask.
    cadence = float(np.median(np.diff(t)[np.diff(t) > 0]))
    effective_gap = max(cfg.gap_days, 2.5 * cadence)
    starts = np.r_[0, np.flatnonzero(np.diff(t) >= effective_gap) + 1, len(t)]
    km = known_mask(t, known)
    chunks = []
    for segment, (a, b) in enumerate(zip(starts[:-1], starts[1:])):
        ts, fs, es = t[a:b], flux[a:b], error[a:b]
        use = km[a:b]
        if use.sum() < 20 or ts[-1] - ts[0] < 2 * cfg.edge_cut_days:
            continue
        original_lo, original_hi = ts[0], ts[-1]
        ts, fs, es = ts[use], fs[use], es[use]
        idx = np.floor((ts - original_lo) / (cfg.bin_minutes / 1440)).astype(int)
        _, inv = np.unique(idx, return_inverse=True)
        weight = 1 / es ** 2
        sw = np.bincount(inv, weights=weight)
        tb = np.bincount(inv, weights=weight * ts) / sw
        fb = np.bincount(inv, weights=weight * fs) / sw
        eb = 1 / np.sqrt(sw)
        dt = max(cfg.bin_minutes / 1440, float(np.median(np.diff(tb))))
        size = max(3, int(cfg.width_days / dt) // 2 * 2 + 1)
        trend = median_filter(fb, size=size, mode='reflect')
        valid = (trend > 0) & np.isfinite(trend)
        fb, eb = fb / trend, eb / trend
        scatter = _mad(fb[valid])
        # Remove positive flares/cosmic rays only: preserve negative event evidence.
        use = valid & (fb < 1 + max(7 * scatter, .001))
        use &= (tb >= original_lo + cfg.edge_cut_days) & (tb <= original_hi - cfg.edge_cut_days)
        if use.sum() < 12:
            continue
        # Conservative diagonal errors; neither red-noise nor a trials correction.
        floor = max(_mad(fb[use]), np.median(eb[use]))
        eb = np.maximum(eb, floor)
        chunks.append((tb[use], fb[use], eb[use], np.full(use.sum(), segment, dtype=int)))
    arrays = [np.concatenate([c[j] for c in chunks]) for j in range(4)] if chunks else [np.array([])] * 4
    return dict(zip(('t', 'f', 'e', 'segment'), arrays), n_clean=n_clean,
                n_known_masked=int((~km).sum()), n_segments=len(starts) - 1,
                cadence_days=cadence, effective_gap_days=effective_gap)


def period_grid(baseline, config=None):
    """Log grid with explicitly bounded period mismatch; no uniform-frequency undersampling."""
    cfg = _config(config)
    maximum = min(cfg.max_period, baseline / 2)
    if maximum <= cfg.min_period:
        return np.array([])
    step = min(cfg.durations) * cfg.grid_drift_durations / baseline
    periods = cfg.min_period * np.exp(np.arange(0, np.log(maximum / cfg.min_period) + step, step))
    return np.unique(np.r_[periods[periods <= maximum], maximum])


def _peak_seeds(periods, scores, baseline, n):
    available = np.isfinite(scores)
    chosen = []
    for _ in range(n):
        if not available.any():
            break
        j = int(np.argmax(np.where(available, scores, -np.inf)))
        chosen.append(j)
        # Remove the same narrow peak, retaining explicitly reviewable aliases.
        available &= np.abs(1 / periods - 1 / periods[j]) > .05 / baseline
        available &= np.abs(periods / periods[j] - 1) > .003
    return chosen


def event_stats(t, f, e, period, epoch, duration, include_details=True):
    """Event depths use local weighted baselines and exclude unobserved transits."""
    cycle = np.round((t - epoch) / period).astype(int)
    inside = np.abs(phase(t, period, epoch)) <= duration / 2
    events = []
    for k in np.unique(cycle[inside]):
        center = epoch + k * period
        ii = inside & (cycle == k)
        oo = (np.abs(t - center) < 2.5 * duration) & (np.abs(t - center) > duration)
        left = oo & (t < center)
        right = oo & (t > center)
        cadence = float(np.median(np.diff(t)))
        min_inside = max(2, min(3, int(duration / cadence * .7)))
        if ii.sum() < min_inside or left.sum() < 2 or right.sum() < 2:
            continue
        wi, wo = 1 / e[ii] ** 2, 1 / e[oo] ** 2
        depth = np.sum(wo * f[oo]) / wo.sum() - np.sum(wi * f[ii]) / wi.sum()
        err = np.sqrt(1 / wo.sum() + 1 / wi.sum())
        events.append(dict(cycle=int(k), time=float(center), n=int(ii.sum()), depth=float(depth),
                           depth_err=float(err), snr=float(depth / err)))
    if not events:
        return dict(n_events=0, n_positive_events=0, events=[], event_snr_max_fraction=None,
                    odd_even_sigma=None, local_snr=0., local_depth=None, local_depth_err=None)
    ds, es = np.array([v['depth'] for v in events]), np.array([v['depth_err'] for v in events])
    w = 1 / es ** 2
    depth, err = np.sum(w * ds) / w.sum(), 1 / np.sqrt(w.sum())
    odd_even = None
    parity = np.array([v['cycle'] % 2 for v in events])
    if (parity == 0).any() and (parity == 1).any():
        means = [np.sum(w[parity == k] * ds[parity == k]) / w[parity == k].sum() for k in (0, 1)]
        odd_even = float(abs(means[0] - means[1]) / np.sqrt(sum(1 / w[parity == k].sum() for k in (0, 1))))
    power = np.maximum(ds / es, 0) ** 2
    return dict(n_events=len(events), n_positive_events=int((ds > 0).sum()), events=events if include_details else [],
                event_snr_max_fraction=float(power.max() / power.sum()) if power.sum() else None,
                odd_even_sigma=odd_even, local_snr=float(depth / err), local_depth=float(depth),
                local_depth_err=float(err))


def local_events(t, f, e, segments, config=None):
    """Linear-time isolated-box scan with local sidebands; retains rare-transit leads."""
    cfg = _config(config)
    leads = []
    for segment in np.unique(segments):
        use = segments == segment
        ts, fs, es = t[use], f[use], e[use]
        if len(ts) < 20:
            continue
        w = 1 / es ** 2
        sw, sf = np.r_[0., np.cumsum(w)], np.r_[0., np.cumsum(w * fs)]
        for duration in cfg.durations:
            limits = [np.searchsorted(ts, ts + offset * duration) for offset in (-2.5, -1., -.5, .5, 1., 2.5)]
            a, b, c, d, g, h = limits
            wi, wo = sw[d] - sw[c], sw[b] - sw[a] + sw[h] - sw[g]
            min_inside = max(2, min(3, int(duration / np.median(np.diff(ts)) * .7)))
            valid = (d - c >= min_inside) & (b - a >= 2) & (h - g >= 2) & (wi > 0) & (wo > 0)
            valid &= (ts >= ts[0] + 2.5 * duration) & (ts <= ts[-1] - 2.5 * duration)
            score = np.full(len(ts), -np.inf)
            depth = np.zeros(len(ts))
            depth[valid] = (sf[b[valid]] - sf[a[valid]] + sf[h[valid]] - sf[g[valid]]) / wo[valid] - (sf[d[valid]] - sf[c[valid]]) / wi[valid]
            err = np.full(len(ts), np.inf)
            err[valid] = np.sqrt(1 / wi[valid] + 1 / wo[valid])
            score[valid] = depth[valid] / err[valid]
            available = valid.copy()
            for _ in range(5):
                if not available.any():
                    break
                j = int(np.argmax(np.where(available, score, -np.inf)))
                if score[j] < cfg.single_event_snr:
                    break
                leads.append(dict(time=float(ts[j]), duration_days=float(duration), depth=float(depth[j]),
                                  depth_err=float(err[j]), snr=float(score[j]), n=int(d[j] - c[j]), segment=int(segment)))
                available &= np.abs(ts - ts[j]) > 3 * duration
    selected = []
    for lead in sorted(leads, key=lambda v: -v['snr']):
        if all(abs(lead['time'] - prev['time']) > max(lead['duration_days'], prev['duration_days']) for prev in selected):
            selected.append(lead)
        if len(selected) >= 10:
            break
    return selected


def screen_curve(t, flux, error, quality=None, known=None, config=None):
    """Return compact JSON-safe ranked periodic/isolated-event screening results."""
    begin = time.perf_counter()
    cfg = _config(config)
    data = preprocess_curve(t, flux, error, quality=quality, known=known, config=cfg)
    t, f, e = data['t'], data['f'], data['e']
    # FFI data may have 30-minute exposures. Never search sub-resolution boxes;
    # this reduces short-transit completeness, which is recorded explicitly.
    effective_cadence = max(cfg.bin_minutes / 1440, data['cadence_days'] or 0.)
    resolved_durations = tuple(d for d in cfg.durations if d >= 1.5 * effective_cadence)
    if not resolved_durations:
        resolved_durations = (2 * effective_cadence, 3 * effective_cadence)
    cfg = replace(cfg, durations=resolved_durations)
    out = dict(status='insufficient_data', n_input=len(flux), n_clean=data['n_clean'], n_searched=len(t),
               n_known_masked=data['n_known_masked'], n_segments=data['n_segments'], peaks=[], single_events=[],
               baseline_days=float(np.ptp(t)) if len(t) else 0.,
               cadence_days=data['cadence_days'], effective_gap_days=data['effective_gap_days'])
    periods = period_grid(out['baseline_days'], cfg) if len(t) > 40 else np.array([])
    if not len(periods):
        out['runtime_seconds'] = time.perf_counter() - begin
        return out
    bls = BoxLeastSquares(t, f, e)
    coarse = bls.power(periods, cfg.durations, oversample=cfg.coarse_oversample, objective='snr')
    if cfg.n_refine > 0:
        seeds = _peak_seeds(periods, coarse.power, out['baseline_days'], cfg.n_refine)
        step = min(cfg.durations) * cfg.grid_drift_durations / out['baseline_days']
        refined_periods = np.unique(np.concatenate([periods[j] * np.exp(np.linspace(-step, step, 25)) for j in seeds]))
        refined_periods = refined_periods[(refined_periods >= cfg.min_period) & (refined_periods <= min(cfg.max_period, out['baseline_days'] / 2))]
        refined = bls.power(refined_periods, cfg.durations, oversample=cfg.refine_oversample, objective='snr')
    else:
        refined_periods, refined = periods, coarse
    picked = _peak_seeds(refined_periods, refined.power, out['baseline_days'], cfg.n_peaks)
    med, mad = float(np.median(coarse.power)), _mad(coarse.power)
    for j in picked:
        p, t0, duration = float(refined.period[j]), float(refined.transit_time[j]), float(refined.duration[j])
        stats = event_stats(t, f, e, p, t0, duration, include_details=cfg.compute_full_stats)
        guards = []
        if stats['n_events'] < 2:
            guards.append('fewer_than_two_sampled_events')
        if stats['n_positive_events'] < 2:
            guards.append('fewer_than_two_positive_local_depths')
        if duration / p > cfg.max_duty:
            guards.append('excessive_duty_cycle')
        if stats['event_snr_max_fraction'] is not None and stats['event_snr_max_fraction'] > .85:
            guards.append('single_event_dominates')
        out['peaks'].append(dict(period=p, t0_btjd=t0, duration_days=duration,
                                 depth=float(refined.depth[j]), depth_err=float(refined.depth_err[j]),
                                 snr=float(refined.depth_snr[j]), duty_cycle=duration / p,
                                 sde_robust=float((refined.power[j] - med) / mad) if mad > 0 else None,
                                 screening_pass=not guards, guard_failures=guards, **stats))
    if cfg.scan_single_events:
        out['single_events'] = local_events(t, f, e, data['segment'], cfg)
    out['status'] = 'screened'
    out['search_metadata'] = dict(config=asdict(cfg), period_grid_count=len(periods),
                                  refined_grid_count=len(refined_periods) if cfg.n_refine > 0 else 0,
                                  search_tier='coarse_to_fine' if cfg.n_refine > 0 else 'broad_coarse',
                                  maximum_period_searched=float(periods[-1]),
                                  nearest_grid_phase_drift_bound_days=.5 * min(cfg.durations) * cfg.grid_drift_durations,
                                  effective_cadence_days=effective_cadence,
                                  minimum_duration_searched_days=min(cfg.durations),
                                  trials_corrected=False, false_alarm_probability=None,
                                  limitations='Diagonal errors; no red-noise or global trials calibration; requires independent vetting.')
    out['runtime_seconds'] = float(time.perf_counter() - begin)
    return out


def screen_fits(path, known=None, config=None, column='PDCSAP_FLUX'):
    """SPOC FITS wrapper; caller must verify downloaded artifacts before use."""
    with fits.open(Path(path), memmap=False) as hdul:
        header, d = hdul[1].header, hdul[1].data
        if header.get('TIMESYS') != 'TDB' or abs(header.get('BJDREFI', 0) + header.get('BJDREFF', 0) - 2457000) > 1e-8:
            raise ValueError('Unexpected time system/reference; do not assume BTJD')
        result = screen_curve(d['TIME'], d[column], d[column + '_ERR'], quality=d['QUALITY'], known=known, config=config)
        result['fits_metadata'] = {key: header.get(key, hdul[0].header.get(key)) for key in ('TICID', 'SECTOR', 'CROWDSAP', 'FLFRCSAP', 'PROCVER')}
        return result
