#!/usr/bin/env python3
"""Independent LS and local-polynomial diagnostics for TIC 214311908.

The initial LS assessment used only the four original QLP files, before the
parent investigator supplied the 5.6915-d BLS result and neighboring EB ephemeris.
This script preserves that assessment and then tests those supplied hypotheses.
No BLS detrended arrays or fitted transit fluxes are used here.
"""
from pathlib import Path
import hashlib
import json
import platform
from datetime import datetime, timezone

import astropy
from astropy.io import fits
from astropy.timeseries import LombScargle
from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.optimize import least_squares
import scipy

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "candidates/tic214311908/independent_variability"
OUT.mkdir(parents=True, exist_ok=True)
RECORDS = json.loads((ROOT / "data/tic214311908/observations.json").read_text())
PERIOD = 5.691497928943136
EPOCH = 1632.6146842062005
DURATION = 0.04958333333333333
EB_PERIOD = 0.31608296204421177
EB_EPOCH = Time(2457558.1303872317, format="jd", scale="tcb").tdb.jd - 2457000


def clean_json(x):
    if isinstance(x, dict):
        return {str(k): clean_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean_json(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean_json(x.tolist())
    if isinstance(x, np.generic):
        return clean_json(x.item())
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def mad(x):
    return float(1.4826 * np.median(abs(x - np.median(x))))


def fit_sinusoid(t, y, period):
    m = np.isfinite(t) & np.isfinite(y)
    t, y = t[m], y[m]
    phase = 2 * np.pi * (t - 3000) / period
    X = np.c_[np.ones(len(t)), np.sin(phase), np.cos(phase)]
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X @ b
    cov = np.linalg.inv(X.T @ X) * np.sum(res**2) / (len(y) - 3)
    amp = np.hypot(*b[1:])
    ae = np.sqrt(b[1:] @ cov[1:, 1:] @ b[1:] / amp**2)
    return dict(amplitude=float(amp), error=float(ae), conditional_snr=float(amp / ae),
                sine=float(b[1]), cosine=float(b[2]), phase_rad_at_btjd3000=float(np.arctan2(b[2], b[1])))


def local_event(t, y, epoch, width=.30, duration=DURATION, relative=True):
    m = np.isfinite(t) & np.isfinite(y) & (abs(t - epoch) <= width)
    x, z = t[m] - epoch, y[m]
    if len(x) < 8:
        return dict(status="insufficient_window", n=len(x))
    oot = abs(x) >= duration * .85
    use = oot.copy()
    X = np.c_[np.ones(len(x)), x, x*x]
    if use.sum() < 6:
        return dict(status="insufficient_baseline", n=len(x))
    for _ in range(4):
        b = np.linalg.lstsq(X[use], z[use], rcond=None)[0]
        resid = z - X @ b
        sigma = mad(resid[use])
        newuse = oot & (abs(resid) < 4 * max(sigma, np.finfo(float).eps))
        if newuse.sum() < 6 or np.array_equal(use, newuse):
            break
        use = newuse
    b = np.linalg.lstsq(X[use], z[use], rcond=None)[0]
    resid = z - X @ b
    inside = abs(x) <= duration / 2
    if not inside.any():
        return dict(status="no_in_window_cadences", n=len(x), nearest_days=float(abs(x).min()))
    sigma = max(mad(resid[use]), np.std(resid[use]) * .7)
    cov = np.linalg.inv(X[use].T @ X[use]) * sigma**2
    xmean = X[inside].mean(axis=0)
    err = float(np.sqrt(sigma**2 / inside.sum() + xmean @ cov @ xmean))
    depth = float(-resid[inside].mean())
    baseline = float(np.mean((X @ b)[inside]))
    if relative:
        depth, err = depth / baseline, err / abs(baseline)
    return dict(status="measured", n=len(x), n_inside=int(inside.sum()), n_baseline=int(use.sum()),
                depth=depth, error=err, conditional_snr=depth/err,
                baseline=baseline, residual_sigma=float(sigma),
                depth_definition="positive is reduction; relative to fitted baseline" if relative else "positive is reduction, native units")


def event_shape(t, y, epoch):
    """Exploratory integrated-trapezoid shape with a simultaneous quadratic.

    This is not a physical transit fit. Individual duration/depth fits are
    noisy and their reported diagonal errors are explicitly conditional.
    """
    m = np.isfinite(t) & np.isfinite(y) & (abs(t-epoch) < .30)
    x, z = t[m]-epoch, y[m]
    if len(x) < 12 or not np.any(abs(x) < DURATION/2):
        return dict(status="insufficient_data")
    z = z / np.median(z)-1
    delta = np.diff(t)
    exposure = float(np.median(delta[(delta > 0) & (delta < .05)]))
    offsets = np.linspace(-.5, .5, 15) * exposure
    X = np.c_[np.ones(len(x)), x, x*x]
    baseline = abs(x) > .10
    if baseline.sum() < 8:
        return dict(status="insufficient_baseline")
    beta = np.linalg.lstsq(X[baseline], z[baseline], rcond=None)[0]
    sigma = mad((z-X@beta)[baseline])
    def template(center, duration):
        # Total duration includes a fixed 20%-of-duration ingress/egress.
        distance = abs(x[:,None]+offsets[None,:]-center)
        return np.clip((duration/2-distance) / (duration*.2), 0, 1).mean(axis=1)
    def fun(par):
        return (X@par[:3] - par[3]*template(par[4],par[5]) - z) / sigma
    starts = [np.r_[beta,.018,shift,.06] for shift in [-.03,0,.03]]
    fits_ = [least_squares(fun, p, bounds=([-np.inf]*3+[0,-.10,.015], [np.inf]*3+[.25,.10,.18]),
                           x_scale="jac", max_nfev=500) for p in starts]
    fit = min(fits_, key=lambda a: np.sum(a.fun**2))
    cov = np.linalg.pinv(fit.jac.T@fit.jac)
    err = np.sqrt(np.maximum(np.diag(cov),0))
    return dict(status="fitted", depth=float(fit.x[3]), depth_diagonal_error=float(err[3]),
                midpoint_btjd=float(epoch+fit.x[4]), midpoint_shift_days=float(fit.x[4]),
                midpoint_diagonal_error_days=float(err[4]), duration_days=float(fit.x[5]),
                duration_diagonal_error_days=float(err[5]),
                chi2=float(np.sum(fit.fun**2)), ndof=len(x)-6,
                fixed_ingress_fraction=.2, exposure_days=exposure,
                boundary_warning=bool(abs(fit.x[4])>.095 or fit.x[5]<.016 or fit.x[5]>.175))


result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
              method="independent LS followed by local polynomial tests on original fluxes",
              method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              versions=dict(python=platform.python_version(), numpy=np.__version__, astropy=astropy.__version__, scipy=scipy.__version__),
              initial_assessment_was_blind_to_parent_BLS=True,
              supplied_hypotheses=dict(period_days=PERIOD, epoch_btjd=EPOCH, duration_days=DURATION,
                                       known_neighbor_period_days=EB_PERIOD, known_neighbor_epoch_btjd_tdb=EB_EPOCH),
              warning="Conditional diagonal errors and Gaussian LS probabilities are diagnostics, not calibrated detection/false-positive probabilities.",
              sectors=[])
rows = []
phasefig, phaseax = plt.subplots(4, 3, figsize=(14, 12))
rawfig, rawax = plt.subplots(4, 2, figsize=(14, 10))
all_t, all_y, all_dy = [], [], []
for ir, record in enumerate(RECORDS):
    with fits.open(record["path"]) as hdul:
        d = hdul[1].data
        t = np.array(d["TIME"], dtype=float)
        q = np.array(d["QUALITY"])
        corrected = "DET_FLUX" if "DET_FLUX" in d.names else "KSPSAP_FLUX"
        names = ["SAP_FLUX", corrected, corrected+"_SML", corrected+"_LAG", "SAP_X", "SAP_Y", "SAP_BKG"]
        channels = {n: np.array(d[n], dtype=float) for n in names}
    m = np.isfinite(t) & (q == 0) & np.isfinite(channels[corrected])
    tt = t[m]
    sector = dict(sector=record["sector"], original=record, n_total=len(t), n_quality_zero=int(m.sum()),
                  time_min=float(tt.min()), time_max=float(tt.max()), channels={}, events=[])
    flux = channels[corrected][m]
    flux = flux / np.median(flux) - 1
    sigma = mad(flux)
    all_t.extend(tt); all_y.extend(flux); all_dy.extend(np.full(len(tt), sigma))
    rawax[ir, 0].plot(tt, channels["SAP_FLUX"][m], ".", ms=1)
    rawax[ir, 1].plot(tt, flux * 1000, ".", ms=1)
    rawax[ir, 0].set_ylabel(f"S{record['sector']} SAP")
    rawax[ir, 1].set_ylabel("Corrected ppt")
    for name, full in channels.items():
        sel = m & np.isfinite(full)
        tx = t[sel]
        y = full[sel]
        isflux = "FLUX" in name
        if isflux:
            y = y / np.median(y) - 1
        ch = dict(median=float(np.median(full[sel])), robust_scatter=mad(y), n=len(y))
        ls = LombScargle(tx, y)
        freq, power = ls.autopower(minimum_frequency=.05, maximum_frequency=20, samples_per_peak=15)
        peaks = find_peaks(power)[0]
        strongest = sorted(peaks, key=lambda k: power[k], reverse=True)[:12]
        ch["initial_unflattened_LS_peaks"] = [dict(period=float(1/freq[k]), power=float(power[k])) for k in strongest]
        ch["known_EB_fundamental"] = fit_sinusoid(tx, y, EB_PERIOD)
        ch["known_EB_half_period"] = fit_sinusoid(tx, y, EB_PERIOD/2)
        ch["five_day_fundamental"] = fit_sinusoid(tx, y, PERIOD)
        if name == corrected:
            sf, sp = ls.autopower(minimum_frequency=6.2, maximum_frequency=6.5, samples_per_peak=100)
            k = sp.argmax()
            ch["best_short_period"] = dict(period=float(1/sf[k]), power=float(sp[k]), fit=fit_sinusoid(tx, y, 1/sf[k]))
            ch["best_short_period"]["gaussian_single_sector_Baluev_FAP_0p05to20cyclesday"] = float(ls.false_alarm_probability(sp[k], minimum_frequency=.05, maximum_frequency=20))
            ch["short_half_samples"] = [fit_sinusoid(tx[sl], y[sl], EB_PERIOD/2) for sl in [tx < np.median(tx), tx >= np.median(tx)]]
        sector["channels"][name] = ch
    for period, epoch, ax in [(EB_PERIOD, EB_EPOCH, phaseax[ir, 0]),
                               (EB_PERIOD/2, EB_EPOCH, phaseax[ir, 1]),
                               (PERIOD, EPOCH, phaseax[ir, 2])]:
        ph = ((tt-epoch)/period + .5) % 1 - .5
        bins = np.linspace(-.5, .5, 61)
        which = np.digitize(ph, bins)-1
        centers = (bins[1:]+bins[:-1])/2
        means = np.array([np.mean(flux[which == j]) if np.any(which == j) else np.nan for j in range(60)])
        errors = np.array([np.std(flux[which == j])/np.sqrt(np.sum(which == j)) if np.sum(which == j)>1 else np.nan for j in range(60)])
        ax.plot(ph, flux*1000, ".", ms=.5, color="gray", alpha=.25)
        ax.errorbar(centers, means*1000, yerr=errors*1000, fmt=".-", color="black", ms=3, lw=.7)
        ax.set_ylim(-25, 25)
        ax.set_title(f"S{record['sector']} P={period:.8f} d")
    first = int(np.floor((tt.min()-EPOCH)/PERIOD))
    last = int(np.ceil((tt.max()-EPOCH)/PERIOD))
    for cycle in range(first, last+1):
        event_time = EPOCH + cycle * PERIOD
        if event_time < tt.min()-.3 or event_time > tt.max()+.3:
            continue
        event = dict(cycle=cycle, predicted_time_btjd=event_time,
                     known_EB_phase=((event_time-EB_EPOCH)/EB_PERIOD) % 1,
                     nearest_quality_zero_days=float(abs(tt-event_time).min()), channels={},
                     inside_quality_counts={str(int(v)):int(np.sum(q[(abs(t-event_time)<DURATION/2)]==v)) for v in np.unique(q[(abs(t-event_time)<DURATION/2)])})
        for name, full in channels.items():
            # All fluxes have their own independent local baseline. Auxiliary
            # variables retain their original units and sign convention.
            e = local_event(t[m], full[m], event_time, relative="FLUX" in name)
            event["channels"][name] = e
        event["corrected_exploratory_shape"] = event_shape(t[m], channels[corrected][m], event_time)
        sector["events"].append(event)
        for name, e in event["channels"].items():
            rows.append(dict(sector=record["sector"], cycle=cycle, predicted_time_btjd=event_time,
                             known_EB_phase=event["known_EB_phase"], channel=name, **e))
        if event["channels"][corrected]["status"] == "measured":
            for ax in rawax[ir]:
                ax.axvline(event_time, color="red", alpha=.5, lw=.6)
    result["sectors"].append(sector)
all_t = np.array(all_t); all_y = np.array(all_y); all_dy = np.array(all_dy)
ls = LombScargle(all_t, all_y, dy=all_dy)
freq, power = ls.autopower(minimum_frequency=6.2, maximum_frequency=6.5, samples_per_peak=15)
peaks = find_peaks(power)[0]
result["combined_short_period_peaks"] = [dict(period=float(1/freq[k]), power=float(power[k])) for k in sorted(peaks, key=lambda k:power[k], reverse=True)[:20]]
summary_rows = []
timings = []
for s in result["sectors"]:
    for name in s["channels"]:
        vals = [e["channels"][name] for e in s["events"] if e["channels"][name]["status"] == "measured"]
        if not vals:
            continue
        weights = np.array([v["error"] for v in vals])**-2
        mean = float(np.sum(weights*np.array([v["depth"] for v in vals]))/weights.sum())
        error = float(1/np.sqrt(weights.sum()))
        summary_rows.append(dict(sector=s["sector"], channel=name, n_events=len(vals),
                                 depth=mean, error=error, conditional_snr=mean/error))
    for event in s["events"]:
        shape = event["corrected_exploratory_shape"]
        if shape["status"] == "fitted":
            timings.append(dict(sector=s["sector"], cycle=event["cycle"], **shape))
# The first S12 event is weak and the free-duration fit broadens to 3.7h.
# It is shown in the table but excluded from this exploratory timing fit.
use_timings = [v for v in timings if not (v["sector"] == 12 and v["cycle"] == 0) and not v["boundary_warning"]]
cycles = np.array([v["cycle"] for v in use_timings],dtype=float)
tm = np.array([v["midpoint_btjd"] for v in use_timings])
te = np.array([v["midpoint_diagonal_error_days"] for v in use_timings])
X = np.c_[np.ones(len(cycles)), cycles]
cov = np.linalg.inv(X.T @ ((1/te**2)[:,None]*X))
beta = cov @ X.T @ (tm/te**2)
chi2 = float(np.sum(((tm-X@beta)/te)**2))
result["exploratory_timing_ephemeris"] = dict(epoch_btjd=float(beta[0]), period_days=float(beta[1]),
        epoch_diagonal_error_days=float(np.sqrt(cov[0,0])), period_diagonal_error_days=float(np.sqrt(cov[1,1])),
        chi2=chi2, ndof=len(tm)-2, n_events=len(tm),
        exclusion="weak first S12 event at cycle0 excluded; all other fitted events retained",
        warning="free-duration per-event timings, diagonal uncertainties; not an independently validated physical ephemeris")
result["local_event_sector_summaries"] = summary_rows
for fig, name in [(rawfig, "original_lightcurves.png"), (phasefig, "independent_phase_folds.png")]:
    fig.tight_layout(); fig.savefig(OUT / name, dpi=150); plt.close(fig)
import pandas as pd
pd.DataFrame(rows).to_csv(OUT / "local_event_measurements.csv", index=False)
pd.DataFrame(summary_rows).to_csv(OUT / "local_event_sector_summaries.csv", index=False)
pd.DataFrame(timings).to_csv(OUT / "exploratory_event_shapes.csv", index=False)
(OUT / "assessment.json").write_text(json.dumps(clean_json(result), indent=2) + "\n")
np.savez_compressed(OUT / "combined_corrected_arrays.npz", time_btjd=all_t, normalized_corrected=all_y,
                    sector_sigma=all_dy, short_frequencies=freq, short_power=power)
print(json.dumps(dict(output=str(OUT), combined_short_peaks=result["combined_short_period_peaks"][:5]), indent=2))
