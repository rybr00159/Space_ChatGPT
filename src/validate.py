"""Real-data controls and pre-detrending injections; no discovery/FAP claims.

Run with ``python src/validate.py``.  Inputs are immutable cached SPOC FITS
products listed in data/observations.json.  No network requests are made.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import json

import numpy as np

from search import ROOT, flatten, known_for, mask_known, phase, read_curve, scan


def local_depth(t, f, e, period, epoch, duration):
    """Separate weighted linear baseline and box depth for each observed event.

    Points between half and one full duration from center are omitted to limit
    sensitivity to imperfect catalog durations.  Diagonal covariance is formal;
    event-to-event scatter is also returned and must not be called a FAP.
    """
    events = []
    cycles = np.round((t - epoch) / period).astype(int)
    for cycle in np.unique(cycles):
        center = epoch + cycle * period
        x = t - center
        inside = np.abs(x) < duration / 2
        outside = (np.abs(x) > duration) & (np.abs(x) < 3 * duration)
        use = inside | outside
        if inside.sum() < 2 or outside.sum() < 6:
            continue
        xx = x[use] / duration
        design = np.column_stack([np.ones(use.sum()), xx, -inside[use].astype(float)])
        weight = 1 / e[use] ** 2
        normal = design.T @ (weight[:, None] * design)
        if np.linalg.cond(normal) > 1e10:
            continue
        covariance = np.linalg.inv(normal)
        fit = covariance @ (design.T @ (weight * f[use]))
        residual = f[use] - design @ fit
        events.append({"cycle": int(cycle), "center_btjd": float(center),
                       "n_in": int(inside.sum()), "n_baseline": int(outside.sum()),
                       "depth": float(fit[2]), "formal_error": float(np.sqrt(covariance[2, 2])),
                       "reduced_chi2": float(np.sum(weight * residual**2) / (use.sum() - 3))})
    if not events:
        return {"n_events": 0, "events": []}
    depth = np.array([event["depth"] for event in events])
    errors = np.array([event["formal_error"] for event in events])
    weight = 1 / errors ** 2
    mean = np.average(depth, weights=weight)
    formal = np.sqrt(1 / weight.sum())
    # Do not shrink formal errors when residual/event scatter happens to be low.
    event_chi2 = np.sum(weight * (depth - mean)**2) / max(1, len(events) - 1)
    inflated = formal * np.sqrt(max(1, event_chi2))
    return {"depth": float(mean), "formal_error": float(formal),
            "formal_snr": float(mean / formal), "event_scatter_error": float(inflated),
            "event_scatter_snr": float(mean / inflated), "event_reduced_chi2": float(event_chi2),
            "n_events": len(events), "events": events}


def noise_beta(t, f, e, minutes=90):
    """Robust binned scatter / independent-error expectation, a red-noise warning.

    This is a descriptive ratio, dependent on bin alignment and detrending.
    It is not a correction proven to yield calibrated significance.
    """
    group = np.floor((t - t.min()) / (minutes / 1440)).astype(int)
    means, errors = [], []
    for block in np.unique(group):
        use = group == block
        if use.sum() < 3:
            continue
        weight = 1 / e[use] ** 2
        means.append(np.average(f[use], weights=weight))
        errors.append(np.sqrt(1 / weight.sum()))
    means = np.array(means)
    observed = 1.4826 * np.median(np.abs(means - np.median(means)))
    expected = np.median(errors)
    return {"minutes": minutes, "n_blocks": len(means),
            "robust_scatter_ppm": float(observed * 1e6),
            "diagonal_expectation_ppm": float(expected * 1e6),
            "beta": float(observed / expected)}


def choose_injection_epoch(t, known, period, duration):
    """Choose coverage, not a noise minimum; maximize unmasked injected points."""
    retained = mask_known(t, known)
    trial_epochs = t.min() + np.linspace(.2, period - .2, 100)
    coverages = []
    for epoch in trial_epochs:
        inside = np.abs(phase(t, period, epoch)) < duration / 2
        coverages.append((inside & retained).sum() / max(1, inside.sum()))
    return float(trial_epochs[int(np.argmax(coverages))]), float(max(coverages))


def controls(records, widths, phase_trials=48):
    out = []
    rng = np.random.default_rng(260101)
    for record in records:
        t, f, e, _ = read_curve(record["local_path"])
        known = known_for(record["tic"])
        for width in widths:
            tb, fb, eb = flatten(t, f, e, width)
            all_mask = mask_known(tb, known)
            noise = [noise_beta(tb[all_mask], fb[all_mask], eb[all_mask], minutes)
                     for minutes in (60, 90, 120)]
            for planet in known:
                other = [k for k in known if k["name"] != planet["name"]]
                keep = mask_known(tb, other)
                period, epoch = planet["period"], planet["epoch_bjd"] - 2457000
                duration = planet["duration_hours"] / 24
                measurement = local_depth(tb[keep], fb[keep], eb[keep], period, epoch, duration)
                null_depths, null_snrs = [], []
                for offset in rng.uniform(.15 * period, .85 * period, phase_trials):
                    null = local_depth(tb[all_mask], fb[all_mask], eb[all_mask],
                                       period, epoch + offset, duration)
                    if null["n_events"]:
                        null_depths.append(null["depth"])
                        null_snrs.append(null["formal_snr"])
                out.append({"tic": record["tic"], "sector": record["sector"], "width_days": width,
                            "planet": planet["name"], "period_days": period, "duration_days": duration,
                            "ephemeris_source": planet["source"], "measurement": measurement,
                            "phase_controls": {"n": len(null_snrs), "formal_snr_std": float(np.std(null_snrs, ddof=1)),
                                               "max_absolute_formal_snr": float(np.max(np.abs(null_snrs))),
                                               "robust_depth_scatter_ppm": float(1.4826 * np.median(np.abs(null_depths - np.median(null_depths))) * 1e6)},
                            "noise": noise, "mask_retained_fraction": float(all_mask.mean()),
                            "data_uri": record["data_uri"]})
                print("CONTROL", record["sector"], width, planet["name"],
                      measurement.get("depth", 0) * 1e6, measurement.get("formal_snr"), flush=True)
    return out


def injections(record, widths, depths=(0, 500, 1000, 2000), period=4.731, duration=1.8/24):
    t, f, e, _ = read_curve(record["local_path"])
    known = known_for(record["tic"])
    epoch, coverage = choose_injection_epoch(t, known, period, duration)
    transit = np.abs(phase(t, period, epoch)) < duration / 2
    rows = []
    for width in widths:
        for depth_ppm in depths:
            injected = f * (1 - depth_ppm * 1e-6 * transit)
            tb, fb, eb = flatten(t, injected, e, width)
            keep = mask_known(tb, known)
            # Scan the entire production range without passing injected ephemeris.
            peaks, _ = scan(tb[keep], fb[keep], eb[keep], n_peaks=5)
            measured = local_depth(tb[keep], fb[keep], eb[keep], period, epoch, duration)
            matching = [i for i, peak in enumerate(peaks)
                        if abs(peak["period"] / period - 1) < .01
                        and abs(float(phase(np.array([peak["t0_btjd"]]), period, epoch)[0])) < duration / 2]
            row = {"tic": record["tic"], "sector": record["sector"], "width_days": width,
                   "injected_depth_ppm": depth_ppm, "period_days": period, "epoch_btjd": epoch,
                   "duration_days": duration, "pre_detrending_injection": True,
                   "unmasked_raw_transit_fraction": coverage,
                   "recovered_in_top5": bool(matching) if depth_ppm else None,
                   "matching_rank": matching[0] + 1 if matching else None,
                   "fixed_ephemeris_measurement": measured, "blind_peaks": peaks,
                   "data_uri": record["data_uri"]}
            rows.append(row)
            print("INJECTION", width, depth_ppm, "rank", row["matching_rank"],
                  "best", peaks[0]["period"], peaks[0]["snr"], flush=True)
    # Differential depth relative to original noise realization isolates filtering
    # attenuation, but is neither a general completeness estimate nor a transit fit.
    for row in rows:
        zero = next(r for r in rows if r["width_days"] == row["width_days"] and r["injected_depth_ppm"] == 0)
        delta = row["fixed_ephemeris_measurement"]["depth"] - zero["fixed_ephemeris_measurement"]["depth"]
        row["differential_recovered_depth_ppm"] = float(delta * 1e6)
        row["depth_retention_fraction"] = float(delta * 1e6 / row["injected_depth_ppm"]) if row["injected_depth_ppm"] else None
    return rows


def write_report(result):
    lines = ["# Pipeline validation on actual TESS data", "",
             "These tests assess basic detection and filtering behavior; they do not establish novelty or calibrated false-alarm probabilities.", "",
             "## Catalog ephemeris controls", "",
             "All known transiting L 98-59 planets are measured at fixed catalog ephemerides, with other known planets excluded. Each event uses a separate weighted linear local baseline. Errors assume diagonal noise; event scatter and off-transit phase controls expose limitations.", "",
             "|Sector|Filter (d)|Planet|Depth (ppm)|Formal S/N|Event scatter S/N|Events|Phase-control S/N std|",
             "|---|---|---|---|---|---|---|---|"]
    for row in result["controls"]:
        m, p = row["measurement"], row["phase_controls"]
        lines.append(f"|{row['sector']}|{row['width_days']}|{row['planet']}|{m['depth']*1e6:.0f}|{m['formal_snr']:.1f}|{m['event_scatter_snr']:.1f}|{m['n_events']}|{p['formal_snr_std']:.2f}|")
    lines += ["", "## Pre-detrending injections and unchanged-data control", "",
              "Box-shaped signals are injected multiplicatively into quality-selected raw PDCSAP flux before median filtering and 10-minute binning. Errors are preserved. A fixed epoch is chosen using coverage only, then the full 0.7–13-day production BLS search runs without the true ephemeris. This favorable epoch has 100% of injected raw cadences outside known-planet masks; it is not a random-phase completeness test. Matching requires period within 1% and phase within half the injected duration. Zero-depth rows are unchanged real residual data, not simulations from a validated noise null.", "",
              "|Filter (d)|Injected (ppm)|Top-five rank|Measured added depth (ppm)|Depth retained|Best blind P (d)|Best formal S/N|",
              "|---|---|---|---|---|---|---|"]
    for row in result["injections"]:
        peak = row["blind_peaks"][0]
        retention = row["depth_retention_fraction"]
        retained_text = f"{retention:.3f}" if retention is not None else "—"
        lines.append(f"|{row['width_days']}|{row['injected_depth_ppm']}|{row['matching_rank'] or '—'}|"
                     f"{row['differential_recovered_depth_ppm']:.0f}|{retained_text}|"
                     f"{peak['period']:.5f}|{peak['snr']:.2f}|")
    lines += ["", "## Interpretation and limitations", "",
              "- The shallow known planet is a more meaningful end-to-end sanity check than only retrieving the deepest peak. Highest peaks include harmonics; top-three output alone cannot recover every known planet.",
              "- Injections use one star, one sector, one period, one epoch, and ideal boxes. They cannot quantify pilot-wide completeness, detectability of longer periods, ingress/egress effects, or global false-positive frequency. There are no repeated noise realizations.",
              "- Broad known-planet masks preferentially remove certain trial phases and reduce observing duty cycle. They should be expanded using ephemeris uncertainties for recent sectors, and verified against predicted events. Flattening occurs before masking, so known dips can affect neighboring trends.",
              "- Median-filter width, incomplete events near data gaps, and 10-minute binning can alter depths and durations. Event scatter and 60/90/120-minute noise ratios in results/validation.json are descriptive diagnostics, not calibrated uncertainty models.",
              "- BLS S/N and robust SDE are ranking scores, not Gaussian significance or discovery probabilities. Multiple stars, periods, epochs, durations and detrending choices create a substantial trial factor. Off-phase tests here cannot calibrate the maximum over that search.",
              "- The shortest searched box is 1.08 hours; genuine shorter events such as L 98-59 d (catalog 0.84 h) are diluted. Single/double transits and periods beyond half-sector baseline are outside the blind periodic screen.",
              "- A positive in-transit count is not evidence of repeated events: per-event depths, independent sectors, raw/SAP comparisons, spacecraft artifacts, source localization and novelty checks remain essential.", ""]
    (ROOT / "reports" / "validation.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tic", type=int, default=307210830)
    parser.add_argument("--width", type=float, nargs="+", default=[1., 2.])
    parser.add_argument("--phase-trials", type=int, default=48)
    args = parser.parse_args()
    records = [r for r in json.loads((ROOT / "data" / "observations.json").read_text()) if r["tic"] == args.tic]
    result = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "seed": 260101, "scope": "Controls and eight blind scans; no calibrated FAP",
              "controls": controls(records, args.width, args.phase_trials),
              "injections": injections(records[0], args.width)}
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "reports").mkdir(exist_ok=True)
    (ROOT / "results" / "validation.json").write_text(json.dumps(result, indent=2, allow_nan=False))
    write_report(result)


if __name__ == "__main__":
    main()
