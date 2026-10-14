"""Independently falsify three metadata-shortlisted public sector-107 TCEs.

Reads immutable native SPOC observations from the separate follow-up manifest.
Writes only candidates/tce107. Formal scores are not survey probabilities.
"""
from pathlib import Path
import argparse
import json
import math
import os

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle

os.environ.setdefault("MPLCONFIGDIR", "/workspace/tess-research-data/matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from search import ROOT, read_curve, flatten, phase
from validate import local_depth
from broad_crossmatch import CatalogIndex

OUT = ROOT / "candidates/tce107"
TICS = [159923666, 178155030, 100777592]


def safe(v):
    if isinstance(v, dict): return {str(k): safe(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)): return [safe(x) for x in v]
    if isinstance(v, np.generic): return safe(v.item())
    if isinstance(v, float) and not math.isfinite(v): return None
    return v


def aggregate(events):
    if not events: return {"n_events": 0, "events": []}
    dep = np.array([r["depth"] for r in events])
    err = np.array([r["formal_error"] for r in events])
    w = 1 / err ** 2
    mean = np.average(dep, weights=w)
    formal = np.sqrt(1 / w.sum())
    chi = np.sum(w * (dep - mean) ** 2) / max(1, len(dep) - 1)
    scat = formal * np.sqrt(max(1, chi))
    return {"depth": mean, "formal_error": formal, "formal_snr": mean / formal,
            "event_scatter_error": scat, "event_scatter_snr": mean / scat,
            "event_reduced_chi2": chi, "n_events": len(events), "events": events}


def local_polynomial(t, f, e, p, t0, d, order):
    """Separate native local baseline per event; preserve negative depths."""
    events = []
    cyc = np.rint((t - t0) / p).astype(int)
    for cycle in np.unique(cyc):
        x = (t - (t0 + cycle * p)) / d
        inside = np.abs(x) < .5
        outside = (np.abs(x) > 1) & (np.abs(x) < min(3, .4 * p / d))
        use = inside | outside
        if inside.sum() < 5 or (outside & (x < 0)).sum() < 5 or (outside & (x > 0)).sum() < 5:
            continue
        xx = x[use]
        design = np.column_stack([xx ** k for k in range(order + 1)] + [-inside[use].astype(float)])
        w = 1 / e[use] ** 2
        normal = design.T @ (w[:, None] * design)
        if np.linalg.cond(normal) > 1e10: continue
        covariance = np.linalg.inv(normal)
        coef = covariance @ (design.T @ (w * f[use]))
        residual = f[use] - design @ coef
        rc = np.sum(w * residual ** 2) / max(1, use.sum() - len(coef))
        events.append({"cycle": int(cycle), "center_btjd": t0 + cycle * p,
                       "depth": coef[-1], "formal_error": np.sqrt(covariance[-1, -1]),
                       "local_residual_scaled_error": np.sqrt(covariance[-1, -1] * max(1, rc)),
                       "local_reduced_chi2": rc, "n_in": int(inside.sum()),
                       "n_baseline": int(outside.sum())})
    return aggregate(events)


def parity(events):
    odd = aggregate([r for r in events if r["cycle"] % 2])
    even = aggregate([r for r in events if not r["cycle"] % 2])
    if not odd["n_events"] or not even["n_events"]:
        return {"available": False, "reason": "one event parity absent"}
    return {"available": True, "odd_depth": odd["depth"], "even_depth": even["depth"],
            "formal_difference_sigma": abs(odd["depth"] - even["depth"]) / np.hypot(odd["formal_error"], even["formal_error"]),
            "odd_n": odd["n_events"], "even_n": even["n_events"]}


def harmonic_outside(t, f, e, p, t0, d, terms):
    """Fit outside selected dip; inspect narrow residual after periodic wave.

    The harmonic model is phenomenological. A remaining dip is not a planet;
    conversely waveform removal can attenuate real transits.
    """
    angle = 2 * np.pi * (t - t0) / p
    design = np.column_stack([np.ones(len(t))] + [fun(k * angle) for k in range(1, terms + 1) for fun in (np.sin, np.cos)])
    keep = np.abs(phase(t, p, t0)) > max(.7 * d, .06 * p)
    w = 1 / e[keep] ** 2
    coef = np.linalg.lstsq(design[keep] * np.sqrt(w[:, None]), f[keep] * np.sqrt(w), rcond=None)[0]
    model = design @ coef
    residual = f - model + 1
    measured = local_polynomial(t, residual, e, p, t0, d, 1)
    return {"harmonics": terms, "coefficients": coef, "excluded_halfwidth_days": max(.7 * d, .06 * p),
            "raw_variance_explained": 1 - np.var(f - model) / np.var(f),
            "remaining_local_primary": measured}, residual


def phase_bins(t, y, p, t0, n=80):
    ph = phase(t, p, t0) / p
    edges = np.linspace(-.5, .5, n + 1)
    mids, vals = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        use = (ph >= lo) & (ph < hi)
        if use.sum() >= 3:
            mids.append((lo + hi) / 2); vals.append(np.median(y[use]))
    return np.array(mids), np.array(vals)


def vet(tic, rows, lead, index):
    p, t0, d = float(lead.period_days_followup), float(lead.epoch_btjd_followup), float(lead.tce_duration) / 24
    result = {"tic": tic, "public_tce": lead.tceid, "prior_signal_recognition": "public SPOC sector-107 TCE, created 2026-10-02",
              "ephemeris": {"period_days": p, "epoch_btjd": t0, "duration_days": d, "source": "full-precision public trapezoidal P/T0, geometric duration"},
              "observations": rows, "sector_results": [], "scientific_claim": "no credible new planet established"}
    fig, axes = plt.subplots(len(rows), 3, figsize=(14, 3.4 * len(rows)), squeeze=False, layout="constrained")
    for ix, obs in enumerate(sorted(rows, key=lambda r: r["sector"])):
        sector_result = {"sector": obs["sector"], "measurements": []}
        with fits.open(obs["local_path"]) as h:
            sector_result["crowdsap"] = h[1].header.get("CROWDSAP")
            sector_result["flfrcsap"] = h[1].header.get("FLFRCSAP")
            ra, dec = h[0].header.get("RA_OBJ"), h[0].header.get("DEC_OBJ")
        for col in ("PDCSAP_FLUX", "SAP_FLUX"):
            t, f, e, header = read_curve(obs["local_path"], col)
            norm = np.median(f); f = f / norm; e = e / norm
            diff = np.diff(f)[np.diff(t) < .01]
            scatter = 1.4826 * np.median(abs(diff - np.median(diff))) / np.sqrt(2)
            e = np.maximum(e, scatter)
            measurements = {"column": col, "native_good_samples": len(t),
                            "native_point_scatter_ppm": scatter * 1e6,
                            "raw_robust_rms_ppm": 1.4826 * np.median(abs(f - np.median(f))) * 1e6,
                            "fixed_linear": local_polynomial(t, f, e, p, t0, d, 1),
                            "fixed_quadratic": local_polynomial(t, f, e, p, t0, d, 2),
                            "half_phase_secondary_linear": local_polynomial(t, f, e, p, t0 + p / 2, d, 1),
                            "alternative_detrending": []}
            measurements["odd_even"] = parity(measurements["fixed_linear"]["events"])
            for width in (1., 2.):
                tb, fb, eb = flatten(t, f, e, width, bin_minutes=10)
                measurements["alternative_detrending"].append({"median_width_days": width,
                    "fixed_local": local_depth(tb, fb, eb, p, t0, d)})
                if col == "PDCSAP_FLUX" and obs["sector"] == 107 and width == 1:
                    # This targeted follow-up searches a declared narrow period
                    # range, not the blind broad survey; all phases are tested.
                    periods = np.linspace(p * .995, p * 1.005, 601)
                    durations = d * np.array([.55, .8, 1., 1.3, 1.7])
                    b = BoxLeastSquares(tb, fb, eb)
                    bl = b.power(periods, durations, objective="snr", oversample=10)
                    k = int(np.argmax(bl.power))
                    result["targeted_fine_bls"] = {"period_days": bl.period[k], "epoch_btjd": bl.transit_time[k],
                        "duration_days": bl.duration[k], "depth": bl.depth[k], "diagonal_snr": bl.depth_snr[k],
                        "phase_offset_from_public_days": phase(np.array([bl.transit_time[k]]), p, t0)[0],
                        "grid_count": len(periods), "period_range_fraction": .005,
                        "limitation": "selected public-TCE neighborhood; not independent discovery or global FAP"}
            if col == "PDCSAP_FLUX":
                freq, power = LombScargle(t, f, e).autopower(minimum_frequency=1 / 15, maximum_frequency=8, samples_per_peak=4)
                j = int(np.argmax(power))
                measurements["lomb_scargle"] = {"maximum_period_days": 1 / freq[j], "power_standard": power[j],
                                                "trial_range_days": [1 / 8, 15], "interpretation": "periodic waveform ranking, not physical rotation proof"}
                harmonic = []
                for terms in (1, 2):
                    hm, residual = harmonic_outside(t, f, e, p, t0, d, terms)
                    harmonic.append(hm)
                    if terms == 2:
                        x, y = phase_bins(t, residual, p, t0)
                        axes[ix, 2].plot(x, (y - 1) * 1e6, "o-", ms=2, label="2-harmonic residual")
                measurements["outside_dip_harmonic"] = harmonic
                axes[ix, 0].plot(t, (f - 1) * 1e6, ".", ms=.8)
                for event in measurements["fixed_linear"]["events"]:
                    axes[ix, 0].axvline(event["center_btjd"], color="red", lw=.5, alpha=.6)
            x, y = phase_bins(t, f, p, t0)
            axes[ix, 1].plot(x, (y - 1) * 1e6, "o-", ms=2, label=col)
            sector_result["measurements"].append(measurements)
            print(tic, obs["sector"], col, "local-linear", measurements["fixed_linear"].get("depth", 0) * 1e6,
                  "event-scatter-SNR", measurements["fixed_linear"].get("event_scatter_snr"), flush=True)
        axes[ix, 0].set(title=f"Sector {obs['sector']} native PDC", xlabel="BTJD", ylabel="Normalized flux (ppm)")
        axes[ix, 1].set(title=f"Public P={p:.6f} d; CROWDSAP={sector_result['crowdsap']}", xlabel="Orbital phase from public epoch", ylabel="Normalized flux (ppm)")
        axes[ix, 1].legend()
        axes[ix, 2].set(title="Out-of-dip periodic wave removed", xlabel="Phase", ylabel="Residual (ppm)")
        axes[ix, 2].axvspan(-d / p / 2, d / p / 2, color="red", alpha=.1)
        axes[ix, 2].axhline(0, color="black", lw=.5)
        result["sector_results"].append(sector_result)
    result["catalog_check"] = index.match_signal(tic, p, ra, dec)
    result["catalog_check"]["scope"] = "dated cached catalogs plus coordinate neighbors; literature/aliases not exhaustive"
    cpv_path = OUT / "cpv2024_table1.csv"
    if cpv_path.exists():
        import hashlib
        provenance = json.loads((OUT / "cpv2024_table1.csv.manifest.json").read_text())
        if hashlib.sha256(cpv_path.read_bytes()).hexdigest() != provenance["sha256"]:
            raise ValueError("CPV prior table changed; preserve and investigate")
        cpv_table = pd.read_csv(cpv_path, sep="|")
        match = cpv_table[cpv_table.ticid == tic]
        result["prior_literature"] = {"paper": "Bouma et al., Transient Corotating Clumps Around Adolescent Low-Mass Stars From Four Years of TESS",
            "arxiv": "https://arxiv.org/abs/2309.06471", "table_provenance": provenance,
            "exact_tic_and_gaia_match": not match.empty,
            "rows": match.to_dict("records"),
            "classification": "known complex periodic variable" if not match.empty else "not in this specific CPV table; no novelty inference"}
        if not match.empty:
            result["scientific_claim"] = "known complex periodic variable; public TCE waveform is not a new planetary detection"
    result["limitations"] = ["Native diagonal errors and event-scatter scores are not calibrated significance or a search FAP.",
        "Archival predictions use a selected public linear ephemeris; long-baseline period uncertainty was not inferred from untrustworthy CSV errors.",
        "Periodic harmonic removal can absorb real transit structure and is a sensitivity diagnostic, not proof of a false positive.",
        "SAP and PDC share photons; their agreement is not independent confirmation.",
        "No pixel source localization or physical planet model is justified without a surviving stable transit interpretation."]
    fig.savefig(OUT / f"tic{tic}_native_waveform_diagnostics.png", dpi=145); plt.close(fig)
    (OUT / f"tic{tic}_measurements.json").write_text(json.dumps(safe(result), indent=2, allow_nan=False))
    return result


def write_report(results):
    lines = ["# Independent follow-up of public sector-107 TCEs", "",
        "These are already-public threshold-crossing signals, independently examined here for possible unrecognized planetary interpretations. **No credible new planet is established.** Read the per-target waveform plots and measurements; a large formal diagonal score is insufficient.", "",
        "Native QUALITY=0 finite positive SAP/PDC flux is preserved and time references are verified. Separate event baselines use linear and quadratic models, central box depths and symmetric nearby out-of-transit coverage. Uncertainties have a native adjacent-sample scatter floor. Alternative one/two-day gap-aware medians, a targeted fine-period BLS, odd/even and half-phase measurements, Lomb–Scargle and out-of-dip one/two-harmonic models expose variability and instability. Raw data and processing provenance are in `data/tce107_followup_observations.json`.", "",
        "| TIC | Sector | Flux | Local linear depth ppm | Formal S/N | Event-scatter S/N | Events | Local quadratic depth ppm | Strong waveform period d |", "|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        for sector in r["sector_results"]:
            for m in sector["measurements"]:
                a, b = m["fixed_linear"], m["fixed_quadratic"]
                def fmt(v, prec=1): return "unavailable" if v is None else f"{v:.{prec}f}"
                depth = a.get("depth"); qdepth = b.get("depth")
                lines.append(f"|{r['tic']}|{sector['sector']}|{m['column']}|{fmt(depth * 1e6 if depth is not None else None)}|{fmt(a.get('formal_snr'))}|{fmt(a.get('event_scatter_snr'))}|{a['n_events']}|{fmt(qdepth * 1e6 if qdepth is not None else None)}|{fmt(m.get('lomb_scargle', {}).get('maximum_period_days'), 6)}|")
    lines += ["", "TIC 159923666: a large approximately 2.734-day waveform explains most variance; the public approximately 8.07-day interpretation groups different waveform phases. Individual sector-107 depths vary strongly, and the predicted earlier-sector primary events are often brightenings. Stable planetary transits are unsupported.", "",
        "TIC 100777592: a strong approximately 0.431-day waveform and highly inconsistent selected-event depths undermine the approximately 5.85-day periodic interpretation. CROWDSAP is approximately 0.397; PDC depths are scaled relative to raw SAP as expected for strong crowding correction. This does not identify the actual variable source. Pixel localization would be needed to attribute the variability; no planet is promoted.", "",
        "TIC 178155030: the approximately 0.487-day modulation repeats in sectors 105–107, with a large asymmetric full-phase waveform and an additional narrow minimum. The dip remains real after one/two-harmonic removal, at approximately 12,163–12,679 ppm for two harmonics; stellar modulation alone is therefore insufficient to dismiss the localized feature. **The decisive prior check identifies it as a previously published complex periodic variable:** the Bouma et al. author table matches both TIC 178155030 and Gaia DR3 5050451154310013056, lists period 11.67 hours, Tucana–Horologium association, age 45 Myr, radius 0.49 solar, and four earlier observed sectors. The paper is [Transient Corotating Clumps Around Adolescent Low-Mass Stars From Four Years of TESS](https://arxiv.org/abs/2309.06471), initially posted 2023-09-12. This recognized waveform and its narrow dips must not be promoted as a newly recognized planetary signal. The public author table, manuscript, source metadata and byte hashes are preserved locally; the table URL is pinned to a Git commit in its manifest. This does not exclude an unrelated additional planet around the star.", "",
        "The coordinate/catalog check remains narrower than a literature novelty investigation. The source TCE table was already public, so none of these signals can be called an original new detection here. Physical planet fitting, injection completeness and centroid/pixel diagnostics cannot repair a waveform interpretation that fails initial screening; any later revival needs independent reduction and targeted physical/literature tests.", "",
        "Reproduce using `python src/triage_tce107.py` then `python src/vet_tce107.py`. The script does not download, overwrite or remove source FITS. All measurements, event-level values, hashes linked through the follow-up manifest and diagnostic plots are retained. These three outcomes do not exclude planets around their hosts or validate a population-level false-positive rate.", ""]
    (OUT / "ASSESSMENT.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/tce107_followup_observations.json")
    parser.add_argument("--tic", type=int, nargs="*", default=TICS)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(args.manifest.read_text())
    triage = pd.read_csv(ROOT / "results/tce107_triage.csv")
    index = CatalogIndex()
    results = []
    for tic in args.tic:
        rows = [r for r in manifest if r["tic"] == tic]
        lead = triage[(triage.tic == tic) & triage.screen_eligible].iloc[0]
        results.append(vet(tic, rows, lead, index))
    write_report(results)


if __name__ == "__main__":
    main()
