"""Rank public sector-107 SPOC TCEs for independent follow-up, not discovery.

Public TCE metadata are prior signal recognition. Scores are explicit heuristic
priorities, not false-alarm probabilities or planetary odds. No paid service.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import shutil

import numpy as np
import pandas as pd

from broad_crossmatch import CatalogIndex

ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = ("https://archive.stsci.edu/missions/tess/catalogs/tce/"
              "tess2026221140108-s0107-s0107_dvr-tcestats.csv")
EXPECTED_SHA256 = "94a19264c671a0d3b6c1857ba4870d596bd882e05517bd599977645189947382"
DEFAULT_CACHE = Path("/workspace/tess-research-data/broad_probe/tce_s107.csv")


def circular_duration(period, rho_solar, ror):
    """Circular central first-to-fourth contact duration, hours.

    TIC density is relative to solar density. An excess is a ranking flag,
    not rejection: eccentricity and erroneous stellar parameters can change it.
    """
    if not all(np.isfinite(x) and x > 0 for x in (period, rho_solar)):
        return np.nan
    rho_si = rho_solar * 1408.0
    scaled_axis = (6.67430e-11 * rho_si * (period * 86400) ** 2 / (3 * math.pi)) ** (1 / 3)
    return period * 24 / math.pi * math.asin(min(1, (1 + max(0, ror)) / scaled_axis))


def sibling_flags(df):
    """Flag repeated-period TCEs and possible secondary eclipses on one star."""
    same = np.zeros(len(df), dtype=bool)
    secondary = np.zeros(len(df), dtype=bool)
    evidence = [[] for _ in range(len(df))]
    for _, group in df.groupby("tic"):
        rows = list(group.itertuples())
        for a in rows:
            for b in rows:
                if a.Index == b.Index or a.tce_period <= 0:
                    continue
                mismatch = abs(a.tce_period / b.tce_period - 1)
                if mismatch > .01:
                    continue
                same[a.Index] = True
                # Use the full-precision trapezoidal epoch; geometric epochs
                # in this public CSV are rounded at about 0.01 day.
                phase = ((a.tcet_time0bt - b.tcet_time0bt) / b.tcet_period) % 1
                is_secondary = (.35 <= phase <= .65 and
                                b.tce_depth > 1.4 * a.tce_depth and
                                b.tce_depth > 2000)
                secondary[a.Index] |= is_secondary
                evidence[a.Index].append({"other_tce": b.tceid,
                                         "fractional_period_mismatch": float(mismatch),
                                         "phase_offset": float(phase),
                                         "other_depth_ppm": float(b.tce_depth),
                                         "possible_secondary": bool(is_secondary)})
    return same, secondary, [json.dumps(x, separators=(",", ":")) for x in evidence]


def run(source, output, report, limit):
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError("Sector-107 source changed; preserve and investigate before replacing expected hash")
    df = pd.read_csv(source, comment="#")
    if len(df) != 846 or df.ticid.nunique() != 623:
        raise ValueError("Pinned source no longer has the expected 846 TCEs / 623 TICs")
    index = CatalogIndex()
    df["tic"] = df.ticid.astype("int64")
    df["period_days"] = df.tce_period
    df = index.annotate_dataframe(df)
    df["pipeline_mes"] = df.tce_max_mult_ev
    df["pipeline_model_snr"] = df.tce_model_snr
    df["epoch_btjd_followup"] = df.tcet_time0bt
    df["period_days_followup"] = df.tcet_period
    df["prior_signal_status"] = "public_SPOC_sector107_threshold_crossing_2026-10-02"
    df["contamination_ratio"] = np.nan  # This table does not contain CROWDSAP.
    df["odd_even_chi2"] = df.tce_bin_oedp_stat
    df["odd_even_sqrt_chi2"] = np.sqrt(np.maximum(df.odd_even_chi2, 0))
    df["weak_secondary_mes"] = df.tce_ws_maxmes
    df["weak_secondary_depth_ratio"] = df.wst_depth / df.tce_depth.where(df.tce_depth > 0)
    df["tic_centroid_offset_arcsec"] = df.tce_ditco_msky
    df["tic_centroid_offset_snr"] = df.tce_ditco_msky / df.tce_ditco_msky_err.where(df.tce_ditco_msky_err > 0)
    df["centroid_diagnostic_available"] = df.tce_ditco_msky_err > 0
    df["stellar_parameters_placeholder"] = ((df.tce_sradius_prov == "Solar") |
                                              (df.tce_steff_prov == "Solar"))
    df["circular_central_duration_h"] = [circular_duration(p, rho, ror)
                                          for p, rho, ror in zip(df.tce_period, df.tce_sdensity, df.tce_ror)]
    df["duration_to_circular_ratio"] = df.tce_duration / df.circular_central_duration_h
    same, secondary, sibling_evidence = sibling_flags(df)
    df["same_period_sibling_tce"] = same
    df["possible_secondary_sibling_tce"] = secondary
    df["sibling_evidence_json"] = sibling_evidence
    df["planet_size_screen"] = df.tce_prad.between(.3, 20)
    df["screen_eligible"] = ((~df.known_period_match) & (~df.known_harmonic_match) &
        (~df.is_eb_host) & (~df.possible_secondary_sibling_tce) &
        (df.pipeline_mes >= 9) & (df.tce_num_transits >= 3) & df.planet_size_screen &
        df.tce_depth.between(30, 20000) & df.tce_duration.between(.5, 8) &
        df.tce_sradius.between(.1, 3) & (df.tce_steff < 7500) & df.tce_full_conv &
        (df.pipeline_model_snr > 0) &
        (df.tce_period > 0))
    # Explicit priority penalties avoid presenting official DV metrics as a
    # calibrated probability. Severe flags remain in the full searchable table.
    score = np.log1p(df.pipeline_mes)
    score -= np.maximum(df.odd_even_sqrt_chi2 - 3, 0)
    score -= np.maximum(df.weak_secondary_mes - 6, 0) / 4
    score -= np.maximum(df.tic_centroid_offset_snr.fillna(3) - 3, 0)
    score -= np.maximum(df.tce_prad - 12, 0) / 5
    score -= np.maximum(df.duration_to_circular_ratio.fillna(1.5) - 1.5, 0) * 2
    score -= df.stellar_parameters_placeholder.astype(float) * .5
    score -= df.same_period_sibling_tce.astype(float) * .5
    df["heuristic_priority"] = score
    df["followup_class"] = np.where(df.screen_eligible, "unmatched_public_TCE_requires_independent_vetting", "deprioritized_metadata_screen")
    flags = []
    for r in df.itertuples():
        f = []
        for yes, name in [(r.known_period_match, "cached_known_period"),
                          (r.known_harmonic_match, "cached_known_harmonic"),
                          (r.is_eb_host, "explicit_EB_member"),
                          (r.possible_secondary_sibling_tce, "possible_secondary_of_same_star_TCE"),
                          (r.same_period_sibling_tce, "duplicate_period_on_same_star"),
                          (r.odd_even_sqrt_chi2 > 3, "odd_even_DV_statistic"),
                          (r.weak_secondary_mes > 6, "weak_secondary_DV_MES"),
                          (r.tic_centroid_offset_snr > 3 and r.tic_centroid_offset_arcsec > 3, "DV_centroid_offset"),
                          (not r.centroid_diagnostic_available, "centroid_unavailable"),
                          (r.stellar_parameters_placeholder, "Solar_stellar_placeholder"),
                          (r.duration_to_circular_ratio > 1.5, "long_duration_for_circular_central_transit"),
                          (r.tce_model_snr < 9, "geometric_fit_SNR_below_9"),
                          (r.tce_num_transits < 3, "fewer_than_three_reported_transits"),
                          (not r.planet_size_screen, "pipeline_radius_outside_screen")]:
            if yes: f.append(name)
        flags.append(";".join(f))
    df["diagnostic_flags"] = flags
    df = df.sort_values(["screen_eligible", "heuristic_priority", "tic"], ascending=[False, False, True])
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False, float_format="%.15g")
    top = df[df.screen_eligible].head(limit)
    report.parent.mkdir(parents=True, exist_ok=True)
    errcols = [c for c in df.columns if c.startswith("tce_") and c.endswith("_err") and df[c].nunique() == 1]
    columns = ["tic", "tceid", "period_days_followup", "epoch_btjd_followup", "pipeline_mes", "tce_num_transits", "tce_prad", "tce_depth", "tce_duration", "diagnostic_flags"]
    body = "# Public sector-107 TCE triage\n\n"
    body += f"Checked 2026-10-13: **846 previously published SPOC TCEs on 623 unique stars**, not 623 new detections and not an independent 160k-star search. {int(df.screen_eligible.sum())} rows meet the initial metadata screen. This is prioritization for independent time-series retrieval, not a credible candidate list.\n\n"
    body += f"Source: [{SOURCE_URL}]({SOURCE_URL}); created 2026-10-02; SHA-256 `{digest}`. The pinned raw CSV is `data/broad_tce107.csv`. Reproduce with `python src/triage_tce107.py`.\n\n"
    body += "The archive MES (`tce_max_mult_ev`) and geometric-model S/N are distinct. No official bootstrap FAP is treated as a survey-level probability. Crossmatching uses confirmed planets including RV detections, TOI, CTOI, DIAmante, explicit TESS EB membership and cached published candidates via `CatalogIndex`. TIC/period/harmonic matches are conservative ranking flags. This CSV contains no sky coordinates or CROWDSAP; position-based contamination screening and quantitative dilution remain unavailable until follow-up. A nonmatch is not proof of novelty. Literature and alias checks remain pending.\n\n"
    body += "Screen: pipeline MES≥9; at least three reported transits; pipeline radius 0.3–20 Earth radii; depth 30–20,000 ppm; duration 0.5–8 hours; stellar radius 0.1–3 solar; Teff<7500 K; converged geometric fit with positive model S/N (excluding missing/failed −1 sentinels); no cached known-period/harmonic or explicit EB host; no likely secondary of another same-star TCE. The ranking uses log(1+MES), with explicit penalties for odd/even, weak-secondary and centroid statistics, very large radius, long duration, Solar placeholders and duplicate-period TCEs. These choices are transparent heuristics, not probabilities, completeness estimates or physical validation.\n\n"
    body += "**Within-star TCE comparison matters:** TIC 146435300's appealing third row has the same period as a deeper first row and is displaced by approximately half an orbit; it is a likely secondary eclipse. TIC 178288581 similarly has a larger-than-planet first TCE plus a weaker half-phase row. Both would be misleading if only one individually planet-sized row were selected. These flags and all 846 rows remain searchable in the output CSV.\n\n"
    body += "**Export caution:** many geometric-fit uncertainty columns have the same value in all 846 rows, while trapezoidal uncertainties are uniformly −1. Do not use these as trustworthy target-specific errors. Geometric transit epochs are rounded (often to 0.01 day). The table supplies full-precision `tcet_period` / `tcet_time0bt` as follow-up starting points, with the geometric values retained; both require native-FITS refitting.\n\n"
    body += "Constant uncertainty columns: " + ", ".join(f"`{c}`" for c in errcols) + ".\n\n"
    body += "| " + " | ".join(columns) + " |\n|" + "|".join(["---"] * len(columns)) + "|\n"
    for r in top[columns].to_dict("records"):
        body += "| " + " | ".join(str(r[c]) for c in columns) + " |\n"
    body += "\nNext: independently retrieve native SPOC sector-107 light curves for the highest-ranked clean targets; measure events in SAP and PDCSAP with local baselines, assess individual recurrence, refit epoch/period/depth/duration and uncertainties, query prior sectors, then inspect pixels/centroids and catalogs/literature for any survivor. Public TCE recognition must remain explicit even if a planetary interpretation later survives. No new planet is claimed from this table.\n"
    report.write_text(body)
    print(top[["tic", "tceid", "period_days_followup", "epoch_btjd_followup", "pipeline_mes", "tce_num_transits", "diagnostic_flags"]].to_string(index=False))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=ROOT / "data/broad_tce107.csv")
    p.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--output", type=Path, default=ROOT / "results/tce107_triage.csv")
    p.add_argument("--report", type=Path, default=ROOT / "reports/tce107_triage.md")
    p.add_argument("--limit", type=int, default=20)
    args = p.parse_args()
    if not args.source.exists():
        if hashlib.sha256(args.cache.read_bytes()).hexdigest() != EXPECTED_SHA256:
            raise ValueError("Cached source checksum mismatch")
        args.source.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.cache, args.source)
    run(args.source, args.output, args.report, args.limit)


if __name__ == "__main__":
    main()
