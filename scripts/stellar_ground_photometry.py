#!/usr/bin/env python3
"""Anonymous ASAS-SN retrieval with byte/provenance caches and fixed-ephemeris checks.

This helper reads public data only. It does not send observing requests, use
credentials, or change catalog/source ephemerides. ASAS-SN's documented public
client currently labels its response ``arrow`` but serves Parquet bytes.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import importlib.metadata
import io
import json
import os
import platform
import re
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data/stellar_systems/deep_context/ground_photometry"
OUT = ROOT / "candidates/stellar_systems/deep_context/ground_photometry"
AUXILIARY_CACHE = Path(os.environ.get("TESS_DATA_DIR", "/workspace/tess-research-data")) / "stellar_systems/deep_context/ground_photometry"
AUXILIARY_FILES = {"asassn_public_http_script0.js"}
API = "http://asassn-lb01.ifa.hawaii.edu:9006"
TARGET = {"tic": 273879360, "gaia_dr3": "4692910966756908288",
          "ra_deg": 32.252462825, "dec_deg": -71.323170344}
EPHEMERIS = {"period_days": 7.586891287, "primary_epoch_btjd_tdb": 1327.899918,
             "secondary_phase": 0.69179, "primary_duration_days": 0.1743613901293776,
             "secondary_duration_days": 0.13900312680749385,
             "primary_ingress_fraction": 0.4232233081844376,
             "secondary_ingress_fraction": 0.4000693414938672}
ARCHIVE_ID = "60ccbd79-2dcc-58ba-8a9b-41b1ac9c8b08"
ARCHIVE_URL = "https://asas-sn.osu.edu/photometry/" + ARCHIVE_ID


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def response_path(name, metadata=None):
    """Large ancillary website assets retain a tracked retrieval manifest."""
    if name in AUXILIARY_FILES or (metadata or {}).get("storage") == "external_auxiliary_cache":
        return AUXILIARY_CACHE / name
    return CACHE / name


def verify_auxiliary_cache():
    records = []
    for name in sorted(AUXILIARY_FILES):
        manifest_path = CACHE / (name + ".manifest.json")
        if not manifest_path.exists():
            continue
        metadata = json.loads(manifest_path.read_text())
        path = response_path(name, metadata)
        state = "missing; restore with --restore-auxiliary"
        if path.exists():
            if sha256(path.read_bytes()) != metadata["sha256"]:
                raise ValueError(f"Auxiliary cache checksum mismatch: {path}")
            state = "retained; SHA256 verified"
        records.append({"name": name, "local_path": str(path), "state": state,
                        "sha256": metadata["sha256"], "bytes": metadata["bytes"],
                        "scientific_requirement": "none; website API provenance only"})
    return records


def restore_auxiliary_cache():
    for name in sorted(AUXILIARY_FILES):
        manifest_path = CACHE / (name + ".manifest.json")
        if not manifest_path.exists():
            continue
        metadata = json.loads(manifest_path.read_text())
        fetch(name, metadata["requested_url"], method=metadata.get("method", "GET"),
              params=metadata.get("query_parameters"), payload=metadata.get("json_payload"))
    return verify_auxiliary_cache()


def fetch(name, url, *, method="GET", params=None, payload=None, refresh=False):
    CACHE.mkdir(parents=True, exist_ok=True)
    meta_path = CACHE / (name + ".manifest.json")
    previous_meta = json.loads(meta_path.read_text()) if meta_path.exists() else None
    path = response_path(name, previous_meta)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not refresh and path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        body = path.read_bytes()
        if sha256(body) != meta["sha256"]:
            raise ValueError(f"Cached byte checksum mismatch: {path}")
        return body, meta
    meta = {"retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "method": method, "requested_url": url,
            "query_parameters": params, "json_payload": payload,
            "authentication": "none", "tls_verification": True}
    try:
        response = requests.request(method, url, params=params, json=payload,
                                    timeout=(10, 45))
        body = response.content
        meta.update({"url": response.url, "http_status": response.status_code,
                     "content_type": response.headers.get("Content-Type"),
                     "bytes": len(body), "sha256": sha256(body)})
    except requests.RequestException as exc:
        body = str(exc).encode()
        meta.update({"exception_type": type(exc).__name__, "exception": str(exc),
                     "http_status": None, "bytes": len(body), "sha256": sha256(body)})
    if name in AUXILIARY_FILES and previous_meta:
        # Recover original provenance without replacing its expected hash or
        # silently accepting changed public website assets. Retain mismatches.
        if sha256(body) != previous_meta["sha256"]:
            mismatch = path.with_name(path.name + ".retrieved-" + sha256(body))
            mismatch.write_bytes(body)
            Path(str(mismatch) + ".manifest.json").write_text(json.dumps(meta, indent=2) + "\n")
            raise ValueError(f"Original auxiliary asset was not recovered; response retained: {mismatch}")
        rehydration = dict(meta)
        meta = dict(previous_meta)
        meta["last_rehydration"] = rehydration
    if name in AUXILIARY_FILES:
        meta.update({"storage": "external_auxiliary_cache", "local_path": str(path),
                     "cache_relative_path": "stellar_systems/deep_context/ground_photometry/" + name})
    path.write_bytes(body)
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return body, meta


def require_ok(body, meta):
    if meta["http_status"] != 200:
        raise RuntimeError(f"Request failed: {meta['requested_url']} HTTP {meta['http_status']}")
    return body


def decode_table(body):
    if body[:4] == b"PAR1":
        return pd.read_parquet(io.BytesIO(body))
    raise ValueError("Unsupported public ASAS-SN response; preserved raw bytes")


def retrieve(refresh=False):
    # Retain the website failure separately from the actual documented API.
    fetch("asassn_public_https_homepage.txt", "https://asas-sn.ifa.hawaii.edu/skypatrol/",
          refresh=refresh)
    schema_body, schema_meta = fetch("asassn_schema.json", API + "/get_schema", refresh=refresh)
    schema = json.loads(require_ok(schema_body, schema_meta))
    fetch("asassn_counts.json", API + "/get_counts", refresh=refresh)
    # Keep the exact released client and authoritative wheel digest that establish the API protocol.
    pypi_body, pypi_meta = fetch("pyasassn_pypi_0.6.4.json",
                               "https://pypi.org/pypi/pyasassn/0.6.4/json", refresh=refresh)
    release = json.loads(require_ok(pypi_body, pypi_meta))
    wheel = next(item for item in release["urls"] if item["filename"].endswith(".whl"))
    wheel_body, wheel_meta = fetch(wheel["filename"], wheel["url"], refresh=refresh)
    require_ok(wheel_body, wheel_meta)
    if sha256(wheel_body) != wheel["digests"]["sha256"]:
        raise ValueError("PyPI release SHA256 differs from downloaded client wheel")
    with zipfile.ZipFile(io.BytesIO(wheel_body)) as archive:
        (CACHE / "pyasassn_0.6.4_client_source.txt").write_bytes(archive.read("pyasassn/client.py"))
    cols = ["asas_sn_id", "ra_deg", "dec_deg", "tic_id", "gaia_id", "gaia_mag"]
    available = {x["col_names"] for x in schema["stellar_main"]}
    if not set(cols).issubset(available):
        raise ValueError("Expected documented catalog columns are absent")
    payload = {"tar_ids": [TARGET["tic"]], "catalog": "stellar_main", "id_col": "tic_id",
               "cols": cols, "format": "arrow", "download": True}
    body, meta = fetch("asassn_target_index.parquet", API + "/lookup_targets/catalog_list",
                       method="POST", payload=payload, refresh=refresh)
    index = decode_table(require_ok(body, meta))
    index.to_csv(CACHE / "asassn_target_index.csv", index=False)
    if len(index) != 1:
        raise ValueError(f"Expected one target, found {len(index)}")
    # The live HTTP frontend (cached JavaScript) documents this JSON route.
    # HTTP200 may carry an explicit maintenance error, so inspect the body.
    body, meta = fetch("asassn_api_status.json", API + "/get_api_status", refresh=refresh)
    body, meta = fetch("asassn_target_lightcurve.json",
                       API + f"/get_lightcurve/{int(index.iloc[0]['asas_sn_id'])}", refresh=refresh)
    current = json.loads(require_ok(body, meta))
    if "light_curve" in current:
        curve = current["light_curve"]
        pd.DataFrame(curve["data"], columns=curve["columns"]).to_csv(
            CACHE / "asassn_target_photometry.csv", index=False)
    # The older public photometry database remains independent of live worker maintenance.
    fetch("asassn_old_photometry_search.html", "https://asas-sn.osu.edu/photometry",
          params={"ra": TARGET["ra_deg"], "dec": TARGET["dec_deg"], "radius": 0.2}, refresh=refresh)
    fetch("asassn_old_target_photometry.html", ARCHIVE_URL, refresh=refresh)
    fetch("asassn_old_target_photometry.json", ARCHIVE_URL + ".json", refresh=refresh)
    body, meta = fetch("asassn_old_target_photometry.csv", ARCHIVE_URL + ".csv", refresh=refresh)
    photometry = pd.read_csv(io.BytesIO(require_ok(body, meta)))
    cone_cols = cols + ["gaia_b_mag", "gaia_r_mag", "tmass_id"]
    cone_url = (API + f"/lookup_cone/radius{120 / 3600}_ra{TARGET['ra_deg']}"
                f"_dec{TARGET['dec_deg']}")
    body, meta = fetch("asassn_field_index.parquet", cone_url, method="POST",
                       payload={"catalog": "stellar_main", "cols": cone_cols,
                                "format": "arrow", "download": False}, refresh=refresh)
    if meta["http_status"] == 200:
        decode_table(body).to_csv(CACHE / "asassn_field_index.csv", index=False)
    print(json.dumps({"target_index": index.to_dict("records"),
                      "photometry_rows": len(photometry), "columns": list(photometry.columns)}, indent=2))


def circular_dt(phase, center):
    return ((np.asarray(phase) - center + 0.5) % 1 - 0.5) * EPHEMERIS["period_days"]


def templates(hjd, offset_seconds=0):
    phase = ((np.asarray(hjd) - 2457000 - EPHEMERIS["primary_epoch_btjd_tdb"]
              - offset_seconds / 86400) / EPHEMERIS["period_days"]) % 1
    result = []
    for kind, center in [("primary", 0), ("secondary", EPHEMERIS["secondary_phase"])]:
        distance = abs(circular_dt(phase, center))
        duration = EPHEMERIS[kind + "_duration_days"]
        ingress = EPHEMERIS[kind + "_ingress_fraction"]
        result.append(np.clip((duration / 2 - distance) / (duration * ingress), 0, 1))
    return phase, np.column_stack(result)


def robust_normalize(frame, mode="camera_season"):
    from scipy.stats import median_abs_deviation
    d = frame.copy()
    phase, shape = templates(d.hjd)
    d["phase"] = phase
    d["primary_dt_days"] = circular_dt(phase, 0)
    d["secondary_dt_days"] = circular_dt(phase, EPHEMERIS["secondary_phase"])
    # Exclude both entire predicted eclipse families plus a 0.01-day timing guard.
    oot = ((abs(d.primary_dt_days) > EPHEMERIS["primary_duration_days"] / 2 + 0.01)
           & (abs(d.secondary_dt_days) > EPHEMERIS["secondary_duration_days"] / 2 + 0.01))
    if mode == "camera":
        groups = d.camera.astype(str)
    else:
        season = np.floor((d.hjd - 2456500) / 365.25).astype(int).astype(str)
        groups = d.camera.astype(str) + "_" + season
    d["baseline_group"] = groups
    d["baseline_flux"] = np.nan
    d["baseline_mag"] = np.nan
    d["normalized_noise_floor"] = np.nan
    for label, rows in d.groupby("baseline_group"):
        indices = rows.index[oot.loc[rows.index]]
        if len(indices) < 8:
            indices = d.index[(d.camera == rows.camera.iloc[0]) & oot]
        flux = d.loc[indices, "flux"]
        baseline = float(np.median(flux))
        noise = float(median_abs_deviation(flux / baseline, scale="normal"))
        d.loc[rows.index, "baseline_flux"] = baseline
        d.loc[rows.index, "baseline_mag"] = float(np.median(d.loc[indices, "mag"]))
        d.loc[rows.index, "normalized_noise_floor"] = noise
    if mode == "local_30d":
        for idx, row in d.iterrows():
            indices = d.index[oot & (d.camera == row.camera) & (abs(d.hjd - row.hjd) < 30)]
            if len(indices) >= 8:
                flux = d.loc[indices, "flux"]
                baseline = float(np.median(flux))
                d.loc[idx, "baseline_flux"] = baseline
                d.loc[idx, "baseline_mag"] = float(np.median(d.loc[indices, "mag"]))
                d.loc[idx, "normalized_noise_floor"] = float(
                    median_abs_deviation(flux / baseline, scale="normal"))
    d["normalized_flux"] = d.flux / d.baseline_flux
    d["normalized_flux_err"] = d.flux_err / d.baseline_flux
    d["fit_sigma"] = np.maximum(d.normalized_flux_err, d.normalized_noise_floor)
    d["delta_mag"] = d.mag - d.baseline_mag
    d["out_of_eclipse"] = oot
    d["night_block"] = np.floor(d.hjd).astype(int)
    d["primary_cycle"] = np.rint((d.hjd - 2457000 - EPHEMERIS["primary_epoch_btjd_tdb"])
                                  / EPHEMERIS["period_days"]).astype(int)
    d["secondary_cycle"] = np.rint((d.hjd - 2457000 - EPHEMERIS["primary_epoch_btjd_tdb"])
                                    / EPHEMERIS["period_days"] - EPHEMERIS["secondary_phase"]).astype(int)
    return d


def fit_depths(d, offset_seconds=0, weights_multiplier=None):
    # Camera/season nuisance offsets, fixed TESS shapes, no period/timing/width search.
    _, shape = templates(d.hjd, offset_seconds)
    group_design = pd.get_dummies(d.baseline_group, dtype=float).to_numpy()
    design = np.column_stack([group_design, shape])
    y = 1 - np.asarray(d.normalized_flux)
    weights = 1 / np.asarray(d.fit_sigma) ** 2
    if weights_multiplier is not None:
        weights *= weights_multiplier
    normal = design.T @ (weights[:, None] * design)
    rhs = design.T @ (weights * y)
    inverse = np.linalg.pinv(normal)
    beta = inverse @ rhs
    residual = y - design @ beta
    return {"primary_depth": float(beta[-2]), "secondary_depth": float(beta[-1]),
            "primary_conditional_error": float(np.sqrt(inverse[-2, -2])),
            "secondary_conditional_error": float(np.sqrt(inverse[-1, -1])),
            "chi2": float(np.sum(weights * residual ** 2)),
            "nominal_dof": int(np.count_nonzero(weights) - design.shape[1])}


def analyze(bootstrap_nights=2000):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.time import Time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    csv_path = CACHE / "asassn_old_target_photometry.csv"
    for name in ["asassn_old_target_photometry.csv", "asassn_old_target_photometry.html"]:
        raw_path = CACHE / name
        manifest = json.loads(Path(str(raw_path) + ".manifest.json").read_text())
        if sha256(raw_path.read_bytes()) != manifest["sha256"]:
            raise ValueError(f"Retained source checksum mismatch: {raw_path}")
    original = pd.read_csv(csv_path).sort_values("hjd").reset_index(drop=True)
    source_html = (CACHE / "asassn_old_target_photometry.html").read_text()
    embedded = re.search(r"window\.__lightCurveData\s*=\s*'(.*?)';", source_html, re.S)
    embedded_table = pd.DataFrame(json.loads(html.unescape(embedded.group(1))))
    # Independently verify all exported rows against the complete embedded data.
    keys = list(original.columns)
    a = original.sort_values(keys).reset_index(drop=True)
    b = embedded_table[keys].sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b, check_dtype=False)
    good = (np.isfinite(original.hjd) & np.isfinite(original.flux) & np.isfinite(original.flux_err)
            & (original.flux_err > 0) & (original.flux > 0) & (original.mag_err < 99))
    d = robust_normalize(original[good].copy().reset_index(drop=True))
    fitted = fit_depths(d)
    phase, shape = templates(d.hjd)
    d["in_primary"] = shape[:, 0] > 0
    d["in_secondary"] = shape[:, 1] > 0
    d.to_csv(OUT / "asassn_phased_measurements.csv", index=False)
    events = []
    for kind in ["primary", "secondary"]:
        selected = d[d["in_" + kind]]
        for cycle, rows in selected.groupby(kind + "_cycle"):
            events.append({"kind": kind, "cycle": int(cycle), "n_observations": len(rows),
                           "mean_hjd": float(rows.hjd.mean()),
                           "median_delta_mag": float(rows.delta_mag.median()),
                           "minimum_delta_mag": float(rows.delta_mag.min()),
                           "cameras": "/".join(sorted(rows.camera.unique()))})
    event_table = pd.DataFrame(events)
    event_table.to_csv(OUT / "asassn_event_checks.csv", index=False)

    # Resample complete observing-night blocks, preserving paired cameras and intranight dependence.
    # Baselines are held fixed; separate baseline-mode sensitivity is reported below.
    rng = np.random.default_rng(273879360)
    nights, night_index = np.unique(d.night_block, return_inverse=True)
    boot = []
    for _ in range(bootstrap_nights):
        counts = np.bincount(rng.integers(0, len(nights), len(nights)), minlength=len(nights))
        row_weights = counts[night_index]
        # A bootstrap lacking a sampled eclipse family does not estimate that depth.
        if (np.sum(row_weights * d.in_primary) == 0) or (np.sum(row_weights * d.in_secondary) == 0):
            continue
        result = fit_depths(d, weights_multiplier=row_weights)
        boot.append([result["primary_depth"], result["secondary_depth"]])
    boot = np.asarray(boot)
    np.savez_compressed(OUT / "asassn_night_block_bootstrap.npz", depths=boot, seed=273879360)
    sensitivity = {}
    for mode in ["camera", "camera_season", "local_30d"]:
        sensitivity[mode] = fit_depths(robust_normalize(original[good].copy().reset_index(drop=True), mode))
    camera_results = {str(camera): fit_depths(rows) for camera, rows in d.groupby("camera")}
    timing_results = {str(offset): fit_depths(d, offset_seconds=offset) for offset in [-300, 0, 300]}
    timing_counts = {}
    for offset in [-300, 0, 300]:
        _, shifted_shape = templates(d.hjd, offset_seconds=offset)
        timing_counts[str(offset)] = {"primary": int((shifted_shape[:, 0] > 0).sum()),
                                      "secondary": int((shifted_shape[:, 1] > 0).sum())}
    leave_primary_out = []
    for cycle in sorted(d.loc[d.in_primary, "primary_cycle"].unique()):
        kept = d[~(d.in_primary & (d.primary_cycle == cycle))]
        leave_primary_out.append({"omitted_cycle": int(cycle), **fit_depths(kept)})
    source_coord = SkyCoord(32.25233 * u.deg, -71.32316 * u.deg)
    target_coord = SkyCoord(TARGET["ra_deg"] * u.deg, TARGET["dec_deg"] * u.deg)
    summary = {
        "target": TARGET, "archive_name": "AP36443078", "archive_source_uuid": ARCHIVE_ID,
        "archive_position_deg": [32.25233, -71.32316],
        "separation_from_target_arcsec": float(source_coord.separation(target_coord).arcsec),
        "fixed_tess_ephemeris": EPHEMERIS,
        "retrieved_rows": len(original), "accepted_rows": len(d),
        "independent_csv_vs_embedded_html_rows_identical": True,
        "time_range_hjd": [float(d.hjd.min()), float(d.hjd.max())],
        "calendar_range_approx_utc": [Time(d.hjd.min(), format="jd").iso, Time(d.hjd.max(), format="jd").iso],
        "primary_observations": int(d.in_primary.sum()), "secondary_observations": int(d.in_secondary.sum()),
        "primary_distinct_cycles": int(d.loc[d.in_primary, "primary_cycle"].nunique()),
        "secondary_distinct_cycles": int(d.loc[d.in_secondary, "secondary_cycle"].nunique()),
        "primary_every_point_fainter_than_robust_baseline": bool((d.loc[d.in_primary, "delta_mag"] > 0).all()),
        "primary_delta_mag_range": [float(d.loc[d.in_primary, "delta_mag"].min()),
                                    float(d.loc[d.in_primary, "delta_mag"].max())],
        "fixed_shape_fit_camera_season": fitted,
        "night_block_bootstrap": {"requested_draws": bootstrap_nights, "valid_draws": len(boot),
            "primary_depth_percentiles_2_5_50_97_5": np.percentile(boot[:, 0], [2.5, 50, 97.5]).tolist(),
            "secondary_depth_percentiles_2_5_50_97_5": np.percentile(boot[:, 1], [2.5, 50, 97.5]).tolist(),
            "note": "Conditional fixed-shape, fixed-normalization resampling of full observing nights; baseline sensitivity separate."},
        "baseline_sensitivity": sensitivity, "camera_sensitivity": camera_results,
        "timing_offsets_seconds": timing_results,
        "timing_offset_window_counts": timing_counts,
        "leave_one_primary_event_out": leave_primary_out,
        "archive_metadata_epochs": 540,
        "archive_metadata_count_note": "HTML summary field says540; full embedded data, export CSV and JSON count all report552.",
        "tess_fixed_shape_source": "candidates/stellar_systems/tic273879360/source_summary.json",
        "tess_fixed_shape_source_sha256": sha256((ROOT / "candidates/stellar_systems/tic273879360/source_summary.json").read_bytes()),
        "source_csv_sha256": sha256(csv_path.read_bytes()),
        "script_sha256": sha256(Path(__file__).read_bytes()),
        "assessment": "Public earlier V-band photometry independently confirms recurring primary eclipses at the fixed TESS ephemeris. The weaker secondary depression is supportive but is not a secure independent detection in both cameras.",
        "caveats": [
            "Archived HJD is calculated for field centers; website explicitly advises recalculation for accuracy better than about 200 seconds.",
            "TESS ephemeris is BJD_TDB whereas archived HJD time scale is not specified here; no precise HJD-to-BJD conversion or timing fit is justified.",
            "ASAS-SN aperture photometry does not exclude unresolved nearby blends. The catalog position matches the source, and the bright foreground star is about65.6 arcseconds away.",
            "Fixed TESS morphology in a different optical band gives descriptive diluted template amplitudes, not component radii or physical models.",
            "Archived data contain no per-epoch quality flags; all finite detections are retained, with an observed out-of-eclipse scatter floor in weights.",
            "The API catalog entry and old public lightcurve establish prior observations, not a published EB classification or first discovery.",
            "Separate current prior-recognition audit identifies published Kostov2025 neural-network eclipse-like-feature flags for this source; no entirely unrecognized or first-recognition claim is warranted."
        ]
    }
    (OUT / "asassn_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    colors = {"be": "#2563eb", "bf": "#e11d48"}
    for camera, rows in d.groupby("camera"):
        axes[0, 0].scatter(rows.hjd - 2450000, rows.delta_mag, s=10, alpha=.65, color=colors[camera], label=camera)
        repeated = pd.concat([rows.assign(plot_phase=rows.phase - 1), rows.assign(plot_phase=rows.phase)])
        axes[0, 1].scatter(repeated.plot_phase, repeated.normalized_flux, s=11, alpha=.6, color=colors[camera])
        for ax, kind, center in [(axes[1, 0], "primary", 0), (axes[1, 1], "secondary", EPHEMERIS["secondary_phase"])]:
            hours = circular_dt(rows.phase, center) * 24
            nearby = abs(hours) < 7
            ax.errorbar(hours[nearby], rows.normalized_flux[nearby], rows.fit_sigma[nearby],
                        fmt="o", ms=4, lw=.8, alpha=.8, color=colors[camera], label=camera)
    axes[0, 0].invert_yaxis()
    axes[0, 0].set(xlabel="Archived field-center HJD − 2450000", ylabel="V − robust camera/season baseline (mag)",
                   title="Public ASAS-SN V-band archive: 552 epochs, 2014–2018")
    axes[0, 0].legend(title="Camera")
    axes[0, 1].set(xlim=(-.05, 1), ylim=(.45, 1.25), xlabel="Phase at fixed TESS ephemeris", ylabel="Normalized V flux")
    for ax, kind, center in [(axes[1, 0], "primary", 0), (axes[1, 1], "secondary", EPHEMERIS["secondary_phase"])]:
        hours = np.linspace(-7, 7, 600)
        curve_time = 2457000 + EPHEMERIS["primary_epoch_btjd_tdb"] + center * EPHEMERIS["period_days"] + hours / 24
        _, curve_templates = templates(curve_time)
        k = 0 if kind == "primary" else 1
        ax.plot(hours, 1 - fitted[kind + "_depth"] * curve_templates[:, k], color="black", lw=1.6, label="Fixed TESS shape")
        ax.axvspan(-EPHEMERIS[kind + "_duration_days"] * 12,
                   EPHEMERIS[kind + "_duration_days"] * 12, color="grey", alpha=.12)
        ax.axhline(1, color="grey", lw=.8)
        ax.set(xlim=(-7, 7), ylim=(.45, 1.2), xlabel="Hours from predicted " + kind,
               ylabel="Normalized V flux", title=f"{kind.title()}: {int(d['in_' + kind].sum())} points in window")
    axes[1, 0].legend(fontsize=8)
    fig.suptitle("TIC 273879360 / Gaia DR3 4692910966756908288 — independent ground photometry", fontsize=13)
    fig.savefig(OUT / "asassn_fixed_ephemeris.png", dpi=170)
    plt.close(fig)
    print(json.dumps({"primary_points": summary["primary_observations"],
                      "primary_events": summary["primary_distinct_cycles"],
                      "fixed_fit": fitted, "bootstrap": summary["night_block_bootstrap"]}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retrieve", action="store_true")
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--restore-auxiliary", action="store_true",
                        help="Restore optional website JS asset outside Git using its immutable SHA256 manifest")
    parser.add_argument("--bootstrap-nights", type=int, default=2000)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.restore_auxiliary:
        restore_auxiliary_cache()
    context = {"target": TARGET, "fixed_tess_ephemeris": EPHEMERIS,
               "python": sys.version, "platform": platform.platform(),
               "package_versions": {name: importlib.metadata.version(name) for name in
                                    ["requests", "numpy", "pandas", "scipy", "astropy", "matplotlib", "pyarrow"]},
               "auxiliary_cache": verify_auxiliary_cache(),
               "note": "Anonymous public catalog and photometry retrieval; website failure does not imply API failure. Optional website JS is retained outside Git and not needed for scientific analysis."}
    (OUT / "retrieval_context.json").write_text(json.dumps(context, indent=2) + "\n")
    if args.retrieve:
        retrieve(args.refresh)
    if (CACHE / "asassn_old_target_photometry.csv").exists():
        analyze(args.bootstrap_nights)


if __name__ == "__main__":
    main()
