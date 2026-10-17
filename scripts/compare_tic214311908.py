"""Verify the pinned external report and reproduce bounded comparison checks.

No external Python source is executed. Published posterior summaries are
checked as arithmetic, not independently calibrated false-positive rates.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import math
import urllib.request

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/tic214311908_external_comparison"
SOURCE = DATA / "source_snapshot/tess-search/tic214311908"
OUT = ROOT / "candidates/tic214311908_external_comparison"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_sources(restore=False):
    manifest = json.loads((DATA / "source_manifest.json").read_text())
    for row in manifest["files"]:
        path = (ROOT / row["local_path"]).resolve()
        if not path.is_relative_to(DATA.resolve()):
            raise ValueError("Source destination outside the comparison directory")
        if not path.exists() and restore:
            with urllib.request.urlopen(row["url"], timeout=60) as response:
                body = response.read()
            temporary = path.with_suffix(path.suffix + ".retrieval_attempt")
            temporary.parent.mkdir(parents=True, exist_ok=True)
            temporary.write_bytes(body)
            if hashlib.sha256(body).hexdigest() != row["sha256"]:
                raise ValueError("Unexpected public bytes; retained retrieval attempt")
            temporary.replace(path)
        body = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(body)).encode() + b"\0" + body).hexdigest()
        if digest(path) != row["sha256"] or blob != row["git_blob_sha1"]:
            raise ValueError("Pinned source mismatch: " + row["path"])
    return manifest


def wide_pair():
    """Read every VizieR table, retaining integer identifiers as strings."""
    path = ROOT / "data/tic214311908/catalogs/elbadry_binary_cone.tsv"
    table, header, rows = None, None, []
    for line in path.read_text().splitlines():
        if line.startswith("#Name:"):
            table, header = line.partition(":")[2].strip(), None
        elif line.startswith("Source1\tSource2\t"):
            header = line.split("\t")
        elif header and line and not line.startswith("#"):
            fields = [s.strip() for s in line.split("\t")]
            if len(fields) == len(header) and fields[0].isdigit() and fields[1].isdigit():
                row = dict(zip(header, fields))
                if {row["Source1"], row["Source2"]} == {"5925373272850784768", "5925373272850782208"}:
                    rows.append(dict(table=table, source1=row["Source1"], source2=row["Source2"],
                                     projected_separation_AU=float(row["sepAU"]),
                                     chance_alignment_statistic_R=float(row["R"]), binary_type=row["BinType"]))
    if len(rows) != 1 or rows[0]["table"].endswith("/shift"):
        raise ValueError("Expected one target pair in the real, not shifted, catalogue")
    return dict(source=str(path.relative_to(ROOT)), source_sha256=digest(path), rows=rows,
                correction="Original report read the empty shifted-control table and incorrectly reported no El-Badry match. A real main-catalogue pair was already in the preserved response.",
                limitation="The catalogue R statistic is not a new measured orbit or a short-period host identification.")


def main(restore=False):
    manifest = verify_sources(restore)
    external = json.loads((SOURCE / "ephemeris_adopted.json").read_text())["adopted_planning"]
    independent = json.loads((ROOT / "candidates/tic214311908/independent_variability/ephemeris_frozen_pre_holdout.json").read_text())
    period = independent["parameters"][1]
    period_error = math.sqrt(independent["conditional_covariance_days_squared"][1][1])
    draws_path = SOURCE / "analysis/panel_fpp/out/p05_draws.csv"
    with draws_path.open() as stream:
        rows = list(csv.DictReader(stream))
    arrays = {name: np.array([float(r[name]) for r in rows]) for name in ["FPP", "FPP_strict", "NFPP"]}
    published = json.loads((SOURCE / "analysis/panel_fpp/out/p05_summary.json").read_text())
    computed = dict(n=len(rows), FPP_pct_5_16_50_84_95=np.percentile(arrays["FPP"], [5, 16, 50, 84, 95]).tolist(),
                    frac_FPP_lt_0p015=float(np.mean(arrays["FPP"] < .015)),
                    FPP_strict_pct_5_50_95=np.percentile(arrays["FPP_strict"], [5, 50, 95]).tolist(),
                    frac_strict_lt_0p015=float(np.mean(arrays["FPP_strict"] < .015)))
    for key, value in computed.items():
        if not np.allclose(value, published[key], rtol=1e-10, atol=1e-12):
            raise ValueError("Published prior-draw summary differs: " + key)
    flux_ratio = 10 ** (-.4 * (16.8469 - 13.0135))
    required_depth = .0158 / flux_ratio
    minimum_k = math.sqrt(required_depth)
    rho_ratio = 1.8 / 9.445379
    g = rho_ratio ** (1 / 3)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), method_sha256=digest(Path(__file__)),
                  source_manifest_sha256=digest(DATA / "source_manifest.json"), source_commit=manifest["commit"],
                  pinned_files_verified=len(manifest["files"]),
                  period_comparison=dict(independent_days=period, external_days=external["P_d"],
                                         difference_seconds=(period - external["P_d"]) * 86400,
                                         difference_over_quadrature_errors=abs(period - external["P_d"]) / math.hypot(period_error, external["P_err_d"]),
                                         limitation="Shared observations; descriptive agreement, not an independent combined measurement."),
                  ordinary_dark_occultor_geometry=dict(nominal_neighbor_to_target_T_flux=flux_ratio,
                      required_neighbor_eclipse_fraction=required_depth, minimum_radius_ratio_uniform_brightness=minimum_k,
                      maximum_flat_to_total_duration_ratio_central_small_angle=(1 - minimum_k) / (1 + minimum_k),
                      limitation="A uniform-surface-brightness, circular/small-angle geometry check for an ordinary dark occultor. Limb darkening changes the bounds. Luminous compact-star occultations need separate SED and binary modelling; this is not a universal exclusion."),
                  photoeccentric_sanity=dict(density_relation="rho_circular = rho_star * g**3", g="(1+e*sin(omega))/sqrt(1-e**2)",
                                           example_rho_circular_g_cm3=1.8, rho_star_g_cm3=9.445379,
                                           minimum_apoastron_e_ignoring_companion_mass=(1-g*g)/(1+g*g)),
                  zero_exceedance_one_sided_95pct_upper={str(n):1-.05**(1/n) for n in [25, 400, 1000]},
                  published_prior_draw_arithmetic=computed,
                  wide_pair_catalogue_correction=wide_pair(),
                  scope="Input identities, geometry, catalogue rereading and published-table arithmetic. No complete external raw-cube PRF extraction, evidence integration or independent FPP calibration.")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "arithmetic_checks.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restore-sources", action="store_true", help="Download only missing pinned public source files, checking expected hashes; never execute them.")
    main(parser.parse_args().restore_sources)
