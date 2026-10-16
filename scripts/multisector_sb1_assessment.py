"""Derive the mass-function implication of a retained published Gaia SB1 row.

This does not fit RV epochs or validate a catalog orbit. The primary mass is
the saved ExoFOP/TIC estimate; sensitivity endpoints are not posterior intervals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from math import pi, sqrt
from pathlib import Path

from astropy.constants import G, M_sun
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]


def sha(body):
    return hashlib.sha256(body).hexdigest()


def assessment(summary, tic, photometric_period):
    summary = Path(summary).resolve()
    raw_summary = summary.read_bytes()
    record = next(r for r in json.loads(raw_summary)["checks"] if r["tic"] == tic)
    nss_path = summary.parent / f"tic{tic}_gaia_nss_orbit.csv"
    raw_nss = nss_path.read_bytes()
    manifest = json.loads(nss_path.with_suffix(".csv.manifest.json").read_text())
    if not manifest["valid"] or sha(raw_nss) != manifest["sha256"]:
        raise ValueError("Invalid or changed retained Gaia orbit response")
    orbit = next(r for r in csv.DictReader(io.StringIO(raw_nss.decode("utf-8-sig")))
                 if r["source_id"] == record["gaia_dr3"] and r["nss_solution_type"] == "SB1")
    p, ep, k, ek, e, ee = [float(orbit[key]) for key in [
        "period", "period_error", "semi_amplitude_primary", "semi_amplitude_primary_error",
        "eccentricity", "eccentricity_error"]]
    f = p * 86400 * (k * 1000)**3 * (1-e*e)**1.5 / (2*pi*G.value*M_sun.value)
    ef = f * sqrt((ep/p)**2 + (3*ek/k)**2 + (3*e*ee/(1-e*e))**2)
    exofop = ROOT / record["exofop_source"]
    raw_exofop = exofop.read_bytes()
    if sha(raw_exofop) != record["exofop_sha256"]:
        raise ValueError("Changed retained primary-star data")
    stellar = next(r for r in json.loads(raw_exofop)["stellar_parameters"] if r.get("mass"))
    mass, mass_error = float(stellar["mass"]), float(stellar["mass_e"])

    def minimum(fm, primary):
        return brentq(lambda secondary: secondary**3/(primary+secondary)**2-fm, 1e-5, 100)

    return dict(tic=tic, gaia_dr3=record["gaia_dr3"],
                source_summary=str(summary.relative_to(ROOT)), source_summary_sha256=sha(raw_summary),
                raw_nss_path=str(nss_path.relative_to(ROOT)), raw_nss_sha256=sha(raw_nss),
                nss_solution=orbit, mass_function_msun=f, mass_function_formal_error_msun=ef,
                primary_mass_msun=mass, primary_mass_error_msun=mass_error,
                primary_mass_provenance="ExoFOP TIC v8.2 2019-04-15",
                exofop_source=str(exofop.relative_to(ROOT)), exofop_sha256=sha(raw_exofop),
                minimum_companion_mass_msun=minimum(f, mass),
                minimum_companion_mass_primary_minus1sigma=minimum(f, mass-mass_error),
                minimum_companion_mass_primary_plus1sigma=minimum(f, mass+mass_error),
                minimum_companion_mass_f_minus1sigma=minimum(f-ef, mass),
                minimum_companion_mass_f_plus1sigma=minimum(f+ef, mass),
                photometric_period_days=photometric_period,
                period_difference_seconds=(photometric_period-p)*86400,
                equations="f(M) = P K1^3 (1-e^2)^(3/2)/(2 pi G) = M2^3 sin(i)^3/(M1+M2)^2; minimum sets sin(i)=1",
                method_sha256=sha(Path(__file__).read_bytes()),
                scope="Derived under the published Gaia SB1 orbital interpretation and TIC primary-mass estimate. Formal mass-function error assumes independent P,K,e errors and excludes their covariance/systematics. Primary-mass endpoints are sensitivity checks, not a posterior confidence interval. This does not independently re-fit Gaia RV epochs, identify the photometric source, or confirm a companion. A matching known spectroscopic-binary period and stellar-scale minimum mass strongly undermine a novel-planet interpretation.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True, help="Immutable deep/<timestamp>/summary.json")
    parser.add_argument("--tic", type=int, required=True)
    parser.add_argument("--photometric-period", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = assessment(args.summary, args.tic, args.photometric_period)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: result[k] for k in ("tic", "mass_function_msun", "minimum_companion_mass_msun", "period_difference_seconds")}))


if __name__ == "__main__":
    main()
