"""Independent stellar properties and bounded recognition checks of two leads.

Public response bytes are immutable and hash-checked. No component mass or
binary classification is inferred from Gaia's assumed-binary MSC fits.
Run --download to fetch missing snapshots; the default only derives retained
data, verifies provenance and validates positive catalog-query controls.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import math
import sys

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import stellar_prior_audit as prior

OUT = ROOT / "data/stellar_systems/deep_context"
RESULT = ROOT / "candidates/stellar_systems/deep_context"
IDS = ["4692910966756908288", "5496090260794096128",
       "5496090260793489024", "5496090256498964096"]
SPECTRO = {"galah3": "J/MNRAS/506/150",
           "galah4": "J/other/PASA/42.51",
           "apogee17": "III/286", "rave6": "III/283"}
SPECTRO_TABLE = {"galah3": "J/MNRAS/506/150/stars",
                "galah4": "J/other/PASA/42.51/allstar",
                "apogee17": "III/286/catalog", "rave6": "III/283/ravedr6"}
TARGETS = {273879360: (32.252462825, -71.323170344),
           260352037: (94.787910153, -57.078218222)}


def utc():
    return datetime.now(timezone.utc).isoformat()


def csv_rows(path):
    """A failed TAP response is not a zero-row catalog result."""
    text = path.read_text()
    if text.lstrip().startswith("<"):
        raise ValueError("Failed TAP response: " + str(path))
    reader = csv.DictReader(io.StringIO(text))
    if not reader.fieldnames:
        raise ValueError("Missing table header: " + str(path))
    return list(reader)


def number(row, key):
    value = row.get(key)
    return None if value in (None, "") else float(value)


def download():
    """Retained fetch implementation writes URLs, UTC, SHA and HTTP status."""
    prior.OUT = OUT
    ids = ",".join(IDS)
    jobs = []
    for table in ["astrophysical_parameters", "astrophysical_parameters_supp",
                  "xp_summary"]:
        name = "gaia_" + table
        query = f"select * from gaiadr3.{table} where source_id in ({ids})"
        jobs.append((prior.tap, (name, query, "gaia")))
    # The historical xp_summary probe used a shorter exact-source list; this
    # reproducible query uses that same list to preserve its immutable bytes.
    jobs[-1] = (prior.tap, ("gaia_xp_summary",
        "select * from gaiadr3.xp_summary where source_id in "
        "(4692910966756908288,5496090260794096128,5496090260793489024)", "gaia"))
    jobs += [
        (prior.tap, ("wide_catalog_targets", f'select * from "J/MNRAS/506/2269/catalog" where Source1 in ({ids}) or Source2 in ({ids})', "vizier")),
        (prior.tap, ("wide_catalog_positive_recovery", 'select * from "J/MNRAS/506/2269/catalog" where Source1 = 3770759051395978240 and Source2 = 3770759051395977472', "vizier")),
        (prior.fetch, ("wide_catalog_readme.txt", "https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/506/2269/ReadMe")),
        (prior.fetch, ("gaia_msc_documentation.html", "https://gea.esac.esa.int/archive/documentation/GDR3/Data_analysis/chap_cu8par/sec_cu8par_apsis/ssec_cu8par_apsis_msc.html")),
        (prior.fetch, ("gaia_source_documentation.html", "https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html")),
        (prior.fetch, ("gaia_astrophysical_parameters_documentation.html", "https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_astrophysical_parameter_tables/ssec_dm_astrophysical_parameters.html")),
        (prior.fetch, ("gaia_astrophysical_parameters_supp_documentation.html", "https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_astrophysical_parameter_tables/ssec_dm_astrophysical_parameters_supp.html")),
    ]
    for label, table in SPECTRO_TABLE.items():
        suffix = "_with_coordinates" if label == "apogee17" else ""
        where = " where RAJ2000 is not null and DEJ2000 is not null" if suffix else ""
        jobs.append((prior.tap, (label + "_positive_sample" + suffix,
                                f'select top 1 * from "{table}"' + where, "vizier")))
    for tic, (ra, dec) in TARGETS.items():
        for label, source in {**SPECTRO, "tmass": "II/246/out", "allwise": "II/328/allwise"}.items():
            params = {"-source": source, "-c": f"{ra} {dec}",
                      "-c.rs": 30 if label in ["tmass", "allwise"] else 10,
                      "-out.max": "unlimited", "-out.all": ""}
            jobs.append((prior.fetch, (f"tic{tic}_{label}.tsv",
                         "https://vizier.cds.unistra.fr/viz-bin/asu-tsv", params)))
    with ThreadPoolExecutor(max_workers=5) as pool:
        for future in [pool.submit(fn, *args) for fn, args in jobs]:
            future.result()
    # Validate the *same* cone endpoint and parameters at an actual positive
    # location from each catalog, not merely the existence of its metadata.
    for label, source in SPECTRO.items():
        suffix = "_with_coordinates" if label == "apogee17" else ""
        rows = csv_rows(OUT / (label + "_positive_sample" + suffix + ".csv"))
        row = rows[0]
        ra = row.get("RA_ICRS") or row.get("RAJ2000")
        dec = row.get("DE_ICRS") or row.get("DEJ2000")
        if not ra or not dec:
            raise ValueError("Positive-control coordinate missing")
        # ASU does not interpret the scientific-notation RA returned by TAP
        # for the APOGEE sample. Preserve that failed cone and use decimals.
        name = label + "_positive_cone_retry1.tsv" if label == "apogee17" else label + "_positive_cone.tsv"
        prior.fetch(name,
                    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
                    {"-source": source, "-c": f"{float(ra):.8f} {float(dec):.8f}", "-c.rs": 10,
                     "-out.max": "unlimited", "-out.all": ""})


def verify():
    records = []
    paths = sorted(OUT.glob("*.manifest.json"))
    paths += [ROOT / "data/stellar_systems/prior_audit/supplement_gaia_sources.csv.manifest.json",
              ROOT / "data/stellar_systems/prior_audit/supplement_tic_sources.json.manifest.json",
              ROOT / "data/stellar_systems/prior_audit/tic273879360_10k2025.tsv.manifest.json",
              ROOT / "data/stellar_systems/scout_novelty/tess10000_paper.html.manifest.json"]
    for path in paths:
        source = Path(str(path)[:-len(".manifest.json")])
        meta = json.loads(path.read_text())
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        if digest != meta["sha256"]:
            raise RuntimeError("Source hash mismatch: " + str(source))
        records.append(dict(file=str(source.relative_to(ROOT)),
                            sha256=digest, bytes=source.stat().st_size,
                            http_status=meta.get("http_status"),
                            retrieved_utc=meta.get("retrieved_utc")))
    return records


def bp_rp_excess_correction(color):
    """Riello et al. 2021, A&A 649 A3 Eq. 6; not a binarity probability."""
    if color < 0.5:
        return 1.154360 + 0.033772 * color + 0.032277 * color**2
    if color < 4:
        return 1.162004 + 0.011464 * color + 0.049255 * color**2 - 0.005879 * color**3
    return 1.057572 + 0.140537 * color


def derive():
    RESULT.mkdir(parents=True, exist_ok=True)
    gaia = {r["source_id"]: r for r in csv_rows(ROOT / "data/stellar_systems/prior_audit/supplement_gaia_sources.csv")}
    aps = {r["source_id"]: r for r in csv_rows(OUT / "gaia_astrophysical_parameters.csv")}
    supp = {r["source_id"]: r for r in csv_rows(OUT / "gaia_astrophysical_parameters_supp.csv")}
    tic_json = json.loads((ROOT / "data/stellar_systems/prior_audit/supplement_tic_sources.json").read_text())
    if tic_json.get("status") != "COMPLETE":
        raise RuntimeError("Incomplete original TIC query")
    tics = {str(r["GAIA"]): r for r in tic_json["data"]}
    context = []
    gfields = ["ra", "dec", "parallax", "parallax_error", "pmra", "pmra_error",
               "pmdec", "pmdec_error", "pmra_pmdec_corr", "phot_g_mean_mag",
               "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp", "ruwe",
               "ipd_frac_multi_peak", "phot_bp_rp_excess_factor", "radial_velocity",
               "radial_velocity_error", "rv_method_used", "rv_nb_transits",
               "rv_chisq_pvalue", "rv_renormalised_gof", "rv_amplitude_robust"]
    for sid in IDS:
        row = gaia[sid]
        p = number(row, "parallax")
        rec = dict(gaia_dr3_source_id=sid, tic=tics[sid]["ID"],
                   gaia={key: number(row, key) for key in gfields},
                   tic_catalog={key: tics[sid].get(key) for key in ["Tmag", "Teff", "rad", "mass", "d", "ebv"]},
                   gaia_astrophysical_parameters=aps[sid],
                   gaia_astrophysical_parameters_supp=supp.get(sid, {}),
                   parallax_inverse_pc=1000 / p,
                   parallax_inverse_error_pc=1000 * number(row, "parallax_error") / p**2,
                   absolute_g_without_extinction=number(row, "phot_g_mean_mag") + 5 * math.log10(p) - 10,
                   bp_rp_corrected_flux_excess=number(row, "phot_bp_rp_excess_factor") - bp_rp_excess_correction(number(row, "bp_rp")))
        context.append(rec)
    giant, faint = gaia[IDS[1]], gaia[IDS[2]]
    gc = SkyCoord(number(giant, "ra") * u.deg, number(giant, "dec") * u.deg)
    fc = SkyCoord(number(faint, "ra") * u.deg, number(faint, "dec") * u.deg)
    separation = gc.separation(fc).arcsec
    d = 1000 / number(giant, "parallax")
    cov = np.zeros((2, 2))
    for row in [giant, faint]:
        se = np.array([number(row, "pmra_error"), number(row, "pmdec_error")])
        cov += np.diag(se**2)
        cov[0, 1] += number(row, "pmra_pmdec_corr") * np.prod(se)
    cov[1, 0] = cov[0, 1]
    delta = np.array([number(faint, "pmra") - number(giant, "pmra"),
                      number(faint, "pmdec") - number(giant, "pmdec")])
    q = float(delta @ np.linalg.solve(cov, delta))
    pair = dict(separation_arcsec=separation, distance_assumed_pc=d,
                projected_separation_au=separation * d,
                delta_parallax_mas=number(faint, "parallax") - number(giant, "parallax"),
                parallax_difference_standard_errors=(number(faint, "parallax") - number(giant, "parallax")) / math.hypot(number(faint, "parallax_error"), number(giant, "parallax_error")),
                delta_proper_motion_masyr=delta.tolist(),
                delta_proper_motion_covariance_masyr2=cov.tolist(),
                proper_motion_mahalanobis_chisq_2d=q,
                proper_motion_zero_difference_chisq_pvalue=float(chi2.sf(q, 2)),
                nominal_transverse_relative_speed_kms=float(4.74047 * d * np.linalg.norm(delta) / 1000),
                projected_separation_escape_speeds_kms={str(m): float(42.1219 * math.sqrt(m / (separation * d))) for m in [1, 2, 3]},
                caveat="Assumes same distance and independent source errors; no binding probability, shared astrometric systematics or line-of-sight velocity.")
    # Differential aperture throughput is unknown. These are conditional
    # flux-budget checks, not measured intrinsic eclipse depths.
    flux_budget = []
    for sid in [IDS[2], IDS[3]]:
        delta_t = tics[sid]["Tmag"] - tics[IDS[1]]["Tmag"]
        ratio = 10**(-0.4 * delta_t)
        fraction = ratio / (1 + ratio)
        flux_budget.append(dict(tic=tics[sid]["ID"], delta_tmag=delta_t,
            faint_to_giant_tess_flux_ratio=ratio, faint_fraction_of_two_source_light=fraction,
            intrinsic_fraction_needed_for_1710ppm_primary=0.001710 / fraction,
            intrinsic_fraction_needed_for_972ppm_window_primary=0.000972 / fraction,
            intrinsic_fraction_needed_for_445ppm_window_secondary=0.000445 / fraction,
            caveat="Equal aperture coupling; approximate catalog T magnitudes; other blend light and calibration uncertainties omitted."))
    controls = {}
    spectra = {}
    for label in SPECTRO:
        name = label + "_positive_cone_retry1.tsv" if label == "apogee17" else label + "_positive_cone.tsv"
        control = prior.parse_tsv(OUT / name, with_table=True)
        if not control:
            raise ValueError("Failed positive cone control: " + label)
        controls[label] = len(control)
        spectra[label] = {str(t): len(prior.parse_tsv(OUT / f"tic{t}_{label}.tsv", with_table=True)) for t in TARGETS}
    known = csv_rows(OUT / "wide_catalog_positive_recovery.csv")
    if len(known) != 1 or known[0]["Source1"] != "3770759051395978240":
        raise ValueError("Wide-catalog exact-source control failed")
    wide = dict(target_exact_id_matches=len(csv_rows(OUT / "wide_catalog_targets.csv")),
                positive_exact_id_matches=1,
                applicability="El-Badry2021 selects parallax>1mas. Source273 is in nominal distance domain; giant and both close neighbors are outside it. Null giant-field matches have no recognition value.")
    prior_flags = prior.parse_tsv(ROOT / "data/stellar_systems/prior_audit/tic273879360_10k2025.tsv", with_table=True)
    prior_flags = [r for r in prior_flags if r["vizier_table"] == "J/ApJS/279/50/table2" and r["TIC"] == "273879360"]
    if len(prior_flags) != 1:
        raise ValueError("Prior exact-TIC eclipse-flag row not recovered")
    result = dict(derived_utc=utc(), sources=verify(), target_context=context,
                  close_pair_astrometry=pair, conditional_flux_budget=flux_budget,
                  positive_cone_controls=controls, spectroscopy_target_cone_counts=spectra,
                  elbadry2021_wide_catalog=wide,
                  prior_eclipse_flag=dict(source="Kostov2025 Table2 / arXiv2506.05631",
                    exact_tic_rows=prior_flags,
                    interpretation="Previously published neural-network eclipse-like-feature flag, score>0.9. Unvetted/unvalidated, no period or EB classification. Not a generic observed-target inventory; first recognition of eclipse-like variability cannot be claimed."))
    (RESULT / "context_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    with (RESULT / "source_properties.csv").open("w", newline="") as stream:
        rows = []
        for r in context:
            g = r["gaia"]
            rows.append(dict(tic=r["tic"], gaia_dr3_source_id=r["gaia_dr3_source_id"],
                gmag=g["phot_g_mean_mag"], tmag=r["tic_catalog"]["Tmag"],
                bp_rp=g["bp_rp"], parallax_mas=g["parallax"], parallax_error_mas=g["parallax_error"],
                inverse_parallax_pc=r["parallax_inverse_pc"],
                absolute_g_no_extinction=r["absolute_g_without_extinction"],
                corrected_bp_rp_flux_excess=r["bp_rp_corrected_flux_excess"],
                ruwe=g["ruwe"], rv_kms=g["radial_velocity"], rv_error_kms=g["radial_velocity_error"],
                rv_method=g["rv_method_used"], rv_chisq_pvalue=g["rv_chisq_pvalue"]))
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps(dict(source_snapshots=len(result["sources"]),
                         controls=controls, pair=pair, flux_budget=flux_budget), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    if args.download:
        download()
    derive()
