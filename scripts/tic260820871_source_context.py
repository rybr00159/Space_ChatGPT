"""Cache exact public stellar/source responses for TIC 260820871.

Original successful response bodies are reused, never overwritten. The
separate variable-context report covers SIMBAD and variable-star catalogs.
"""
import concurrent.futures
import csv
import hashlib
import io
import json
import math
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "broad_catalogs"
TIC = 260820871
RA, DEC = 348.054595469913, -70.1131990097929


def fetch(label, url, body=None):
    path = BASE / f"tic{TIC}_{label}"
    if path.exists():
        return path.read_bytes()
    with urllib.request.urlopen(urllib.request.Request(url, data=body), timeout=90) as r:
        data = r.read()
        metadata = {"retrieved_utc": datetime.now(timezone.utc).isoformat(),
                    "url": url, "method": "POST" if body else "GET",
                    "http_status": r.status, "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "method_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if body:
        metadata["request_form"] = urllib.parse.parse_qs(body.decode())
    path.write_bytes(data)
    Path(str(path) + ".manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return data


def gaia(label, query):
    params = {"request": "doQuery", "lang": "ADQL", "format": "csv", "query": query}
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode(params)
    return fetch(label + ".csv", url)


def separation(ra, dec):
    r, d, tr, td = map(math.radians, (ra, dec, RA, DEC))
    v = math.sin((d - td) / 2) ** 2 + math.cos(d) * math.cos(td) * math.sin((r - tr) / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(min(v, 1)))) * 3600


def main():
    req = {"service": "Mast.Catalogs.Tic.Cone", "params": {"ra": RA, "dec": DEC, "radius": 1 / 30},
           "format": "json", "pagesize": 1000, "page": 1}
    q = "select source_id,ra,dec,phot_g_mean_mag,bp_rp,parallax,parallax_error,pmra,pmdec," \
        "ruwe,duplicated_source,astrometric_params_solved,astrometric_excess_noise," \
        "astrometric_excess_noise_sig,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude," \
        "phot_bp_rp_excess_factor,non_single_star,phot_variable_flag,radial_velocity," \
        "radial_velocity_error,rv_nb_transits,rv_amplitude_robust,rv_chisq_pvalue," \
        "rv_renormalised_gof from gaiadr3.gaia_source where 1=CONTAINS(" \
        f"POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.0333333333333))"
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        fs = [pool.submit(fetch, "exofop_current.json", f"https://exofop.ipac.caltech.edu/tess/target.php?id={TIC}&json"),
              pool.submit(fetch, "tic_cone120.json", "https://mast.stsci.edu/api/v0/invoke",
                          urllib.parse.urlencode({"request": json.dumps(req)}).encode()),
              pool.submit(gaia, "gaia_neighbors_current", q)]
        for f in fs:
            f.result()
    rows = list(csv.DictReader(io.StringIO((BASE / f"tic{TIC}_gaia_neighbors_current.csv").read_text())))
    target = min(rows, key=lambda r: separation(float(r["ra"]), float(r["dec"])))
    if separation(float(target["ra"]), float(target["dec"])) > 2:
        raise RuntimeError("No Gaia source within 2 arcsec of TIC position")
    source_id = target["source_id"]
    print("Gaia target", source_id, target, "cone rows", len(rows), flush=True)
    queries = {
        "gaia_stellar_current": "select source_id,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper," \
            "logg_gspphot,radius_gspphot,radius_flame,mass_flame from gaiadr3.astrophysical_parameters " \
            f"where source_id={source_id}",
        "gaia_class_current": f"select * from gaiadr3.vari_classifier_result where source_id={source_id}",
        "gaia_eb_current": f"select * from gaiadr3.vari_eclipsing_binary where source_id={source_id}",
        "gaia_rotation_current": f"select * from gaiadr3.vari_rotation_modulation where source_id={source_id}",
        "gaia_ms_oscillator_current": f"select * from gaiadr3.vari_ms_oscillator where source_id={source_id}",
        "gaia_short_timescale_current": f"select * from gaiadr3.vari_short_timescale where source_id={source_id}",
    }
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        fs = {pool.submit(gaia, label, q): label for label, q in queries.items()}
        for f in concurrent.futures.as_completed(fs):
            data = f.result()
            print(fs[f], data.decode()[:1000], flush=True)
    sources = []
    for p in sorted(BASE.glob(f"tic{TIC}_*.manifest.json")):
        raw = Path(str(p)[:-len(".manifest.json")])
        metadata = json.loads(p.read_text())
        digest = hashlib.sha256(raw.read_bytes()).hexdigest()
        if digest != metadata["sha256"]:
            raise RuntimeError("Source integrity failure: " + str(raw))
        sources.append({"file": str(raw.relative_to(ROOT)), "sha256": digest,
                        "bytes": raw.stat().st_size, "manifest": str(p.relative_to(ROOT))})
    output = ROOT / "reports" / f"tic{TIC}_source_context_sources.json"
    output.write_text(json.dumps({"checked_utc": datetime.now(timezone.utc).isoformat(),
                                 "verified_source_count": len(sources), "sources": sources}, indent=2) + "\n")
    print("Verified", len(sources), "public source bodies", flush=True)


if __name__ == "__main__":
    main()
