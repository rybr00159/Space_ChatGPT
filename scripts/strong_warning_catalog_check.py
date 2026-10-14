"""Finite, reproducible public source/prior checks for three warning leads.

No paid calls. Frozen survey/model/ranking inputs are only read. Every response
is retained with URL/query/UTC/checksum; an existing snapshot is not replaced.
Catalog and finite alias-search absence cannot establish novelty.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
import csv
import hashlib
import io
import json
import math
import re
import sys
import xml.etree.ElementTree as ET

import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from broad_crossmatch import CatalogIndex

OUT = ROOT / "data" / "broad_catalogs"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
TARGETS = [
    (300814564, 116.922122, -68.62146, 1.861107),
    (273879358, 32.307882, -71.320039, 7.58187),
    (270312735, 32.629476, -83.14817, 8.373627),
]
GAIA_COLUMNS = "source_id,ra,dec,phot_g_mean_mag,bp_rp,ruwe,non_single_star,phot_variable_flag,astrometric_params_solved,parallax,parallax_error,astrometric_excess_noise,astrometric_excess_noise_sig,duplicated_source,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,phot_bp_rp_excess_factor,phot_g_mean_flux_over_error,phot_g_n_obs,radial_velocity,radial_velocity_error,rv_nb_transits,rv_amplitude_robust,rv_chisq_pvalue,rv_renormalised_gof,rv_template_teff,rv_expected_sig_to_noise"


def utc():
    return datetime.now(timezone.utc).isoformat()


def save(name, body, meta):
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / name
    if p.exists():
        p = p.with_name(p.stem + "_" + STAMP + p.suffix)
    if p.exists():
        raise RuntimeError(f"Refuse overwrite: {p}")
    p.write_bytes(body)
    meta.update(sha256=hashlib.sha256(body).hexdigest(), bytes=len(body))
    p.with_suffix(p.suffix + ".manifest.json").write_text(json.dumps(meta, indent=2) + "\n")
    return p


def fetch(tic, label, url, query=None):
    meta = dict(retrieved_utc=utc(), url=url)
    if query is not None:
        meta["query"] = query
    try:
        r = requests.get(url, timeout=(20, 80), headers={"User-Agent": "Space_ChatGPT-public-research/0.4"})
        meta.update(http_status=r.status_code, content_type=r.headers.get("Content-Type"))
        p = save(f"tic{tic}_{label}", r.content, meta)
        r.raise_for_status()
        print(tic, label, "OK", r.status_code, len(r.content), flush=True)
        return p
    except Exception as e:
        print(tic, label, "ERROR", type(e).__name__, str(e)[:140], flush=True)
        return None


def tap(tic, label, query, service="gaia"):
    base = {"gaia": "https://gea.esac.esa.int/tap-server/tap/sync", "simbad": "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"}[service]
    url = base + "?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    return fetch(tic, label + ".csv", url, query)


def handle(target, ix):
    tic, ra, dec, period = target
    match = dict(checked_utc=utc(), tic=tic, ra=ra, dec=dec, period_days=period,
                 scope="Cached all-method catalog target/120-arcsec position and exact/2–5x integer harmonic/subharmonic periods within 1%. Does not establish novelty.",
                 catalog_metadata=ix.metadata(), signal=ix.match_signal(tic, period, ra, dec))
    mpath = save(f"tic{tic}_cached_matches_current.json", (json.dumps(match, indent=2) + "\n").encode(),
                dict(derived_utc=utc(), sources=ix.metadata()))
    if match["signal"]["matches"]:
        print("URGENT CATALOG", tic, json.dumps(match["signal"]), flush=True)
    exopath = fetch(tic, "exofop_current.json", f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json")
    if not exopath:
        return dict(tic=tic, cached_matches=str(mpath), error="No ExoFOP response")
    j = json.loads(exopath.read_text())
    aliases = j["basic_info"].get("star_names", "")
    gaia = re.search(r"Gaia DR3 (\d+)", aliases)
    ra, dec = float(j["coordinates"]["ra"]), float(j["coordinates"]["dec"])
    jobs = [
        (fetch, (tic, "exofop_nearby_current.csv", f"https://exofop.ipac.caltech.edu/tess/download_nearbytarget.php?id={tic}&output=csv")),
        (fetch, (tic, "tesseb_live.html", f"https://tessebs.villanova.edu/{tic:010d}")),
        (fetch, (tic, "vsx_current.tsv", "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source":"B/vsx/vsx", "-c": f"{ra} {dec}", "-c.rs":120,"-out.max":"unlimited","-out.all":"","-out.form":"tsv"}))),
        (tap, (tic, "simbad_current", f"select oid,main_id,ra,dec,otype,nbref from basic where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333))", "simbad")),
    ]
    if gaia:
        sid = gaia.group(1)
        jobs += [(tap, (tic, label, query)) for label, query in [
            ("gaia_target_current", f"select {GAIA_COLUMNS} from gaiadr3.gaia_source where source_id={sid}"),
            ("gaia_stellar_current", f"select source_id,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,radius_gspphot,radius_gspphot_lower,radius_gspphot_upper,radius_flame,mass_flame from gaiadr3.astrophysical_parameters where source_id={sid}"),
            ("gaia_class_current", f"select source_id,best_class_name,best_class_score from gaiadr3.vari_classifier_result where source_id={sid}"),
            ("gaia_eb_current", f"select * from gaiadr3.vari_eclipsing_binary where source_id={sid}"),
            ("gaia_neighbors_current", f"select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.bp_rp,g.ruwe,g.non_single_star,g.phot_variable_flag,v.in_vari_eclipsing_binary from gaiadr3.gaia_source g left outer join gaiadr3.vari_summary v on g.source_id=v.source_id where 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333))"),
        ]]
    files = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [(args[1], pool.submit(fn, *args)) for fn, args in jobs]
        for label, fut in futures:
            p = fut.result()
            files[label] = str(p) if p else None
    # Closest SIMBAD association only for aliases/reference follow-up; all cone
    # rows retained. The positional match itself is not signal attribution.
    sfile = files.get("simbad_current")
    simrows = list(csv.DictReader(Path(sfile).open())) if sfile else []
    def separation(r):
        return math.hypot((float(r["ra"])-ra)*math.cos(math.radians(dec)), float(r["dec"])-dec)*3600
    simrows.sort(key=separation)
    simhost = simrows[0] if simrows and separation(simrows[0]) <= 2 else None
    if simhost:
        oid = int(simhost["oid"])
        follow = [
            ("simbad_current_aliases", f"select id from ident where oidref={oid}"),
            ("simbad_current_references", f"select r.bibcode,r.doi,r.title,r.abstract from has_ref h join ref r on h.oidbibref=r.oidbib where h.oidref={oid}"),
        ]
        with ThreadPoolExecutor(max_workers=2) as pool:
            for label, fut in [(lab, pool.submit(tap,tic,lab,q,"simbad")) for lab,q in follow]:
                p = fut.result()
                files[label] = str(p) if p else None
    finite_aliases = [a.strip() for a in aliases.split(",") if a.strip()]
    if files.get("simbad_current_aliases"):
        finite_aliases += [r["id"] for r in csv.DictReader(Path(files["simbad_current_aliases"]).open())]
    finite_aliases = list(dict.fromkeys(finite_aliases))
    query = " OR ".join('all:"' + a + '"' for a in finite_aliases)
    p = fetch(tic, "arxiv_current.xml", "https://export.arxiv.org/api/query?" + urlencode(dict(search_query=query, max_results=100)), query)
    files["arxiv_current"] = str(p) if p else None
    result = dict(tic=tic, period_days=period, checked_utc=utc(), files=files, cached_matches=str(mpath),
                  exofop_file=str(exopath), finite_aliases=finite_aliases, simbad_host=simhost,
                  scope="Finite SIMBAD-reference and arXiv metadata alias search, not exhaustive journal/full-text/supplementary/literature review. Source and EB checks do not assign photometric signal.")
    if p:
        ns={"op":"http://a9.com/-/spec/opensearch/1.1/", "a":"http://www.w3.org/2005/Atom"}
        aroot=ET.fromstring(p.read_bytes())
        result["arxiv_total_results"] = aroot.findtext("op:totalResults", "?", ns)
        result["arxiv_entries"] = [dict(id=e.findtext("a:id","",ns), title=e.findtext("a:title","",ns)) for e in aroot.findall("a:entry",ns)]
    save(f"tic{tic}_finite_source_check_current.json",(json.dumps(result,indent=2)+"\n").encode(),dict(derived_utc=utc(),sources=files))
    return result


def derive_neighbors(result):
    """Preserve exact angular distances/G ratios without assigning signal source."""
    tic = result["tic"]
    p = result.get("files", {}).get("gaia_neighbors_current")
    target_path = result.get("files", {}).get("gaia_target_current")
    if not p or not target_path:
        return
    j = json.loads(Path(result["exofop_file"]).read_text())
    center = SkyCoord(float(j["coordinates"]["ra"]) * u.deg, float(j["coordinates"]["dec"]) * u.deg)
    host = next(csv.DictReader(Path(target_path).open()))
    rows = list(csv.DictReader(Path(p).open()))
    for row in rows:
        row["separation_arcsec"] = float(center.separation(SkyCoord(float(row["ra"])*u.deg, float(row["dec"])*u.deg)).arcsec)
        row["delta_g"] = float(row["phot_g_mean_mag"])-float(host["phot_g_mean_mag"]) if row["phot_g_mean_mag"] else None
        row["g_flux_ratio"] = 10**(-0.4*row["delta_g"]) if row["delta_g"] is not None else None
    rows.sort(key=lambda row: row["separation_arcsec"])
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    save(f"tic{tic}_gaia_neighbors_derived_current.csv",buffer.getvalue().encode(),dict(derived_utc=utc(),
         sources={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest()},
         scope="120arcsec Gaia cone distances from current ExoFOP coordinates; G flux ratios are not TESS/aperture contamination or source attribution."))


def follow_close_bright_sources():
    """Inspect two real Gaia neighbors suggested by the preserved source cones."""
    jobs = [
        (tap,(300814564,"gaia_close_neighbor_current",f"select {GAIA_COLUMNS},pmra,pmra_error,pmdec,pmdec_error from gaiadr3.gaia_source where source_id=5273775091332140800")),
        (tap,(300814564,"gaia_target_motion_current","select source_id,pmra,pmra_error,pmdec,pmdec_error from gaiadr3.gaia_source where source_id=5273775087036017024")),
        (fetch,(300814564,"close_neighbor765053491_exofop_current.json","https://exofop.ipac.caltech.edu/tess/target.php?id=765053491&json")),
        (fetch,(300814564,"close_neighbor765053491_tesseb_current.html","https://tessebs.villanova.edu/0765053491")),
        (tap,(270312735,"gaia_bright_neighbor_current",f"select {GAIA_COLUMNS} from gaiadr3.gaia_source where source_id=4617279616331585408")),
        (fetch,(270312735,"bright_neighbor270312732_exofop_current.json","https://exofop.ipac.caltech.edu/tess/target.php?id=270312732&json")),
        (fetch,(270312735,"bright_neighbor270312732_tesseb_current.html","https://tessebs.villanova.edu/0270312732")),
    ]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in [pool.submit(fn,*args) for fn,args in jobs]:
            future.result()


if __name__ == "__main__":
    ix = CatalogIndex(position_radius_arcsec=120)
    with ThreadPoolExecutor(max_workers=3) as pool:
        for result in pool.map(lambda target: handle(target, ix), TARGETS):
            derive_neighbors(result)
            print("DONE", result["tic"], json.dumps(result.get("simbad_host")), result.get("arxiv_total_results"), flush=True)
    follow_close_bright_sources()
