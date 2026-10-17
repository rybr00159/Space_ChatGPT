"""Finite public source/catalogue checks for TIC 214311908, no paid calls.

Retains exact responses and URL/query/UTC/checksum manifests. Successful cached
responses are verified and reused; failure bodies are preserved on retries.
Catalogue absence and a finite alias search cannot establish novelty.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError
import csv
import hashlib
import io
import json
import math
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/tic214311908/catalogs"
SUMMARY = ROOT / "candidates/tic214311908/catalogs"
TIC = 214311908
SID = "5925373272850784768"
RA, DEC = 260.96969903043, -51.3238331382385
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
NASA = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


def utc():
    return datetime.now(timezone.utc).isoformat()


def fetch(name, url, query=None):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    meta_path = OUT / (name + ".manifest.json")
    # A preserved original failure can have a later successful immutable retry.
    for candidate in sorted(OUT.glob(path.stem + "*.manifest.json"), reverse=True):
        old = json.loads(candidate.read_text())
        cached = candidate.with_name(candidate.name.removesuffix(".manifest.json"))
        if old.get("http_status")==200 and old.get("url")==url and cached.exists():
            if hashlib.sha256(cached.read_bytes()).hexdigest()!=old["sha256"]:
                raise ValueError(f"Cached checksum mismatch: {cached}")
            return cached
    if path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest() != meta["sha256"]:
            raise ValueError(f"Cached checksum mismatch: {path}")
        if meta["http_status"] == 200 and meta["url"] == url:
            return path
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        path = path.with_name(path.stem + "_" + stamp + path.suffix)
        meta_path = path.with_suffix(path.suffix + ".manifest.json")
    try:
        response = urlopen(Request(url, headers={"User-Agent": "Space_ChatGPT-public-research/0.6"}), timeout=100)
    except HTTPError as error:
        response = error
    with response:
        body = response.read()
        status = response.status
        content_type = response.headers.get("Content-Type")
    path.write_bytes(body)
    meta_path.write_text(json.dumps(dict(url=url, query=query, retrieved_utc=utc(),
        http_status=status, bytes=len(body), sha256=hashlib.sha256(body).hexdigest(),
        content_type=content_type), indent=2) + "\n")
    if status != 200:
        raise RuntimeError(f"HTTP {status} for {name}; exact failure retained")
    print(name, status, len(body), flush=True)
    return path


def tap(name, query, endpoint=GAIA):
    params = dict(query=query, format="csv")
    if endpoint != NASA:
        params.update(request="doQuery", lang="ADQL")
    return fetch(name, endpoint + "?" + urlencode(params), query)


def rows(path):
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))


def main():
    SUMMARY.mkdir(parents=True, exist_ok=True)
    cols = "source_id,ra,dec,phot_g_mean_mag,bp_rp,ruwe,non_single_star,phot_variable_flag,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,astrometric_excess_noise,astrometric_excess_noise_sig,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,duplicated_source,radial_velocity,radial_velocity_error,rv_nb_transits,rv_amplitude_robust,rv_chisq_pvalue"
    cone = f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.0333333333333333))"
    queries = {
        "gaia_target.csv": f"select {cols} from gaiadr3.gaia_source where source_id={SID}",
        "gaia_neighbors.csv": f"select {cols} from gaiadr3.gaia_source where {cone} order by phot_g_mean_mag",
        "gaia_stellar.csv": f"select source_id,teff_gspphot,teff_gspphot_lower,teff_gspphot_upper,logg_gspphot,logg_gspphot_lower,logg_gspphot_upper,radius_gspphot,radius_gspphot_lower,radius_gspphot_upper,radius_flame,mass_flame from gaiadr3.astrophysical_parameters where source_id={SID}",
        "gaia_class.csv": f"select * from gaiadr3.vari_classifier_result where source_id={SID}",
        "gaia_eb.csv": f"select * from gaiadr3.vari_eclipsing_binary where source_id={SID}",
        "gaia_rotation.csv": f"select * from gaiadr3.vari_rotation_modulation where source_id={SID}",
        "gaia_nss_orbit.csv": f"select * from gaiadr3.nss_two_body_orbit where source_id={SID}",
        "gaia_nss_acceleration.csv": f"select * from gaiadr3.nss_acceleration_astro where source_id={SID}",
        "gaia_nss_spectroscopic.csv": f"select * from gaiadr3.nss_non_linear_spectro where source_id={SID}",
    }
    for name, table in [("eb", "vari_eclipsing_binary"),
                        ("rotation", "vari_rotation_modulation"),
                        ("lpv", "vari_long_period_variable"),
                        ("class", "vari_classifier_result")]:
        queries[f"gaia_field_{name}.csv"] = f"select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,v.* from gaiadr3.gaia_source g join gaiadr3.{table} v on g.source_id=v.source_id where 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{RA},{DEC},0.0333333333333333))"
    queries["gaia_eb_epoch_definition.csv"] = "select table_name,column_name,description,unit from tap_schema.columns where table_name='gaiadr3.vari_eclipsing_binary' and column_name in ('reference_time','frequency','geom_model_gaussian1_phase','derived_primary_ecl_phase')"
    jobs = [(name, tap, (name, query)) for name, query in queries.items()]
    for name, query, endpoint in [
        ("nasa_confirmed.csv", f"select pl_name,hostname,tic_id,gaia_dr2_id,ra,dec,discoverymethod,pl_orbper,pl_tranmid,pl_trandur,pl_rade,pl_orbper_reflink from pscomppars where tic_id='TIC {TIC}' or {cone}", NASA),
        # The NASA TOI backend currently returns ORA-00904 HTM20 for CONTAINS.
        # Preserve that failed response and use a conservative bounding box.
        ("nasa_toi.csv", f"select * from toi where tid={TIC} or (ra between {RA-.054} and {RA+.054} and dec between {DEC-.034} and {DEC+.034})", NASA),
        ("simbad_cone.csv", f"select oid,main_id,ra,dec,otype,nbref from basic where {cone}", SIMBAD),
    ]:
        jobs.append((name, tap, (name, query, endpoint)))
    urls = {
        "exofop.json": f"https://exofop.ipac.caltech.edu/tess/target.php?id={TIC}&json",
        "exofop_nearby.csv": f"https://exofop.ipac.caltech.edu/tess/download_nearbytarget.php?id={TIC}&output=csv",
        "exofop_ctoi_all.csv": "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?output=csv",
        "tesseb.html": f"https://tessebs.villanova.edu/{TIC:010d}",
        "vsx.tsv": "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source":"B/vsx/vsx", "-c":f"{RA} {DEC}","-c.rs":120,"-out.max":"unlimited","-out.all":"","-out.form":"tsv"}),
        "raven_zenodo.json": "https://zenodo.org/api/records/19661443",
        "gaia_eb_model_documentation.html": "https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_variability_tables/ssec_dm_vari_eclipsing_binary.html",
    }
    for name, source in [("elbadry_binary_cone.tsv", "J/MNRAS/506/2269"),
                         ("wds_cone.tsv", "B/wds/wds")]:
        urls[name] = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source":source,"-c":f"{RA} {DEC}","-c.rs":30,"-out.max":"unlimited","-out.all":"","-out.form":"tsv"})
    jobs += [(name, fetch, (name, url)) for name,url in urls.items()]
    paths, failures = {}, {}
    with ThreadPoolExecutor(max_workers=5) as pool:
        future_map = {pool.submit(fn,*args):name for name,fn,args in jobs}
        for future in as_completed(future_map):
            name = future_map[future]
            try:
                paths[name] = str(future.result().relative_to(ROOT))
            except Exception as error:
                failures[name] = repr(error)
                print(name, repr(error), flush=True)
    target = rows(ROOT / paths["gaia_target.csv"])[0] if "gaia_target.csv" in paths else None
    neighbors = []
    if target and "gaia_neighbors.csv" in paths:
        for row in rows(ROOT / paths["gaia_neighbors.csv"]):
            dra = math.radians(float(row["ra"]) - RA)
            a = math.sin(math.radians(float(row["dec"])-DEC)/2)**2 + math.cos(math.radians(DEC))*math.cos(math.radians(float(row["dec"])))*math.sin(dra/2)**2
            row["separation_arcsec"] = math.degrees(2*math.asin(math.sqrt(a)))*3600
            row["delta_g"] = float(row["phot_g_mean_mag"]) - float(target["phot_g_mean_mag"]) if row["phot_g_mean_mag"] else None
            row["g_flux_ratio_to_target"] = 10**(-.4*row["delta_g"]) if row["delta_g"] is not None else None
            row["is_target"] = row["source_id"] == SID
            neighbors.append(row)
        neighbors.sort(key=lambda row: row["separation_arcsec"])
        with (SUMMARY / "gaia_neighbors_derived.csv").open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=list(neighbors[0]));writer.writeheader();writer.writerows(neighbors)
    simrows = rows(ROOT / paths["simbad_cone.csv"]) if "simbad_cone.csv" in paths else []
    for row in simrows:
        row["separation_arcsec"] = math.hypot((float(row["ra"])-RA)*math.cos(math.radians(DEC)),float(row["dec"])-DEC)*3600
    simhost = min(simrows,key=lambda row: row["separation_arcsec"]) if simrows else None
    follow = []
    if simhost and simhost["separation_arcsec"] < 3:
        oid = simhost["oid"]
        follow = [("simbad_aliases.csv",f"select id from ident where oidref={oid}"),
                  ("simbad_references.csv",f"select r.bibcode,r.doi,r.title,r.abstract from has_ref h join ref r on h.oidbibref=r.oidbib where h.oidref={oid}")]
    for name,query in follow:
        try: paths[name] = str(tap(name,query,SIMBAD).relative_to(ROOT))
        except Exception as error: failures[name] = repr(error)
    aliases = ["TIC 214311908","214311908","2MASS J17235292-5119271","5925373272850784768","UCAC4 194-162330"]
    if "simbad_aliases.csv" in paths:
        aliases += [r["id"] for r in rows(ROOT / paths["simbad_aliases.csv"])]
    aliases = list(dict.fromkeys(aliases))
    query = " OR ".join('all:"'+alias+'"' for alias in aliases)
    try:
        path = fetch("arxiv_aliases.xml","https://export.arxiv.org/api/query?"+urlencode(dict(search_query=query,max_results=100)),query)
        paths["arxiv_aliases.xml"] = str(path.relative_to(ROOT))
    except Exception as error: failures["arxiv_aliases.xml"] = repr(error)
    summary = dict(tic=TIC,gaia_source_id=SID,checked_utc=utc(),ra=RA,dec=DEC,
        paths=paths,failures=failures,target=target,simbad_host=simhost,finite_aliases=aliases,
        neighbor_count=len(neighbors),nearest_non_target_neighbors=[n for n in neighbors if not n["is_target"]][:12],
        scope="Exact public snapshots, 120arcsec source cone and finite SIMBAD references/arXiv alias checks. Does not establish novelty, locate TESS signal, validate a planet, or rule out binarity. G flux ratios are not TESS aperture coupling.")
    for name in ["nasa_confirmed.csv","nasa_toi.csv"]:
        summary[name] = rows(ROOT/paths[name]) if name in paths else None
    if "exofop_ctoi_all.csv" in paths:
        allrows = rows(ROOT/paths["exofop_ctoi_all.csv"])
        summary["ctoi_snapshot_rows"] = len(allrows)
        summary["ctoi_exact_tic"] = [r for r in allrows if r.get("TIC ID")==str(TIC)]
    summary["gaia_target_table_counts"] = {name:len(rows(ROOT/path)) for name,path in paths.items()
        if name.startswith("gaia_") and name.endswith(".csv") and not name.startswith("gaia_field")}
    if "raven_zenodo.json" in paths:
        record = json.loads((ROOT/paths["raven_zenodo.json"]).read_text())
        summary["raven_release"] = dict(doi=record.get("doi"),record_id=record.get("id"),tables={})
        for item in record["files"]:
            if item["key"].endswith(".csv") and item["size"]<5_000_000:
                name = "raven_"+item["key"]
                try:
                    path = fetch(name,item["links"]["self"])
                    if "md5:"+hashlib.md5(path.read_bytes()).hexdigest()!=item["checksum"]:
                        raise ValueError("RAVEN file differs from Zenodo record MD5")
                    paths[name] = str(path.relative_to(ROOT))
                    rr = rows(path)
                    summary["raven_release"]["tables"][item["key"]] = dict(rows=len(rr),exact_tic_matches=[r for r in rr if r.get("ticid")==str(TIC)],md5_verified=True)
                except Exception as error: failures[name] = repr(error)
    if "arxiv_aliases.xml" in paths:
        feed = ET.fromstring((ROOT/paths["arxiv_aliases.xml"]).read_bytes())
        ns = {"a":"http://www.w3.org/2005/Atom", "op":"http://a9.com/-/spec/opensearch/1.1/"}
        summary["arxiv_total_results"] = feed.findtext("op:totalResults",default=None,namespaces=ns)
        summary["arxiv_entries"] = [dict(id=e.findtext("a:id",namespaces=ns),title=e.findtext("a:title",namespaces=ns)) for e in feed.findall("a:entry",ns)]
    if "gaia_field_eb.csv" in paths:
        lookup = {r["source_id"]:r for r in neighbors}
        ephemerides = []
        for r in rows(ROOT/paths["gaia_field_eb.csv"]):
            frequency = float(r["frequency"])
            period = 1/frequency
            ephemerides.append(dict(source_id=r["source_id"],ra=r["ra"],dec=r["dec"],
                separation_arcsec=lookup[r["source_id"]]["separation_arcsec"],gmag=r["phot_g_mean_mag"],
                period_days=period,period_error_days=float(r["frequency_error"])/frequency**2,
                primary_reference_epoch_jd_tcb=2455197.5+float(r["reference_time"])+float(r["derived_primary_ecl_phase"])*period,
                reference_time_jd_tcb_minus_2455197_5=r["reference_time"],
                primary_phase=r["derived_primary_ecl_phase"],primary_phase_error=r["derived_primary_ecl_phase_error"],
                secondary_phase=r["derived_secondary_ecl_phase"],primary_depth_magnitudes=r["derived_primary_ecl_depth"],
                secondary_depth_magnitudes=r["derived_secondary_ecl_depth"],
                note="Published Gaia geometric EB model; convert TCB to TDB for exact TESS comparison and propagate long-baseline ephemeris uncertainty. This source is not TIC214311908."))
        (SUMMARY/"known_neighbor_eb_ephemerides.json").write_text(json.dumps(ephemerides,indent=2)+"\n")
    (SUMMARY/"catalog_check.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(dict(summary=str(SUMMARY/"catalog_check.json"),failures=failures),indent=2))


if __name__ == "__main__":
    main()
