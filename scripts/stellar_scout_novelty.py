"""Finite, immutable recognition audit for eleven raw-light-curve binary leads.

Public read-only catalog queries preserve exact bytes, URLs, errors and hashes.
Absence from these finite searches does not establish discovery or attribution.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import requests
import xml.etree.ElementTree as ET
import numpy as np
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/stellar_systems/scout_novelty"
RESULTS = ROOT / "results/stellar_systems"
REPORT = ROOT / "reports/stellar_systems/scout_novelty.md"
TICS = [294272014,355544861,38825979,141866704,141411137,
        281581467,388181335,394136010,197764088,52285789,403150111]
VIZIER = {
    "vsx": "B/vsx/vsx", "asassn": "II/366/catv2021",
    "shi2022": "J/ApJS/259/50/table1",
    "tang2025": "J/AJ/170/356/table1",
    "youngster2022": "J/MNRAS/511/4285/table2",
}


def utc():
    return datetime.now(timezone.utc).isoformat()


def fetch(name, url, params=None):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    manifest = OUT / (name + ".manifest.json")
    if manifest.exists():
        metadata = json.loads(manifest.read_text())
        path = Path(metadata.get("local_path",path))
    if path.exists() and manifest.exists():
        if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["sha256"]:
            raise RuntimeError("Snapshot checksum mismatch: " + str(path))
        return path
    if path.exists() or manifest.exists():
        raise RuntimeError("Incomplete snapshot: " + str(path))
    try:
        response = requests.get(url, params=params, timeout=(15, 75),
                    headers={"User-Agent":"Space_ChatGPT-public-stellar-research/1.1"})
        body = response.content
        if len(body)>10_000_000:
            path=Path(os.environ.get("TESS_DATA_DIR","/workspace/tess-research-data"))/"stellar_systems/catalogs"/name
            path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(body)
        metadata = dict(retrieved_utc=utc(), requested_url=url,
            parameters=params, url=response.url, http_status=response.status_code,
            content_type=response.headers.get("Content-Type"), bytes=len(body),
            sha256=hashlib.sha256(body).hexdigest(),local_path=str(path))
        manifest.write_text(json.dumps(metadata, indent=2) + "\n")
        print("FETCH", name, response.status_code, len(body), flush=True)
        return path
    except Exception as exc:
        error = OUT / (name + ".error.json")
        if not error.exists():
            error.write_text(json.dumps(dict(attempted_utc=utc(), url=url,
                parameters=params, error=type(exc).__name__+": "+str(exc)), indent=2)+"\n")
        print("FAILURE", name, type(exc).__name__, str(exc)[:100], flush=True)
        return None


def tap(name, query, service):
    base = {"gaia":"https://gea.esac.esa.int/tap-server/tap/sync",
            "simbad":"https://simbad.cds.unistra.fr/simbad/sim-tap/sync"}[service]
    return fetch(name + ".csv", base,
        dict(request="doQuery", lang="ADQL", format="csv", query=query))


def read_csv(path):
    if not path or not path.exists(): return []
    text = path.read_text()
    if text.lstrip().startswith("<"): return []
    return list(csv.DictReader(io.StringIO(text)))


def valid_path(name):
    """Select a checksum-verified successful immutable retry, preserving errors."""
    original=OUT/name
    candidates=[original]+sorted(OUT.glob(original.stem+"_retry*"+original.suffix))
    for path in reversed(candidates):
        m=OUT/(path.name+".manifest.json")
        if not m.exists():continue
        metadata=json.loads(m.read_text())
        local=Path(metadata.get("local_path",path))
        if metadata.get("http_status")!=200 or not local.exists():continue
        if hashlib.sha256(local.read_bytes()).hexdigest()!=metadata["sha256"]:
            raise RuntimeError("Snapshot hash mismatch: "+str(local))
        if local.suffix==".csv" and local.read_text().lstrip().startswith("<"):continue
        return local
    return None


def read_vizier(path):
    if not path or not path.exists(): return []
    lines = [line for line in path.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    if not lines or "\t" not in lines[0]: return []
    # ASU-TSV gives a header, units row, dash ruler, then scientific rows.
    fields = [x.strip() for x in lines[0].split("\t")]
    rows = []
    for line in lines[1:]:
        values = [x.strip() for x in line.split("\t")]
        if len(values) != len(fields): continue
        if all(not v or set(v) <= {"-"} for v in values): continue
        row = dict(zip(fields, values))
        try:
            float(row.get("_r", ""))
        except ValueError:
            continue
        rows.append(row)
    return rows


def separation(ra, dec, a, d):
    ra,dec,a,d = map(math.radians,(ra,dec,a,d))
    h = math.sin((dec-d)/2)**2 + math.cos(dec)*math.cos(d)*math.sin((ra-a)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1, h))))*3600


def targets():
    rows = list(csv.DictReader((RESULTS/"scout_raw_measurements.csv").open()))
    return [r for tic in TICS for r in rows if int(r["tic"]) == tic]


def download():
    jobs = []
    for row in targets():
        tic = int(row["tic"]); ra = float(row["ra"]); dec = float(row["dec"])
        for label, catalog in VIZIER.items():
            # Old internal filename is preserved; scientific label is corrected.
            filelabel = "cdips2022" if label == "youngster2022" else label
            jobs.append((fetch, (f"tic{tic}_{filelabel}.tsv",
                "https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
                {"-source":catalog, "-c":f"{ra} {dec}", "-c.rs":120,
                 "-out.max":"unlimited", "-out.all":"", "-out.add":"_r,_RAJ2000,_DEJ2000", "-out.form":"tsv"})))
        jobs.append((fetch, (f"tic{tic}_exofop.json",
            f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json")))
        jobs.append((tap, (f"tic{tic}_gaia_neighbors",
            f"select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.ruwe,g.non_single_star,g.ipd_frac_multi_peak,v.in_vari_eclipsing_binary,c.best_class_name,c.best_class_score,e.frequency,e.global_ranking from gaiadr3.gaia_source g left outer join gaiadr3.vari_summary v on g.source_id=v.source_id left outer join gaiadr3.vari_classifier_result c on g.source_id=c.source_id left outer join gaiadr3.vari_eclipsing_binary e on g.source_id=e.source_id where 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333)) order by g.phot_g_mean_mag", "gaia")))
        jobs.append((tap, (f"tic{tic}_simbad",
            f"select oid,main_id,ra,dec,otype,nbref from basic where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333))", "simbad")))
    with ThreadPoolExecutor(max_workers=6) as pool:
        for future in as_completed([pool.submit(fn,*args) for fn,args in jobs]):
            future.result()
    # Preserve aliases/references of actual <=2arcsec SIMBAD counterparts.
    jobs = []
    for row in targets():
        tic = int(row["tic"]); ra = float(row["ra"]); dec = float(row["dec"])
        for s in read_csv(OUT/f"tic{tic}_simbad.csv"):
            if separation(ra,dec,float(s["ra"]),float(s["dec"])) <= 2:
                oid = s["oid"]
                jobs.append((tap, (f"tic{tic}_simbad_aliases_{oid}",
                    f"select id from ident where oidref={oid}", "simbad")))
                jobs.append((tap, (f"tic{tic}_simbad_refs_{oid}",
                    f"select r.bibcode,r.title,r.doi from has_ref h join ref r on h.oidbibref=r.oidbib where h.oidref={oid}", "simbad")))
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(fn,*args) for fn,args in jobs]):
            future.result()


def derive():
    all_records = []
    sources = []
    for m in sorted(OUT.glob("*.manifest.json")):
        metadata = json.loads(m.read_text())
        p = Path(metadata.get("local_path",OUT/m.name.removesuffix(".manifest.json")))
        if hashlib.sha256(p.read_bytes()).hexdigest() != metadata["sha256"]:
            raise RuntimeError("Snapshot hash mismatch: "+str(p))
        sources.append(dict(file=str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p), **metadata))
    for target in targets():
        tic=int(target["tic"]); ra=float(target["ra"]); dec=float(target["dec"])
        record = dict(tic=tic, ra=ra, dec=dec, period_days=float(target["period_days"]),
            catalog_rows={}, gaia_neighbors=[], simbad_counterparts=[])
        for label in VIZIER:
            filelabel = "cdips2022" if label == "youngster2022" else label
            record["catalog_rows"][label] = read_vizier(OUT/f"tic{tic}_{filelabel}.tsv")
        for row in read_csv(valid_path(f"tic{tic}_gaia_neighbors.csv")):
            row["separation_arcsec"] = separation(ra,dec,float(row["ra"]),float(row["dec"]))
            if row.get("frequency"):
                row["eb_period_days"] = 1/float(row["frequency"])
            record["gaia_neighbors"].append(row)
        record["gaia_neighbors"].sort(key=lambda r:r["separation_arcsec"])
        for row in read_csv(OUT/f"tic{tic}_simbad.csv"):
            row["separation_arcsec"] = separation(ra,dec,float(row["ra"]),float(row["dec"]))
            record["simbad_counterparts"].append(row)
        record["simbad_counterparts"].sort(key=lambda r:r["separation_arcsec"])
        all_records.append(record)
        print("TARGET",tic,"CATS",{k:len(v) for k,v in record["catalog_rows"].items()},
              "GAIA EB",[(r["source_id"],round(r["separation_arcsec"],2),r.get("eb_period_days"),r.get("best_class_name")) for r in record["gaia_neighbors"] if r.get("frequency") or r.get("best_class_name")=="ECL"],flush=True)
    result=dict(checked_utc=utc(),method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope="Finite 120arcsec positional cones in VSX, ASAS-SN 2021/2024 snapshot, Shi2022EA, Tang2025eccentric, YOUNGSTER2022, GaiaDR3 EB/classification, and SIMBAD. Position matches do not attribute TESS signals; absence is not proof of novelty.",
        records=all_records,sources=sources,
        failures=[json.loads(p.read_text()) for p in OUT.glob("*.error.json")])
    (RESULTS/"scout_novelty.json").write_text(json.dumps(result,indent=2)+"\n")


def lead_checks():
    jobs=[]
    for row in targets():
        tic=int(row["tic"])
        if tic not in [38825979,197764088,52285789,403150111]:continue
        ra=float(row["ra"]);dec=float(row["dec"])
        exo=json.loads((OUT/f"tic{tic}_exofop.json").read_text())
        aliases=[a.strip() for a in exo["basic_info"].get("star_names","").split(",")]
        for p in OUT.glob(f"tic{tic}_simbad_aliases_*.csv"):
            aliases += [r["id"] for r in read_csv(p)]
        aliases=list(dict.fromkeys(aliases))
        query=" OR ".join('all:"'+a+'"' for a in aliases)
        jobs.append((fetch,(f"tic{tic}_arxiv_aliases.xml",
            "https://export.arxiv.org/api/query",dict(search_query=query,max_results=100))))
        for label,catalog in [("sb9","B/sb9/main"),("wds","B/wds/wds"),
                              ("traven_sb2","J/A+A/638/A145/catalog"),
                              ("galah_dr3","J/MNRAS/506/150")]:
            jobs.append((fetch,(f"tic{tic}_{label}.tsv",
                "https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
                {"-source":catalog,"-c":f"{ra} {dec}","-c.rs":120,
                 "-out.max":"unlimited","-out.all":"","-out.add":"_r,_RAJ2000,_DEJ2000","-out.form":"tsv"})))
        cols="source_id,ra,dec,phot_g_mean_mag,bp_rp,ruwe,non_single_star,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,astrometric_excess_noise,ipd_frac_multi_peak,radial_velocity,radial_velocity_error,rv_nb_transits,rv_amplitude_robust,rv_chisq_pvalue,rv_renormalised_gof,rv_template_teff,rv_expected_sig_to_noise"
        jobs.append((tap,(f"tic{tic}_gaia_context",
            f"select {cols} from gaiadr3.gaia_source where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.0333333333333333)) order by phot_g_mean_mag","gaia")))
    sids="4675527035087606144,6410698972335936512,6410698633034526976"
    for tic in [52285789,403150111]:
        e=json.loads((OUT/f"tic{tic}_exofop.json").read_text())
        m=re.search(r"Gaia DR3 (\d+)",e["basic_info"].get("star_names",""))
        if m:sids+=","+m.group(1)
    for table in ["nss_two_body_orbit","nss_acceleration_astro","nss_non_linear_spectro","nss_vim_fl"]:
        jobs.append((tap,("extended_leads_gaia_"+table,f"select * from gaiadr3.{table} where source_id in ({sids})","gaia")))
    jobs.append((fetch,("traven2020_full.tsv","https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
        {"-source":"J/A+A/638/A145/catalog","-out.all":"","-out.max":"unlimited","-out.form":"tsv"})))
    jobs.append((tap,("tic197764088_neighbor_simbad_refs_7735619",
        "select r.bibcode,r.title,r.doi from has_ref h join ref r on h.oidbibref=r.oidbib where h.oidref=7735619","simbad")))
    jobs.append((tap,("tic197764088_neighbor_simbad_aliases_7735619",
        "select id from ident where oidref=7735619","simbad")))
    query="select table_name,description from TAP_SCHEMA.tables where (description like '%GALAH%' or description like '%APOGEE%' or description like '%TESS-HERMES%') and (description like '%binary%' or description like '%Binary%' or description like '%binaries%' or description like '%Binaries%' or description like '%double-lined%' or description like '%ouble-Lined%' or description like '%double-star%' or description like '%SB2%' or description like '%Joker%')"
    jobs.append((fetch,("spectroscopic_binary_catalogs.csv",
        "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
        dict(request="doQuery",lang="ADQL",format="csv",query=query))))
    with ThreadPoolExecutor(max_workers=5) as pool:
        for future in as_completed([pool.submit(fn,*args) for fn,args in jobs]):future.result()


def primitive(row):
    result={}
    for k in row.colnames:
        v=row[k]
        if np.ma.is_masked(v):v=None
        elif isinstance(v,np.generic):v=v.item()
        result[k]=v
    return result


def full_crossmatch():
    """Match exact identifiers and real sky positions in two complete sources."""
    p=valid_path("traven2020_full.tsv")
    lines=[line for line in p.read_text().splitlines() if line.strip() and not line.startswith("#")]
    binary_rows=[r for r in csv.DictReader(lines,delimiter="\t") if r["GaiaDR2"].strip().isdigit()]
    binary_rows=[{k:v.strip() for k,v in r.items()} for r in binary_rows]
    binary_sky=SkyCoord([float(r["RAJ2000"]) for r in binary_rows]*u.deg,
                        [float(r["DEJ2000"]) for r in binary_rows]*u.deg)
    tables={n:ascii.read(valid_path(f"tess10000_table{n}_mrt.txt"),format="cds") for n in [3,4]}
    if len(tables[3])!=7936 or len(tables[4])!=2065:
        raise RuntimeError("Unexpected TESS Ten Thousand table counts")
    skies={n:SkyCoord(np.asarray(t["RAdeg"],dtype=float)*u.deg,
                     np.asarray(t["DEdeg"],dtype=float)*u.deg) for n,t in tables.items()}
    target_list=[]
    for row in targets():
        tic=int(row["tic"])
        exo=json.loads((OUT/f"tic{tic}_exofop.json").read_text())
        m=re.search(r"Gaia DR2 (\d+)",exo["basic_info"].get("star_names",""))
        target_list.append(dict(tic=tic,ra=float(row["ra"]),dec=float(row["dec"]),
            gaia_dr2_id=m.group(1) if m else None,group="eleven_scouts"))
    # Earlier retained stellar leads: coordinates and identifiers are frozen
    # from scripts/stellar_prior_audit.py, not fitted to catalog detections.
    earlier={300814564:("5273775087036017024",116.92191756103277,-68.62140624817651),
        273879358:("4692911684016155008",32.30851048028379,-71.32001533152403),
        270312735:("4617279547612110208",32.62890012752531,-83.14812279702222),
        370247818:("4625697305555206912",55.45564700638868,-79.11323002464255),
        238202207:("4629516356114834944",58.37129206255529,-73.78568182519616),
        269994642:("6806240691223695744",311.4615173986208,-24.650049533874235)}
    for tic,(sid,ra,dec) in earlier.items():
        target_list.append(dict(tic=tic,ra=ra,dec=dec,gaia_dr2_id=sid,
            gaia_id_scope="Earlier retained Gaia identifier; unchanged DR2/DR3 correspondence is not assumed beyond coordinate matches.",group="earlier_six"))
    supplement=ROOT/"data/stellar_systems/prior_audit/supplement_gaia_sources.csv"
    if supplement.exists():
        sid_to_tic={"4692910966756908288":273879360,"5496090260794096128":260352037,
                    "5496090260793489024":737466840,"5496090256498964096":260352041}
        for row in read_csv(supplement):
            if row["source_id"] in sid_to_tic:
                target_list.append(dict(tic=sid_to_tic[row["source_id"]],ra=float(row["ra"]),dec=float(row["dec"]),
                    gaia_id=row["source_id"],gaia_id_scope="Gaia DR3; SB2 exact identifier comparison is supplemented by position.",group="hidden_source_supplement"))
    records=[]
    for target in target_list:
        center=SkyCoord(target["ra"]*u.deg,target["dec"]*u.deg)
        distances=center.separation(binary_sky).arcsec
        sb2=[]
        for i in np.flatnonzero(distances<=120):
            r=dict(binary_rows[i]);r["separation_arcsec"]=float(distances[i])
            r["exact_gaia_identifier"]=r["GaiaDR2"]==target.get("gaia_dr2_id",target.get("gaia_id"))
            sb2.append(r)
        catalog_matches={}
        for n,t in tables.items():
            distances=center.separation(skies[n]).arcsec
            matches=[]
            for i in np.flatnonzero((distances<=120)|(t["TIC"]==target["tic"])):
                r=primitive(t[i]);r["separation_arcsec"]=float(distances[i]);r["exact_tic"]=r["TIC"]==target["tic"]
                matches.append(r)
            catalog_matches[f"table{n}"]=matches
        records.append(dict(**target,traven2020_sb2_matches=sb2,tess10000_matches=catalog_matches))
    result=dict(checked_utc=utc(),method_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        table_counts={"traven2020_successfully_analyzed_sb2_spectra":len(binary_rows),
                      "kostov2025_table3_new_ebs":len(tables[3]),"kostov2025_table4_previously_known_ebs":len(tables[4])},
        scope="Exact TIC and <=120arcsec positions in both vetted TESS Ten Thousand tables; complete successfully analyzed GALAH SB2 catalog. Nearby matches do not attribute TESS signals. Excludes unvetted TESS neural table2; absence does not establish novelty.",
        records=records,
        sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/"traven2020_full.tsv",OUT/"tess10000_table3_mrt.txt",OUT/"tess10000_table4_mrt.txt"]})
    (RESULTS/"scout_fullcatalog_matches.json").write_text(json.dumps(result,indent=2)+"\n")
    print("VERIFIED COMPLETE CATALOG COUNTS",result["table_counts"],flush=True)
    for r in records:
        print("COMPLETE MATCH",r["tic"],"SB2",[(x["GaiaDR2"],round(x["separation_arcsec"],2)) for x in r["traven2020_sb2_matches"]],"TESS",{k:[(x["TIC"],round(x["separation_arcsec"],2)) for x in v] for k,v in r["tess10000_matches"].items()},flush=True)


def report():
    provenance=json.loads((RESULTS/"scout_novelty.json").read_text())
    complete=json.loads((RESULTS/"scout_fullcatalog_matches.json").read_text())
    decisions={
        294272014:("known_nearby_eb","Gaia DR3 5486571994791177344 / VSX3356349",2.33,.260680,"Known eclipsing star 2.33arcsec away at approximately half the detected period; TESS source attribution remains unresolved."),
        355544861:("known_nearby_eb","Gaia DR3 6463554115486620544 / VSX7940975",11.39,.492630,"Known nearby eclipsing star with the same period; probable contaminant, requiring pixels for attribution."),
        38825979:("known_host_sb2","GALAH160124002601234 / Gaia4675527035087606144",.038,None,"Known double-lined spectroscopic binary from Traven2020. TESS eclipsing characterization may add information, but a second star is already recognized."),
        141866704:("known_nearby_eb","TIC141866705 / ASASSN-V J061547.24-751510.2",79.79,1.4232,"Known eclipsing neighbor with same or twice-period ambiguity; candidate target CTOI is a separate existing recognition."),
        141411137:("known_host_eb","Gaia DR3 4650316195543230592 / VSX4313616",.17,1.757750,"Host is already a recognized eclipsing binary in Gaia and VSX."),
        281581467:("known_host_eb","ASASSN-V J074356.39-602347.2 / VSX838559",.02,3.990950,"Host classified EA+BY in VSX; catalog binary period approximately twice the scout's detected period."),
        388181335:("known_nearby_eb","TIC388181338 / ASASSN-V J070310.13-735512.1",53.37,1.5886,"Known eclipsing neighbor with matching period; probable contaminant pending pixels."),
        394136010:("known_host_eb","UCAC4 096-000233 / PATHOS-I2019",.148,None,"SIMBAD positively recognizes host as EB* from Nardiello2019, although several other EB catalogs have no match."),
        197764088:("known_nearby_eb","TIC197764089 / TYC8820-40-1",63.4730487,2.682147,"Kostov2025 vetted EB nearby has matching period and41/16ppt depths. Host Gaia RV is stable, whereas the neighbor is strongly variable; likely known-star contamination."),
        52285789:("known_host_sb2","GALAH170911004701327 / Gaia4692405049673289472",.065,None,"Known double-lined spectroscopic binary from Traven2020; companion is already recognized."),
        403150111:("known_visual_pair_unclassified_eclipsing_subsystem","UCAC4 156-214940 / Gaia6411065522024603392+6411065526320580736",None,None,"Existing SIMBAD double-star recognition and two Gaia sources about1.6arcsec apart. One component has strong RV variability. A previously unclassified eclipsing subsystem remains possible; a new third star is not established."),
    }
    rows=[]
    for r in provenance["records"]:
        tic=r["tic"];status,source,offset,period,assessment=decisions[tic]
        row=dict(tic=tic,scout_period_days=r["period_days"],recognition_status=status,
            recognized_source=source,source_separation_arcsec=offset,
            known_period_days=period,assessment=assessment)
        rows.append(row)
    # Real positive controls prevent treating missing fields as catalog absence.
    r=next(r for r in provenance["records"] if r["tic"]==141411137)
    if not any(v["Type"]=="E" and v["OID"]=="4313616" for v in r["catalog_rows"]["vsx"]):
        raise RuntimeError("Known VSX host positive control failed")
    r=next(r for r in complete["records"] if r["tic"]==38825979)
    if not any(v["GaiaDR2"]=="4675527035087606144" for v in r["traven2020_sb2_matches"]):
        raise RuntimeError("Known SB2 positive control failed")
    r=next(r for r in complete["records"] if r["tic"]==197764088)
    if not any(v["TIC"]==197764089 for v in r["tess10000_matches"]["table3"]):
        raise RuntimeError("Known neighboring EB positive control failed")
    out=io.StringIO();writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (RESULTS/"scout_novelty.csv").write_text(out.getvalue())
    table=["| Scout TIC | Scout period (d) | Prior recognition | Assessment |",
           "| --- | ---: | --- | --- |"]
    for r in rows:
        table.append(f"| {r['tic']} | {r['scout_period_days']:.6f} | {r['recognized_source']} | {r['assessment']} |")
    body=f"""# Stellar scout recognition audit

Checked {provenance['checked_utc']} (UTC). Eleven selected repeatable eclipse-like scout signals were investigated against multiple source catalogs. **Ten have a positively recognized binary host or a nearby known eclipsing binary with a matching period/harmonic. The remaining target, TIC403150111, is already a recognized visual double, with an eclipsing subsystem not yet fully classified. These results do not establish a previously unknown star.** Nearby matches are contamination hypotheses, not pixel-based source assignments.

The important falsifications were TIC38825979 and TIC52285789: they had no match in the original eclipsing-binary checks, but both are published double-lined spectroscopic binaries in [Traven et al.2020](https://doi.org/10.1051/0004-6361/202037484). TIC197764088 has a matching known EB, TIC197764089, 63.473arcsec away in [Kostov et al.2025](https://doi.org/10.3847/1538-4365/ade2d8). TIC394136010 is already recognized from [Nardiello et al.2019 PATHOS-I](https://doi.org/10.1093/mnras/stz2878). These cases show why absence from one EB catalog cannot establish a new companion.

{chr(10).join(table)}

## Catalog scope and interpretation

The checks include 120arcsec cones in current VizieR VSX, ASAS-SN II/366/catv2021 (a2021 catalog version archived2024), Shi et al.2022 EA binaries, Tang et al.2025 eccentric binaries, Battley et al.2022 **YOUNGSTER** EA/EB catalog, Gaia DR3 EB/variability classification and SIMBAD. The earlier internal filenames `cdips2022.tsv` are preserved unchanged; J/MNRAS/511/4285/table2 is YOUNGSTER, not CDIPS. Exact bytes, URLs, parameters, timestamps and SHA256 hashes are retained alongside every response. One initial Gaia HTTP503 for TIC281581467 was preserved and a successful immutable retry reused.

Additional checks for TIC38825979,197764088,52285789,403150111 include SB9, WDS, GALAH DR3, GALAH double-lined binaries, all four Gaia DR3 NSS tables and finite arXiv alias metadata searches. Zero results from arXiv metadata searches do not exclude mentions in full text or supplementary tables. NSS=0 does not imply a single star: TIC38825979 has Gaia RV amplitude174.9km/s while its binary is positively recognized in GALAH.

The complete GALAH Traven2020 table has12,760 successfully analyzed entries. The complete **TESS Ten Thousand** publication comprises **7,936 new EBs in Table3 and2,065 previously known EBs in Table4**, verified using Astropy's CDS/MRT parser. Both were checked by exact TIC and <=120arcsec separation; they must not be confused with Table2, which contains about900,000 unvetted/unvalidated neural-network detections. Position matches alone do not establish that an aperture's eclipses originate on the nearby source.

TIC403150111 has two Gaia sources about1.6arcsec apart. Both have large IPD multi-peak fractions (47%,33%), and one shows RV amplitude119.7km/s with a zero rounded chi-square p-value while the other has a5.37km/s amplitude and p=.982. This favors investigating an unresolved eclipsing subsystem, but Gaia summary RV scatter is not a solved orbital model, masses, or proof of a third star. Existing visual-pair recognition must remain explicit.

## Retained hidden-source leads

The separate full-table audit includes the six earlier stellar leads and proposed sources TIC273879360,260352037,737466840,260352041. They have no exact or120arcsec-position match in the complete Traven2020 successfully analyzed SB2 table or either vetted TESS Ten Thousand table. This is a bounded negative recognition check. It does not establish novelty, and any presence in unvetted Table2 must be recorded separately as a prior pipeline flag. Source attribution, stellar classification and further literature checks remain necessary.

## Reproduction and evidence

Run `/workspace/openrouter-check-venv/bin/python scripts/stellar_scout_novelty.py --download --lead-checks --full-crossmatch --report` from the repository. Existing successful snapshots are hash-verified and reused. [Searchable assessments](../../results/stellar_systems/scout_novelty.csv), [full target/cone rows and provenance](../../results/stellar_systems/scout_novelty.json), and [complete-table matches](../../results/stellar_systems/scout_fullcatalog_matches.json) retain the scientific scope. The two TESS machine-readable supplements and complete Traven table are preserved under `data/stellar_systems/scout_novelty/`, with immutable request manifests. No authenticated or paid API was used.
"""
    REPORT.parent.mkdir(parents=True,exist_ok=True);REPORT.write_text(body)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download",action="store_true")
    parser.add_argument("--lead-checks",action="store_true")
    parser.add_argument("--full-crossmatch",action="store_true")
    parser.add_argument("--report",action="store_true")
    args=parser.parse_args()
    if args.download:download()
    if args.lead_checks:lead_checks()
    derive()
    if args.full_crossmatch:full_crossmatch()
    if args.report:report()
