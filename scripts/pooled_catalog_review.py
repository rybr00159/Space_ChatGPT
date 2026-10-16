"""Lead-only, finite recognition for the completed pooled sector search.

Calls the frozen recognition module without changing old evidence or outputs.
The target and nearby-source identities are kept distinct. Catalog absence is
never a claim of novelty. All HTTP bodies and manifests remain in timestamped
live/deep directories belonging to this run.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from urllib.parse import urlencode

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import multisector_recognition as recognition

TICS = [382190230, 353807236, 359679661, 375090561, 353933809,
        366490373, 293268396, 287086419, 349972227, 356821129]
PERIODS = [1.282392, 6.224718, 1.913410, 8.136209, 18.978402,
           2.810239, 2.596856, 12.613388, 8.276215, 3.871952]
OUT = ROOT / "results/multisector_recognition/pooled_context"


def separation(ra, dec, row):
    r2 = row.get("ra", row.get("RAJ2000"))
    d2 = row.get("dec", row.get("DEJ2000"))
    if not r2 or not d2:
        return None
    r1, d1, r2, d2 = map(math.radians, [ra, dec, float(r2), float(d2)])
    return math.degrees(2 * math.asin(min(1, math.sqrt(
        math.sin((d2-d1)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((r2-r1)/2)**2))))*3600


def controls():
    prior = ROOT / "data/new_stars/rv_recognition/4690783166942834816_vsx.tsv"
    rows = recognition.parse_vsx_tsv(prior.read_bytes())
    assert any(row["Name"] == "ASASSN-V J005926.01-712933.7" for row in rows)
    failed = False
    try:
        recognition.parse_vsx_tsv(b"#INFO\tError=Table or Catalog not found\n")
    except ValueError:
        failed = True
    assert failed
    index = recognition.FreshCatalogIndex()
    assert index.match_signal(4206066, 3.182785)["known_period_match"]
    assert index.match_signal(71316387, 2*.6518369465110034)["known_harmonic_match"]
    result = {"known_rotator_recovered": True, "service_error_rejected": True,
              "reported_manuscript_recovered": True, "raven_harmonic_recovered": True,
              "scope": "Retained positive inputs, not simultaneous remote positive requests."}
    recognition.write_json(OUT / "positive_controls.json", result)


def fetch():
    OUT.mkdir(parents=True, exist_ok=True)
    recognition.RESULTS = OUT
    with ThreadPoolExecutor(max_workers=4) as executor:
        live = list(executor.map(lambda tic: str(recognition.capture_live_context([tic])), TICS))
    recognition.write_json(OUT / "run_sources.json", {"live": live, "tics": TICS})
    deep = str(recognition.capture_deep_context(TICS))
    recognition.write_json(OUT / "run_sources.json", {"live": live, "deep": deep, "tics": TICS})


def summarize():
    sources = json.loads((OUT / "run_sources.json").read_text())
    deep = {int(row["tic"]): row for row in json.loads((Path(sources["deep"]) / "summary.json").read_text())["checks"]}
    index = recognition.FreshCatalogIndex()
    records = []
    for live_path in sources["live"]:
        directory = Path(live_path)
        for live in json.loads((directory / "summary.json").read_text())["checks"]:
            tic = int(live["tic"])
            target_path = directory / f"tic{tic}_exofop.json"
            body = target_path.read_bytes()
            manifest = json.loads(target_path.with_suffix(".json.manifest.json").read_text())
            assert manifest["valid"] and recognition.sha(body) == manifest["sha256"]
            target = json.loads(body)
            record = {"tic": tic, "live": live, "deep": deep[tic],
                      "target_parameters": {k:v for k,v in target.items() if k not in ("files",)},
                      "exofop_source": str(target_path.relative_to(ROOT)),
                      "exofop_sha256": recognition.sha(body)}
            ra, dec = float(target["coordinates"]["ra"]), float(target["coordinates"]["dec"])
            aliases = target["basic_info"].get("star_names", "")
            record["identity"] = {key: re.findall(pattern, aliases) for key,pattern in [
                ("gaia_dr3", r"Gaia DR3 (\d+)"), ("gaia_dr2", r"Gaia DR2 (\d+)"),
                ("twomass", r"2MASS ([^,]+)")]}
            record["recognition"] = index.match_signal(tic, PERIODS[TICS.index(tic)], ra, dec)
            for key in ("vsx_rows", "simbad_rows"):
                for row in record["live"].get(key, []) or []:
                    row["separation_arcsec"] = separation(ra, dec, row)
            for key in ("asassn", "gaia_neighbors", "simbad_bibliography", "simbad_aliases"):
                for row in record["deep"].get(key, []) or []:
                    row["separation_arcsec"] = separation(ra, dec, row)
                    row["association_scope"] = "same_position_within_2arcsec" if row["separation_arcsec"] is not None and row["separation_arcsec"] < 2 else "nearby_source_do_not_assign_to_target"
            records.append(record)
    result = {"created_utc": recognition.utc(), "sources": sources, "checks": records,
              "method_sha256": recognition.sha(Path(__file__).read_bytes()),
              "frozen_recognition_sha256": recognition.sha((ROOT / "src/multisector_recognition.py").read_bytes()),
              "limitations": "Finite public ExoFOP target JSON, bulk prior reports, Gaia DR3 classifiers/NSS, VSX 30arcsec, ASAS-SN 2021 60arcsec, SIMBAD 60arcsec bibliography/aliases and exact-TIC arXiv. Catalog nulls and out-of-selection lists provide no novelty proof. Nearby variables require source localization. Neither independent all-TCE inventory nor exhaustive literature priority is established."}
    recognition.write_json(OUT / "review.json", result)
    print(json.dumps({"tics": len(records), "output": str(OUT / "review.json")}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("fetch", "summarize", "controls", "remote-controls", "supplement", "dvxml", "dvsummary", "rotation", "aliases", "extra238", "tce238", "report", "all"))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.command in ("controls", "all"):
        controls()
    if args.command in ("fetch", "all"):
        fetch()
    if args.command in ("summarize", "all"):
        summarize()
    if args.command == "supplement":
        supplemental()
    if args.command == "dvxml":
        dvxml()
    if args.command == "dvsummary":
        dvsummary()
    if args.command == "rotation":
        rotation()
    if args.command == "aliases":
        aliases()
    if args.command == "extra238":
        extra238()
    if args.command == "tce238":
        supplemental([238182368], OUT / "extra_tic238182368")
    if args.command == "remote-controls":
        remote_controls()
    if args.command == "report":
        report()


def supplemental(tics=None, destination=None):
    """Finite CAOM SPOC/DV inventory; no all-TCE completeness implied."""
    directory = recognition.BASE / "deep" / (recognition.utc().replace("-", "").replace(":", "").replace(".", "") + "_pooled_supplement")
    directory.mkdir(parents=True, exist_ok=False)
    endpoint = "https://mast.stsci.edu/api/v0/invoke"

    def mast(name, service, params):
        request = dict(service=service, params=params, format="json", pagesize=100000, page=1)
        path = directory / (name + ".json")
        metadata = {"url": endpoint, "request": request, "retrieved_utc": recognition.utc(), "valid": False}
        try:
            response = requests.post(endpoint, data={"request": json.dumps(request)}, timeout=(15, 120))
            path.write_bytes(response.content)
            metadata.update(http_status=response.status_code, bytes=len(response.content), sha256=recognition.sha(response.content))
            response.raise_for_status()
            body = response.json()
            if body.get("status") != "COMPLETE" or "data" not in body or len(body["data"]) >= 100000:
                raise ValueError("Incomplete/paginated MAST response")
            metadata["valid"] = True
            return body["data"]
        except Exception as exc:
            metadata["error"] = str(exc)
            return None
        finally:
            recognition.write_json(path.with_suffix(".json.manifest.json"), metadata)

    ids = (TICS if tics is None else tics) + [55525572]
    observations = mast("spoc_observations", "Mast.Caom.Filtered", {
        "columns": "obsid,obs_id,target_name,sequence_number,dataproduct_type,provenance_name",
        "filters": [{"paramName": "obs_collection", "values": ["TESS"]},
                    {"paramName": "provenance_name", "values": ["SPOC"]},
                    {"paramName": "target_name", "values": [str(x) for x in ids]}]})
    products = []
    failures = []
    if observations:
        obsids = sorted({str(row["obsid"]) for row in observations})
        for start in range(0, len(obsids), 30):
            batch = mast(f"spoc_products_{start//30:03d}", "Mast.Caom.Products", {"obsid": ",".join(obsids[start:start+30])})
            if batch is None:
                failures.append(start)
            else:
                products.extend(batch)
    result = {"created_utc": recognition.utc(), "directory": str(directory.relative_to(ROOT)),
              "observations": observations, "product_batches_failed": failures,
              "dv_products": [row for row in products if row.get("productSubGroupDescription", "").startswith("DV") or "_dvt." in row.get("productFilename", "")],
              "positive_control_tic": 55525572,
              "limitations": "Finite target-name CAOM SPOC observations and associated products only. Products are prior public data, not a complete inventory of every pipeline TCE. Lack of these products or target pages does not establish novelty."}
    recognition.write_json((OUT if destination is None else destination) / "spoc_dv_inventory.json", result)
    print(json.dumps({"observations": None if observations is None else len(observations), "products": len(products), "dv_products": len(result["dv_products"])}), flush=True)


def dvxml():
    """Read numerical TCE parameters from actual public SPOC DV reports."""
    import xml.etree.ElementTree as ET
    inventory = json.loads((OUT / "spoc_dv_inventory.json").read_text())
    directory = recognition.BASE / "deep" / (recognition.utc().replace("-", "").replace(":", "").replace(".", "") + "_pooled_dvxml")
    directory.mkdir(parents=True, exist_ok=False)
    products = {p["productFilename"]: p for p in inventory["dv_products"] if p["productFilename"].endswith("_dvr.xml")}
    # The control's longest baseline provides a positive known-planet recovery.
    controls = [p for p in products.values() if "0000000055525572" in p["productFilename"]]
    chosen = [p for p in products.values() if "0000000055525572" not in p["productFilename"]]
    if controls:
        chosen.append(sorted(controls, key=lambda x: x["productFilename"])[-1])

    def task(product):
        def validator(body):
            root = ET.fromstring(body)
            if "dvTargetResults" not in root.tag:
                raise ValueError("Unexpected SPOC DV XML root")
            return parse_dv_xml(root)
        url = "https://mast.stsci.edu/api/v0.1/Download/file?" + urlencode({"uri": product["dataURI"]})
        parsed = recognition.retain_public_response(directory, product["productFilename"], url, validator)
        return {"product": product, "parsed": parsed, "valid": parsed is not None}
    with ThreadPoolExecutor(max_workers=3) as executor:
        checks = list(executor.map(task, chosen))
    recognition.write_json(OUT / "spoc_dv_xml.json", {"directory": str(directory.relative_to(ROOT)), "checks": checks,
        "limitations": "Public SPOC TCEs are prior detections, not validated planets; finite product inventory and pipeline versions retained. Exact period identity and source localization must be assessed separately."})
    print(json.dumps({"xml_reports": len(checks), "valid": sum(c["valid"] for c in checks)}), flush=True)


def parse_dv_xml(root):
    ns = {"d": root.tag.split("}")[0][1:]}
    def convert(element):
        return {"tag": element.tag.split("}")[-1], "attributes": element.attrib,
                "children": [convert(child) for child in element]}
    planets = []
    for planet in root.findall("d:planetResults", ns):
        result = {"planet_number": planet.attrib["planetNumber"]}
        for tag in ["planetCandidate", "binaryDiscriminationResults", "centroidResults", "secondaryEventResults", "bootstrapResults"]:
            item = planet.find("d:" + tag, ns)
            if item is not None:
                # Bootstrap distribution arrays are kept in original XML only.
                result[tag] = dict(item.attrib) if tag in ("planetCandidate", "bootstrapResults") else convert(item)
        for tag in ["allTransitsFit", "oddTransitsFit", "evenTransitsFit", "trapezoidalFit"]:
            fit = planet.find("d:" + tag, ns)
            if fit is not None:
                result[tag] = {"attributes": fit.attrib, "parameters": {
                    item.attrib["name"]: item.attrib for item in fit.findall(".//d:modelParameter", ns)}}
        planets.append(result)
    return {"root_attributes": root.attrib, "planets": planets}


def dvsummary():
    """Re-derive compact summary from retained, checksum-verified source XML."""
    import xml.etree.ElementTree as ET
    source = OUT / "spoc_dv_xml.json"
    original = json.loads(source.read_text())
    directory = ROOT / original["directory"]
    for check in original["checks"]:
        path = directory / check["product"]["productFilename"]
        body = path.read_bytes()
        manifest = json.loads(path.with_suffix(".xml.manifest.json").read_text())
        assert manifest["valid"] and recognition.sha(body) == manifest["sha256"]
        check["parsed"] = parse_dv_xml(ET.fromstring(body))
        check["raw_sha256"] = manifest["sha256"]
    source.write_text(json.dumps(original, indent=2) + "\n")
    known = [check for check in original["checks"] if "0000000055525572" in check["product"]["productFilename"]]
    assert any(abs(float(p["allTransitsFit"]["parameters"]["orbitalPeriodDays"]["value"]) - 83.8911) < .02
               for check in known for p in check["parsed"]["planets"])
    print("DV XML positive known-planet period control passed; compact summary regenerated.")


def rotation():
    """Actual published classification rows behind a target bibliography hit."""
    directory = recognition.BASE / "deep" / (recognition.utc().replace("-", "").replace(":", "").replace(".", "") + "_pooled_rotation")
    directory.mkdir(parents=True, exist_ok=False)
    result = {"directory": str(directory.relative_to(ROOT)), "catalog": "J/ApJ/930/7; Avallone et al.2022 ApJ9307"}
    for table, required in [("table1", ("Mission", "ID", "Class")), ("table2", ("Mission", "ID", "PRot"))]:
        name = "J/ApJ/930/7/" + table
        url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source": name, "-out.all": "", "-out.max": "unlimited"})
        rows = recognition.retain_public_response(directory, table + ".tsv", url, lambda b, n=name, r=required: recognition.parse_asu_tsv(b, n, r))
        if rows is not None:
            if table == "table1":
                assert any(row["ID"] == "25063396" and row["Class"] == "RUWE Binary" for row in rows)
            result[table] = {"valid": True, "row_count": len(rows), "target_rows": [row for row in rows if row["Mission"] == "TESS" and int(row["ID"]) in TICS]}
        else:
            result[table] = {"valid": False}
    result["scope"] = "Published 2022 finite sample and classification labels, not independently confirmed multiplicity/planet classification. Positive named RUWE Binary control validates table1 identity."
    recognition.write_json(OUT / "rotation_literature.json", result)
    print(json.dumps(result), flush=True)


def aliases():
    directory = recognition.BASE / "deep" / (recognition.utc().replace("-", "").replace(":", "").replace(".", "") + "_pooled_aliases")
    directory.mkdir(parents=True, exist_ok=False)
    review = json.loads((OUT / "review.json").read_text())
    ids = [record["identity"]["gaia_dr3"][0] for record in review["checks"]]
    query = "select dr3_source_id,dr2_source_id,angular_distance,magnitude_difference from gaiadr3.dr2_neighbourhood where dr3_source_id in (" + ",".join(ids) + ")"
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    rows = recognition.retain_public_response(directory, "gaia_dr2_neighborhood.csv", url,
        lambda b: recognition.parse_csv(b, required_columns=("dr3_source_id", "dr2_source_id", "angular_distance"), min_rows=0)[0])
    result = {"directory": str(directory.relative_to(ROOT)), "dr2_neighborhood": rows, "dr2_valid": rows is not None,
              "aliases": [{"tic": record["tic"], "exofop_aliases": record["identity"]} for record in review["checks"]],
              "scope": "Official Gaia DR3-DR2 candidate neighborhood mapping, not every row automatically the same physical source; use angular distance and ambiguity. ExoFOP aliases retained separately."}
    for record in review["checks"]:
        if record["live"].get("vsx_valid"):
            continue
        tic = record["tic"]
        manifest = json.loads((ROOT / record["exofop_source"]).with_name(f"tic{tic}_vsx.tsv.manifest.json").read_text())
        retried = recognition.retain_public_response(directory, f"tic{tic}_vsx_retry.tsv", manifest["url"], recognition.parse_vsx_tsv)
        result.setdefault("vsx_retries", []).append({"tic": tic, "valid": retried is not None, "rows": retried})
    recognition.write_json(OUT / "alias_resolution.json", result)
    print(json.dumps(result), flush=True)


def extra238():
    """Requested companion vetting team's additional mixed photometric lead."""
    tic = 238182368
    recognition.RESULTS = OUT / "extra_tic238182368"
    recognition.RESULTS.mkdir(parents=True, exist_ok=True)
    live = recognition.capture_live_context([tic])
    deep = recognition.capture_deep_context([tic])
    recognition.write_json(recognition.RESULTS / "sources.json", {
        "live": str(live.relative_to(ROOT)), "deep": str(deep.relative_to(ROOT)),
        "recognition": recognition.FreshCatalogIndex().match_signal(tic, 14.444722),
        "scope": "Finite extra lead context, added after independent photometric review; not a discovery or exhaustive literature check."})
    print(json.dumps({"extra_tic": tic, "live": str(live), "deep": str(deep)}), flush=True)


def remote_controls():
    import xml.etree.ElementTree as ET
    directory = recognition.BASE / "deep" / (recognition.utc().replace("-", "").replace(":", "").replace(".", "") + "_pooled_controls")
    directory.mkdir(parents=True, exist_ok=False)
    columns = ("ASASSN-V", "RAJ2000", "DEJ2000", "Vmag", "Amp", "Per", "Type", "Prob", "TIC", "GaiaDR3")
    url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source": "II/366/catv2021", "-c": "45.09515 0.47689", "-c.rs": "10", "-out": ",".join(columns), "-out.max": "unlimited"})
    def validate_asassn(body):
        rows = recognition.parse_asu_tsv(body, "II/366/catv2021", columns)
        if not any(row["ASASSN-V"] == "J030022.84+002836.8" and row["Type"] == "VAR" for row in rows):
            raise ValueError("Known published ASAS-SN 2021 positive object was not recovered")
        return rows
    def validate_arxiv(body):
        root = ET.fromstring(body)
        ns = {"a": "http://www.w3.org/2005/Atom"}
        rows = [{"id": e.findtext("a:id", namespaces=ns), "title": e.findtext("a:title", namespaces=ns)} for e in root.findall("a:entry", ns)]
        assert any("2603.22597" in (r["id"] or "") for r in rows)
        return rows
    jobs = [("asassn_positive.tsv", url, validate_asassn), ("arxiv_positive.xml", "https://export.arxiv.org/api/query?id_list=2603.22597", validate_arxiv)]
    result = {"directory": str(directory.relative_to(ROOT)), "checks": []}
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [(name, executor.submit(recognition.retain_public_response, directory, name, url, validator)) for name, url, validator in jobs]
        for name, future in futures:
            parsed = future.result()
            result["checks"].append({"name": name, "valid": parsed is not None, "rows": parsed})
    recognition.write_json(OUT / "remote_positive_controls.json", result)
    print(json.dumps(result), flush=True)


def report():
    import csv
    review = json.loads((OUT / "review.json").read_text())
    aliases_result = json.loads((OUT / "alias_resolution.json").read_text())
    vsx_retry = {r["tic"]: r for r in aliases_result.get("vsx_retries", [])}
    notes = {
        382190230: "Published Avallone2022 RUWE Binary, vsini75.60km/s; public SPOC TCEs near0.8427d. P1.2824d is not directly identical to those TCEs; variable/binary alternatives need photometric tests.",
        353807236: "No checked prior planet/variable period; no novelty conclusion. Gaia astrometry relatively clean.",
        359679661: "High RUWE and significant Gaia RV variability; no known period identity in finite checked sources.",
        375090561: "Yu2019 reported5.423940d; SPOCfit5.42428305d. Pooled8.136209d is near3:2 alias. GaiaAcceleration9 and TICcontamination1.3083.",
        353933809: "Published GaiaAstroSpectroSB1 P1124.79934+/-38.11657d, highRUWE, nearby comoving source; does not by itself identify the P18.9784d dips.",
        366490373: "Very highRUWE17.33 and significant GaiaRVvariability; known multiplicity warning, no orbital period from checkedNSS tables.",
        293268396: "GaiaRVrobustrange81.13km/s,p0; bright neighbor47arcsec. RV warning is not a calibrated companion mass or identified eclipse source.",
        287086419: "TICcontamination.7654; adjacent catalogued redgiants have literature, not target assignments. FirstVSX503failure retained; retryvalidzero.",
        349972227: "Same-position VSXRSvariable, GaiaRSclassifier score.6723, GaiaRVrobustrange150.29km/s,p0. Priorvariable, not novel star recognition.",
        356821129: "TICcontamination.1602; no targetprior planet/variable period in finite checked sources; nearby SIMBAD bibliography belongs to neighbor, not target."
    }
    rows = []
    for record in review["checks"]:
        tic = record["tic"]
        identity = record["identity"]
        sid = identity["gaia_dr3"][0]
        g = next(row for row in record["deep"]["gaia_neighbors"] if row["source_id"] == sid)
        target = record["target_parameters"]
        params = next(row for row in target["stellar_parameters"] if "srad" in row)
        vsx = record["live"].get("vsx_rows")
        if tic in vsx_retry and vsx_retry[tic]["valid"]:
            vsx = vsx_retry[tic]["rows"]
        rows.append({"tic": tic, "searched_period_days": PERIODS[TICS.index(tic)],
                     "ra": record["deep"]["ra"], "dec": record["deep"]["dec"],
                     "gaia_dr3": sid, "gaia_dr2": ";".join(identity["gaia_dr2"]),
                     "twomass": ";".join(identity["twomass"]),
                     "tic_teff_K": params["teff"], "tic_radius_Rsun": params["srad"],
                     "tic_mass_Msun": params["mass"], "tic_contamination_ratio": target["basic_info"]["tic_contamination_ratio"],
                     "gaia_ruwe": g["ruwe"], "gaia_nss_flag": g["non_single_star"],
                     "gaia_rv_amplitude_robust_kms": g["rv_amplitude_robust"],
                     "gaia_rv_chisq_pvalue": g["rv_chisq_pvalue"],
                     "gaia_best_class": g["best_class_name"], "gaia_best_class_score": g["best_class_score"],
                     "vsx_entries_within30arcsec": "failed" if vsx is None else len(vsx),
                     "recognition_status": record["recognition"]["recognition_status"], "assessment": notes[tic]})
    with (OUT / "target_context.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    recognition.write_json(OUT / "assessments.json", {"created_utc": recognition.utc(), "targets": rows,
        "scope": review["limitations"], "method_sha256": recognition.sha(Path(__file__).read_bytes())})
    print(json.dumps({"target_rows": len(rows), "output": str(OUT / "target_context.csv")}))


if __name__ == "__main__":
    main()
