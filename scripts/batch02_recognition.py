"""Finite current recognition and correlated-width follow-up for batch b02.

The producer-error and empirical-error streams remain separately identified.
Target coordinates join on the exact TIC-sector pair. The two filter widths
share their data and are not independent evidence or false-alarm calibration.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import sys
from urllib.parse import urlencode

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import multisector_recognition as m

INPUT = ROOT / "results/multisector_batches/b02_096_105"
OUT = ROOT / "results/multisector_recognition/b02_096_105"


def annotate():
    OUT.mkdir(parents=True, exist_ok=True)
    main = pd.read_csv(INPUT / "signals.csv")
    recovery = pd.read_csv(INPUT / "error_recovery/signals.csv")
    targets = pd.read_csv(INPUT / "targets.csv")
    assert not targets.duplicated(["tic", "sector"]).any()
    recovery = recovery.merge(targets[["tic", "sector", "ra", "dec"]], on=["tic", "sector"], how="left", validate="many_to_one")
    assert recovery[["ra", "dec"]].notna().all().all()
    main["uncertainty_stream"] = "producer_reported"
    main["error_assumption"] = "producer_reported_uncertainties"
    recovery["uncertainty_stream"] = "empirical_adjacent_difference_MAD"
    recovery["label"] = "unmasked"
    frame = pd.concat([main, recovery], ignore_index=True)
    assert not frame.duplicated(["tic", "sector", "width_days", "label", "rank", "uncertainty_stream"]).any()
    ix = m.FreshCatalogIndex()
    assert ix.match_signal(4206066, 3.182785)["known_period_match"]
    assert ix.match_signal(71316387, 2*.6518369465110034)["known_harmonic_match"]
    evidence = []
    columns = ["known_period_match", "known_harmonic_match", "possible_nearby_period_contamination", "recognition_status", "is_confirmed_host", "is_toi_host", "is_ctoi_host", "is_eb_host", "is_raven_host"]
    annotations = []
    for row in frame.itertuples():
        result = ix.match_signal(int(row.tic), float(row.period_days), float(row.ra), float(row.dec))
        annotations.append({key: result.get(key) for key in columns})
        evidence.append({"tic": int(row.tic), "sector": int(row.sector), "width_days": float(row.width_days), "rank": int(row.rank), "uncertainty_stream": row.uncertainty_stream, "recognition": result})
    frame = pd.concat([frame, pd.DataFrame(annotations)], axis=1)
    frame.to_csv(OUT / "signals_annotated.csv.gz", index=False)
    with gzip.open(OUT / "recognition_evidence.json.gz", "wt") as stream:
        json.dump(evidence, stream, allow_nan=False)
    valid = frame[(frame.period_days > .8) & (frame.local_depth_ppm >= 100) & (frame.local_depth_ppm <= 15000) & (frame.n_events >= 3) & (frame.local_snr >= 6)].copy()
    queue = []
    for _, group in valid.groupby(["tic", "sector", "label", "uncertainty_stream"]):
        for _, row in group.iterrows():
            match = group[(group.period_days / row.period_days - 1).abs() < .01]
            if match.width_days.nunique() != 2:
                continue
            result = row.to_dict()
            result.update(matching_detrend_widths=2,
                          minimum_width_local_snr=float(match.groupby("width_days").local_snr.max().min()),
                          same_data_correlated_widths=True,
                          duration_fraction=float(row.duration_days / row.period_days),
                          priority_scope="Operational queue only; correlated filter widths do not calibrate significance; known hosts with new periods remain eligible.")
            queue.append(result)
    ranked = pd.DataFrame(queue)
    if len(ranked):
        ranked = ranked.sort_values(["minimum_width_local_snr", "local_snr"], ascending=False).drop_duplicates(["tic", "sector"])
    ranked.to_csv(OUT / "ranked_width_robust_all.csv", index=False)
    unknown = ranked[~ranked.known_period_match & ~ranked.known_harmonic_match].drop_duplicates("tic") if len(ranked) else ranked
    unknown.head(6).to_csv(OUT / "priority_unknown_top6.csv", index=False)
    source_paths = [INPUT / "signals.csv", INPUT / "error_recovery/signals.csv", INPUT / "targets.csv", ROOT / "src/multisector_recognition.py", ROOT / "src/broad_crossmatch.py", ROOT / "scripts/pooled_catalog_review.py", Path(__file__)]
    summary = {"created_utc": m.utc(), "input_sources": {str(path.relative_to(ROOT)): m.sha(path.read_bytes()) for path in source_paths},
               "producer_rows": len(main), "empirical_rows": len(recovery), "total_signal_rows": len(frame),
               "unique_stars": int(frame.tic.nunique()), "unique_target_sector_observations": len(frame[["tic", "sector"]].drop_duplicates()),
               "eligible_width_robust_observations": len(ranked), "unmatched_width_robust_stars": len(unknown),
               "catalog_provenance": ix.metadata(), "positive_controls_passed": True,
               "scope": "Finite current frozen catalog index; absence is not novelty. Error streams tagged, no independent significance from widths/ranks. Known-host different periods retained. Only direct/integer-harmonic reported periods removed from unknown queue; rational aliases still need vetting."}
    m.write_json(OUT / "summary.json", summary)
    print(json.dumps({k: summary[k] for k in ["total_signal_rows", "unique_stars", "unique_target_sector_observations", "eligible_width_robust_observations", "unmatched_width_robust_stars"]}), flush=True)
    print(unknown.head(6)[["tic", "sector", "period_days", "local_depth_ppm", "minimum_width_local_snr", "uncertainty_stream", "raw_path"]].to_string(index=False), flush=True)


def live(tics):
    m.RESULTS = OUT / "live_context" / ("tics_" + "_".join(map(str, tics)))
    m.RESULTS.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as executor:
        paths = list(executor.map(lambda tic: str(m.capture_live_context([tic])), tics))
    m.write_json(m.RESULTS / "sources.json", {"live": paths, "tics": tics})
    deep = m.capture_deep_context(tics)
    m.write_json(m.RESULTS / "sources.json", {"live": paths, "deep": str(deep), "tics": tics})
    print(json.dumps({"live": paths, "deep": str(deep)}), flush=True)


def spoc(tics):
    sys.path.insert(0, str(ROOT / "scripts"))
    import pooled_catalog_review as helper
    helper.OUT = OUT / "spoc_context" / ("tics_" + "_".join(map(str, tics)))
    helper.OUT.mkdir(parents=True, exist_ok=True)
    helper.supplemental(tics, helper.OUT)
    helper.dvxml()
    helper.dvsummary()


def neighbor790():
    directory = m.BASE / "deep" / (m.utc().replace("-", "").replace(":", "").replace(".", "") + "_b02_neighbor790")
    directory.mkdir(parents=True, exist_ok=False)
    sid = "6476196643780163456"
    jobs = []
    for label, table in [("gaia_eb", "vari_eclipsing_binary"), ("gaia_nss", "nss_two_body_orbit")]:
        query = "select * from gaiadr3." + table + " where source_id=" + sid
        url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
        jobs.append((label, url, lambda b: m.parse_csv(b, required_columns=("source_id",), min_rows=0)[0]))
    url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source": "B/vsx/vsx", "-c": "314.9530591491028 -53.988804937631365", "-c.rs": "10", "-out": "OID,Name,RAJ2000,DEJ2000,Type,Period", "-out.max": "unlimited"})
    jobs.append(("vsx", url, m.parse_vsx_tsv))
    result = {"directory": str(directory.relative_to(ROOT)), "tic_aperture": 79004991, "neighbor_gaia_dr3": sid,
              "scope": "Neighbor source, not assigned to target; period agreement still requires pixel localization."}
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [(label, executor.submit(m.retain_public_response, directory, label + (".tsv" if label == "vsx" else ".csv"), url, validator)) for label, url, validator in jobs]
        for label, future in futures:
            result[label] = future.result()
            result[label + "_valid"] = result[label] is not None
    m.write_json(OUT / "neighbor790_context.json", result)
    print(json.dumps(result), flush=True)


def focus229():
    directory = m.BASE / "deep" / (m.utc().replace("-", "").replace(":", "").replace(".", "") + "_b02_focus229")
    directory.mkdir(parents=True, exist_ok=False)
    sid = "4913366449740781696"
    result = {"directory": str(directory.relative_to(ROOT)), "tic": 229146531, "gaia_dr3": sid}
    query = "select * from gaiadr3.vari_rotation_modulation where source_id=" + sid
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    rows = m.retain_public_response(directory, "gaia_rotation_modulation.csv", url, lambda b: m.parse_csv(b, required_columns=("source_id",), min_rows=0)[0])
    result["gaia_rotation_modulation"] = rows
    result["gaia_rotation_modulation_valid"] = rows is not None
    # WDS target cone is a recognition check, not complete close-companion imaging.
    url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source": "B/wds/wds", "-c": "27.118916 -52.138940", "-c.rs": "60", "-out.all": "", "-out.max": "unlimited"})
    rows = m.retain_public_response(directory, "wds_cone.tsv", url, lambda b: m.parse_asu_tsv(b, "B/wds/wds", ("WDS",)))
    result["wds"] = rows
    result["wds_valid"] = rows is not None
    result["scope"] = "Direct Gaia rotational-modulation period and WDS60arcsec recognition only. Existing target Gaia321day orbit does not identify the0.935day photometry. Neither variable classification nor a binary catalog entry alone excludes a planet."
    m.write_json(OUT / "focus229_context.json", result)
    print(json.dumps(result), flush=True)


def literature229():
    import xml.etree.ElementTree as ET
    from astropy.io import ascii
    import numpy as np
    sys.path.insert(0, str(ROOT / "scripts"))
    from pooled_catalog_review import separation
    directory = m.BASE / "deep" / (m.utc().replace("-", "").replace(":", "").replace(".", "") + "_b02_literature229")
    directory.mkdir(parents=True, exist_ok=False)
    result = {"directory": str(directory.relative_to(ROOT)), "tic": 229146531, "retained_published_sources": {}, "retained_published_matches": {}}
    # Two complete published EB lists with their original publisher manifests.
    for number, count in [(3, 7936), (4, 2065)]:
        path = ROOT / f"data/stellar_systems/scout_novelty/tess10000_table{number}_mrt.txt"
        body = path.read_bytes()
        manifest = json.loads(path.with_suffix(".txt.manifest.json").read_text())
        assert m.sha(body) == manifest["sha256"]
        table = ascii.read(path, format="cds")
        assert len(table) == count
        if number == 3:
            assert 197764089 in table["TIC"]
        matches = []
        for row in table:
            sep = separation(27.118943046382501, -52.138937892653601, {"ra": row["RAdeg"], "dec": row["DEdeg"]})
            if int(row["TIC"]) == 229146531 or sep < 120:
                record = {key: None if np.ma.is_masked(row[key]) else row[key].item() if isinstance(row[key], np.generic) else row[key] for key in table.colnames}
                matches.append({"row": record, "separation_arcsec": sep})
        result["retained_published_sources"][str(path.relative_to(ROOT))] = {"sha256": manifest["sha256"], "rows": count, "source_url": manifest["url"]}
        result["retained_published_matches"][f"kostov2025_table{number}"] = matches
    jobs = []
    for label, table, query, required in [
        ("kostov_unvetted_table2", 'J/ApJS/279/50/table2', 'select * from "J/ApJS/279/50/table2" where TIC in (229146531,197764089)', ("TIC",)),
        ("shi2022_ea", 'J/ApJS/259/50/table1', 'select * from "J/ApJS/259/50/table1" where TIC=229146531', ("TIC",)),
        ("tang2025_eccentric", 'J/AJ/170/356/table1', 'select * from "J/AJ/170/356/table1" where TIC=229146531', ("TIC",))]:
        url = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
        jobs.append((label, ".csv", url, lambda b, cols=required: m.parse_csv(b, required_columns=cols, min_rows=0)[0]))
    query = "select dr3_source_id,dr2_source_id,angular_distance,magnitude_difference from gaiadr3.dr2_neighbourhood where dr3_source_id=4913366449740781696"
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    jobs.append(("dr2_neighborhood", ".csv", url, lambda b: m.parse_csv(b, required_columns=("dr3_source_id", "dr2_source_id"), min_rows=0)[0]))
    query = "select * from gaiadr3.astrophysical_parameters where source_id=4913366449740781696"
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    jobs.append(("gaia_astrophysical", ".csv", url, lambda b: m.parse_csv(b, required_columns=("source_id",), min_rows=0)[0]))
    # Actual table rows, not a null inferred from the surrounding paper.
    url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source": "J/ApJ/930/7/table1", "ID": "229146531", "-out.all": "", "-out.max": "unlimited"})
    jobs.append(("avallone2022_exact_tic", ".tsv", url, lambda b: m.parse_asu_tsv(b, "J/ApJ/930/7/table1", ("ID", "Class"))))
    def parse_atom(body):
        root = ET.fromstring(body)
        assert root.tag == "{http://www.w3.org/2005/Atom}feed"
        ns = {"a": "http://www.w3.org/2005/Atom"}
        rows = [{"id": e.findtext("a:id", namespaces=ns), "title": e.findtext("a:title", namespaces=ns)} for e in root.findall("a:entry", ns)]
        assert not any("api/errors" in (row["id"] or "") for row in rows)
        return rows
    aliases_query = ' OR '.join('all:"' + alias + '"' for alias in ["TIC 229146531", "TYC 8047-270-1", "UCAC4 190-001692", "2MASS J01482854-5208201", "4913366449740781696"])
    url = "https://export.arxiv.org/api/query?" + urlencode(dict(search_query=aliases_query, max_results=100))
    jobs.append(("alias_arxiv", ".xml", url, parse_atom))
    def html(body):
        text = body.decode("utf-8", errors="replace")
        if "<html" not in text.lower() or "cloudflare" in text.lower():
            raise ValueError("Unusable public HTML response")
        return {"bytes": len(body), "content_scope": "Retained HTML only. Membership/notes require inspecting actual table content; retrieval is not a null."}
    jobs.append(("tesseb_live_html", ".html", "https://tessebs.villanova.edu/0229146531", html))
    jobs.append(("exofop_target_html", ".html", "https://exofop.ipac.caltech.edu/tess/target.php?id=229146531", html))
    jobs.append(("vsx_detail_html", ".html", "https://www.aavso.org/vsx/index.php?view=detail.top&oid=2748095", html))
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [(label, executor.submit(m.retain_public_response, directory, label + suffix, url, validator)) for label, suffix, url, validator in jobs]
        for label, future in futures:
            parsed = future.result()
            result[label] = parsed
            result[label + "_valid"] = parsed is not None
    result["scope"] = "Finite alias and actual published-catalog checks, not all literature, all TCEs or a discovery assessment. Unvetted Kostov inventory membership, if any, is a prior flag, not planet/EB validation. HTML retrieval alone is not a catalog-null test."
    m.write_json(OUT / "literature229_context.json", result)
    print(json.dumps(result), flush=True)


def catalog229_controls():
    from bs4 import BeautifulSoup
    directory = m.BASE / "deep" / (m.utc().replace("-", "").replace(":", "").replace(".", "") + "_b02_catalog229_controls")
    directory.mkdir(parents=True, exist_ok=False)
    result = {"directory": str(directory.relative_to(ROOT)), "tic": 229146531}
    query = 'select TOP 1 TIC from "J/ApJS/279/50/table2"'
    base = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync?"
    url = base + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    sample = m.retain_public_response(directory, "kostov_table2_positive_sample.csv", url, lambda b: m.parse_csv(b, required_columns=("TIC",), min_rows=1)[0])
    result["kostov_table2_positive_sample"] = sample
    if sample:
        source_id = sample[0]["TIC"]
        query = 'select * from "J/ApJS/279/50/table2" where TIC=' + source_id
        url = base + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
        rows = m.retain_public_response(directory, "kostov_table2_positive_recovery.csv", url, lambda b: m.parse_csv(b, required_columns=("TIC",), min_rows=1)[0])
        result["kostov_table2_positive_recovery"] = rows
        result["kostov_positive_control_passed"] = rows is not None and any(row["TIC"] == source_id for row in rows)
    literature = json.loads((OUT / "literature229_context.json").read_text())
    source_directory = ROOT / literature["directory"]
    pages = {}
    for name in ["tesseb_live_html.html", "vsx_detail_html.html", "exofop_target_html.html"]:
        path = source_directory / name
        body = path.read_bytes()
        manifest = json.loads(path.with_suffix(".html.manifest.json").read_text())
        assert manifest["valid"] and m.sha(body) == manifest["sha256"]
        soup = BeautifulSoup(body, "html.parser")
        for tag in soup(["script", "style"]):
            tag.decompose()
        text = soup.get_text(" ", strip=True)
        (directory / (name + ".text.txt")).write_text(text + "\n")
        pages[name] = {"title": soup.title.get_text() if soup.title else None, "source_sha256": manifest["sha256"], "source": str(path.relative_to(ROOT))}
        if name == "tesseb_live_html.html":
            assert "TIC 229146531" in text and "This TESS ID is not in the catalog." in text
            pages[name]["reviewed_statement"] = "Actual target page explicitly says this TESS ID is not in this live catalog. Finite membership scope only."
        elif name == "vsx_detail_html.html":
            assert "Gaia DR3 4913366449740781696" in text and "Period --" in text
            pages[name]["reviewed_statement"] = "ROT; period blank; G11.39-11.44; reference GaiaDR3variability2022; no public remarks/supporting documents. No independently measured rotation period supplied."
        else:
            assert "TIC 229146531" in text and "Observing Notes 0" in text
            pages[name]["reviewed_statement"] = "Actual public target page reports ObservingNotes0,Files0,TOIs0,CTOIs0. This does not inventory every team's TCEs or unpublished work."
    result["reviewed_public_pages"] = pages
    result["scope"] = "Catalog positive control plus actual target page content, not exhaustive priority proof; private/unpublished notes are outside scope."
    m.write_json(OUT / "catalog229_controls.json", result)
    print(json.dumps({"positive_control_passed": result.get("kostov_positive_control_passed"), "reviewed_pages": len(pages)}), flush=True)


def context():
    sys.path.insert(0, str(ROOT / "scripts"))
    from pooled_catalog_review import separation
    sources, records = [], {}
    for path in sorted((OUT / "live_context").rglob("sources.json")):
        source = json.loads(path.read_text())
        if "deep" not in source:
            raise ValueError("Incomplete lead context: " + str(path))
        sources.append(source)
        live = {}
        for directory in source["live"]:
            for item in json.loads((Path(directory) / "summary.json").read_text())["checks"]:
                live[int(item["tic"])] = item
        for deep in json.loads((Path(source["deep"]) / "summary.json").read_text())["checks"]:
            tic = int(deep["tic"])
            live_record = live[tic]
            target_path = ROOT / deep["exofop_source"]
            body = target_path.read_bytes()
            assert m.sha(body) == deep["exofop_sha256"]
            target = json.loads(body)
            target_gaia = next(g for g in deep["gaia_neighbors"] if g["source_id"] == deep["gaia_dr3"])
            record = {"tic": tic, "live": live_record, "deep": deep, "target_parameters": target,
                      "target_gaia": target_gaia}
            for collection in [live_record.get("vsx_rows") or [], deep.get("asassn") or [], deep.get("simbad_bibliography") or [], deep.get("gaia_neighbors") or []]:
                for row in collection:
                    row["separation_arcsec"] = separation(deep["ra"], deep["dec"], row)
                    row["association_scope"] = "same_position_within2arcsec_at_catalog_coordinates" if row["separation_arcsec"] is not None and row["separation_arcsec"] < 2 else "neighbor_not_assigned_to_target"
            records[tic] = record
    notes = {
        271640189: "P2.459959d matches nearby known EA ASASSN J071536.11-721806.9 / TIC271640184 P2.4603d, ~38arcsec; localization required.",
        79004991: "P1.186706d near twice known neighbor TIC219298873 period0.59322d. Independent ASAS-SN/VSX/GaiaEB/GaiaSB1 agree; ~59.7arcsec; localization required.",
        139526205: "Same-position known KELT KS30C13480 EA P2.922027d; actual 2018 KELT transit false-positive-catalog bibliography recovered.",
        229112973: "Same-position known ASASJ011031-4929.4 ROT P1.051d; GaiaRSclassifier and significant RV variability. SPOCTCE12.017d also prior public context, different period.",
        229146531: "Same-position VSXROT/GaiaSOLAR_LIKE, published GaiaOrbitalP321.3894d; no assignment of P0.934911d to that long-period orbit.",
        71274048: "PublicSPOCTCE158.828374d is a different period. Known bright comoving pair ~50arcsec; contamination.4605. P1.359244d is not identified by catalog context.",
        147271811: "Same-position VSXROT/GaiaSOLAR_LIKE. Two prior public SPOCTCEs near0.28746d with opposite epochs; new0.57445d peak is near their2xharmonic.",
        389752955: "Target not known variable in these finite queries; nearby23.4arcsec G14.06 GaiaSOLAR_LIKE/VSXROT neighbor6526819464795267328 requires localization.",
        80204331: "No matched prior period in checked sources; GaiaIPDmulti_peak22percent indicates source/blend warning. No novelty conclusion.",
        80427865: "Known CD-45216A accelerationbinary/RUWE7.13; nearby CD-45216B/TIC80427866 is knownROT1.4742d. That period is not identified with follow-up2.165521d peak."
    }
    rows = []
    for tic, record in records.items():
        g, target = record["target_gaia"], record["target_parameters"]
        rows.append({"tic": tic, "aliases": record["deep"]["aliases"], "gaia_dr3": record["deep"]["gaia_dr3"],
                     "ra": record["deep"]["ra"], "dec": record["deep"]["dec"],
                     "gaia_ruwe": g["ruwe"], "gaia_nss_flag": g["non_single_star"], "gaia_best_class": g["best_class_name"],
                     "gaia_rv_amplitude_robust_kms": g["rv_amplitude_robust"], "gaia_rv_chisq_pvalue": g["rv_chisq_pvalue"],
                     "tic_contamination_ratio": target["basic_info"]["tic_contamination_ratio"],
                     "assessment": notes[tic]})
    pd.DataFrame(rows).to_csv(OUT / "lead_catalog_context.csv", index=False)
    m.write_json(OUT / "lead_catalog_context.json", {"created_utc": m.utc(), "sources": sources, "checks": list(records.values()),
        "method_sha256": m.sha(Path(__file__).read_bytes()),
        "limitations": "Selectively checked prioritized targets only; prior-variable classifications do not themselves disprove planets or localize a TESS signal. Nearby-source references and ephemerides must not be assigned to the selected target. No null proves novelty; exact-TIC arXiv searches are incomplete. All HTTP responses and failures retained."})
    print(json.dumps({"lead_context_stars": len(records)}), flush=True)


def verify():
    summary = json.loads((OUT / "summary.json").read_text())
    for name, expected in summary["input_sources"].items():
        assert m.sha((ROOT / name).read_bytes()) == expected, name
    context = json.loads((OUT / "lead_catalog_context.json").read_text())
    assert context["method_sha256"] == m.sha(Path(__file__).read_bytes())
    directories = []
    for source in context["sources"]:
        directories.extend(Path(x) for x in source["live"])
        directories.append(Path(source["deep"]))
    for path in list((OUT / "spoc_context").rglob("spoc_dv_inventory.json")) + list((OUT / "spoc_context").rglob("spoc_dv_xml.json")) + [OUT / "neighbor790_context.json"] + list(OUT.glob("focus229_context.json")) + list(OUT.glob("literature229_context.json")) + list(OUT.glob("catalog229_controls.json")):
        directories.append(ROOT / json.loads(path.read_text())["directory"])
    checked = {}
    for directory in directories:
        for manifest in directory.glob("*.manifest.json"):
            meta = json.loads(manifest.read_text())
            if "sha256" not in meta:
                continue
            body = manifest.with_name(manifest.name.removesuffix(".manifest.json"))
            assert m.sha(body.read_bytes()) == meta["sha256"], body
            checked[str(body.relative_to(ROOT))] = meta["sha256"]
    evidence = {"verified_utc": m.utc(), "input_sources_valid": True, "retained_http_bodies_checked": len(checked), "response_hashes": checked,
                "scope": "Integrity checks, not scientific confirmation. Raw observation hashes are separately retained in annotated rows; public query bodies include failures."}
    m.write_json(OUT / "verification.json", evidence)
    print(json.dumps({"passed": True, "http_bodies": len(checked)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["annotate", "live", "spoc", "neighbor790", "focus229", "literature229", "catalog229-controls", "context", "verify"])
    parser.add_argument("--tic", type=int, action="append")
    args = parser.parse_args()
    if args.command == "annotate":
        annotate()
    elif args.command == "neighbor790":
        neighbor790()
    elif args.command == "focus229":
        focus229()
    elif args.command == "literature229":
        literature229()
    elif args.command == "catalog229-controls":
        catalog229_controls()
    elif args.command == "context":
        context()
    elif args.command == "verify":
        verify()
    elif args.command == "live" and args.tic:
        live(args.tic)
    elif args.command == "spoc" and args.tic:
        spoc(args.tic)
    else:
        parser.error("live needs --tic")


if __name__ == "__main__":
    main()
