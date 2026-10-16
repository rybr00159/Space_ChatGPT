"""Immutable fresh catalog snapshots and recognition for the multi-sector search.

Catalog absence is not novelty. The index retains all dispositions, RV planets,
published candidates, and nearby contaminants. Period/harmonic matches are
screening flags and require ephemeris and pixel checks before identification.
The earlier broad_crossmatch.py and all its snapshots remain unchanged.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

import requests

from broad_crossmatch import CatalogIndex, numeric

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/multisector/catalogs"
RESULTS = ROOT / "results/multisector_recognition"
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
CONFIRMED_COLUMNS = (
    "tran_flag,rv_flag,pl_name,hostname,tic_id,gaia_dr2_id,ra,dec,sy_pnum,"
    "sy_snum,discoverymethod,pl_orbper,pl_orbpererr1,pl_orbpererr2,"
    "pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandur,pl_controv_flag,"
    "pl_orbper_reflink,pl_tranmid_reflink,pl_trandur_reflink"
).split(",")
TOI_COLUMNS = (
    "tid,toi,toidisplay,tfopwg_disp,ra,dec,pl_tranmid,pl_tranmiderr1,"
    "pl_tranmiderr2,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_trandurh,"
    "sectors,rowupdate,release_date"
).split(",")


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(body):
    return hashlib.sha256(body).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def parse_csv(body, expected_columns=None, required_columns=(), min_rows=1):
    text = body.decode("utf-8-sig")
    if text.lstrip().startswith(("<", "upstream ", "Internal Server")):
        raise ValueError("HTTP response is not a CSV table")
    reader = csv.DictReader(io.StringIO(text))
    columns = reader.fieldnames
    if not columns or any(x not in columns for x in required_columns):
        raise ValueError(f"Missing required CSV columns: {required_columns}; got {columns}")
    if expected_columns is not None and columns != list(expected_columns):
        raise ValueError(f"Changed table schema; expected {expected_columns}, got {columns}")
    rows = list(reader)
    if len(rows) < min_rows or any(None in r for r in rows):
        raise ValueError(f"CSV is truncated or too short: {len(rows)} rows")
    return rows, columns


def fetch_csv(directory, name, url, expected_columns, min_rows, query=None):
    """Retain even failed HTTP responses; never accept them as catalog absence."""
    path = directory / name
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite catalog snapshot {path}")
    response = requests.get(url, timeout=(20, 180),
                            headers={"User-Agent": "Space_ChatGPT-public-research/1.0"})
    body = response.content
    path.write_bytes(body)
    metadata = dict(retrieved_utc=utc(), url=url, query=query,
                    http_status=response.status_code, bytes=len(body),
                    sha256=sha(body), content_type=response.headers.get("Content-Type"),
                    valid=False)
    try:
        response.raise_for_status()
        rows, columns = parse_csv(body, expected_columns, min_rows=min_rows)
        metadata.update(valid=True, row_count=len(rows), columns=columns)
    except Exception as exc:
        metadata["validation_error"] = str(exc)
        raise
    finally:
        write_json(path.with_suffix(path.suffix + ".manifest.json"), metadata)
    return metadata


def retrieve_snapshot():
    """Retrieve new core tables; an incomplete snapshot never becomes latest."""
    directory = BASE / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    directory.mkdir(parents=True, exist_ok=False)
    # ExoFOP's entire public schema is checked against the prior successful
    # snapshot. A changed field requires a deliberate parser update.
    ctoi_columns = json.loads((ROOT / "data/catalogs/ctoi.csv.manifest.json").read_text())["columns"]
    jobs = []
    for name, columns, table, minimum in [
            ("confirmed_all.csv", CONFIRMED_COLUMNS, "pscomppars", 5000),
            ("toi.csv", TOI_COLUMNS, "toi", 6000)]:
        query = "select " + ",".join(columns) + " from " + table
        jobs.append((name, TAP + "?" + urlencode(dict(query=query, format="csv")),
                     columns, minimum, query))
    jobs.append(("ctoi.csv", "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?output=csv",
                 ctoi_columns, 1000, None))
    failures, catalogs = {}, {}
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [(job[0], pool.submit(fetch_csv, directory, *job)) for job in jobs]
        for name, future in futures:
            try:
                catalogs[name] = future.result()
            except Exception as exc:
                failures[name] = str(exc)
    summary = dict(created_utc=utc(), snapshot=str(directory.relative_to(ROOT)),
                   catalogs=catalogs, failures=failures,
                   complete=not failures,
                   limitations="Fresh public snapshots are finite prior-report checks, not a complete literature or all-TCE inventory.")
    write_json(directory / "snapshot.json", summary)
    if failures:
        raise RuntimeError(f"Incomplete snapshot retained in {directory}: {failures}")
    write_json(BASE / "latest.json", summary)
    return directory


def verified_csv(path):
    body = path.read_bytes()
    manifest = json.loads(path.with_suffix(path.suffix + ".manifest.json").read_text())
    if manifest.get("sha256") != sha(body) or manifest.get("valid") is False:
        raise ValueError(f"Invalid or changed catalog snapshot: {path}")
    rows, columns = parse_csv(body, manifest.get("columns"))
    return rows, {**manifest, "file": str(path), "columns": columns}


def parse_asu_tsv(body, table, required):
    """Validate ASU service/table metadata before allowing a zero-result cone."""
    text = body.decode("utf-8-sig")
    if "#INFO\tError=" in text or "Table or Catalog not found" in text or "<html" in text.lower():
        raise ValueError("VizieR returned a service/catalog error")
    if "#Name: " + table not in text:
        raise ValueError("Missing expected ASU table metadata")
    if not all("#Column\t" + name + "\t" in text for name in required):
        raise ValueError("Changed VSX metadata schema")
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return []  # Valid ASU no-row response carries its table/column metadata.
    columns = lines[0].split("\t")
    if not all(name in columns for name in required) or len(lines) < 3:
        raise ValueError("Malformed VSX TSV header")
    rows = list(csv.DictReader([lines[0]] + lines[3:], delimiter="\t"))
    return [{k: v.strip() for k, v in row.items()} for row in rows]


def parse_vsx_tsv(body):
    return parse_asu_tsv(body, "B/vsx/vsx", ("OID", "Name", "RAJ2000", "DEJ2000", "Type", "Period"))


def retain_public_response(directory, name, url, validator):
    """Lead-only fetch preserving response bytes, errors and validation status."""
    path = directory / name
    if path.exists():
        raise FileExistsError(f"Refusing to replace immutable live response {path}")
    for old_manifest in sorted(BASE.rglob(name + ".manifest.json"), reverse=True):
        prior = json.loads(old_manifest.read_text())
        if prior.get("url") != url or not prior.get("valid"):
            continue
        old_path = old_manifest.with_name(old_manifest.name.removesuffix(".manifest.json"))
        old_body = old_path.read_bytes()
        if sha(old_body) != prior["sha256"]:
            raise ValueError(f"Changed retained public response {old_path}")
        parsed = validator(old_body)
        path.write_bytes(old_body)
        write_json(path.with_suffix(path.suffix + ".manifest.json"),
                   {**prior, "reused_utc": utc(), "reused_from": str(old_path.relative_to(ROOT))})
        return parsed
    metadata = dict(url=url, retrieved_utc=utc(), valid=False)
    try:
        response = requests.get(url, timeout=(15, 75),
                                headers={"User-Agent": "Space_ChatGPT-public-research/1.0"})
        body = response.content
        path.write_bytes(body)
        metadata.update(http_status=response.status_code, bytes=len(body), sha256=sha(body))
        response.raise_for_status()
        parsed = validator(body)
        metadata["valid"] = True
        return parsed
    except Exception as exc:
        metadata["error"] = str(exc)
        return None
    finally:
        write_json(path.with_suffix(path.suffix + ".manifest.json"), metadata)


def raven_sources(root=ROOT):
    """Use the versioned RAVEN release, verifying both SHA and record MD5.

    This release was already downloaded for the TIC 71316387 audit. Keeping its
    version/date explicit avoids presenting a reused release as a live query.
    Failed 91-byte download attempts are deliberately not eligible inputs.
    """
    parent = root / "candidates/tic71316387/catalog"
    record_path = parent / "zenodo.json"
    record_body = record_path.read_bytes()
    record_manifest = json.loads((parent / "zenodo.json.manifest.json").read_text())
    if sha(record_body) != record_manifest["sha256"]:
        raise ValueError("Changed RAVEN release record")
    record = json.loads(record_body)
    sources = []
    for item in record["files"]:
        if not item["key"].endswith("_table.csv"):
            continue
        expected_md5 = item["checksum"].removeprefix("md5:")
        valid = []
        for path in sorted(parent.glob(Path(item["key"]).stem + "*.csv")):
            body = path.read_bytes()
            if hashlib.md5(body).hexdigest() == expected_md5:
                valid.append((path, body))
        if not valid:
            raise ValueError(f"Missing checksum-valid RAVEN release file {item['key']}")
        path, body = valid[0]
        rows, columns = parse_csv(body, required_columns=("ticid", "bls_per", "bls_t0", "bls_tdur"))
        sources.append((item["key"], rows, dict(
            file=str(path), sha256=sha(body), bytes=len(body), row_count=len(rows),
            url=item["links"]["self"], release_record_url=record_manifest["url"],
            release_record_sha256=sha(record_body),
            retrieved_utc=record_manifest["retrieved_utc"], md5=expected_md5,
            columns=columns, scope="RAVEN versioned public release 19661443; reused immutable 2026-10-14 snapshot, not a new live query.")))
    return sources


def capture_live_context(tics):
    """Small lead-only public target/alias checks; retain failures explicitly.

    ExoFOP's target JSON exposes notes and existing files that a bulk period
    table may omit. SIMBAD is a finite cone, arXiv an exact TIC-text search;
    neither is a comprehensive literature/alias or novelty determination.
    """
    directory = BASE / "live" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    directory.mkdir(parents=True, exist_ok=False)
    checks = []

    def fetch(name, url, validator):
        path = directory / name
        metadata = dict(url=url, retrieved_utc=utc(), valid=False)
        try:
            response = requests.get(url, timeout=(15, 75),
                                    headers={"User-Agent": "Space_ChatGPT-public-research/1.0"})
            body = response.content
            path.write_bytes(body)
            metadata.update(http_status=response.status_code, bytes=len(body), sha256=sha(body))
            response.raise_for_status()
            parsed = validator(body)
            metadata["valid"] = True
            return parsed
        except Exception as exc:
            metadata["error"] = str(exc)
            return None
        finally:
            write_json(path.with_suffix(path.suffix + ".manifest.json"), metadata)

    for tic in dict.fromkeys(int(t) for t in tics):
        def exofop_validator(body):
            parsed = json.loads(body)
            if int(parsed["basic_info"]["tic_id"]) != tic:
                raise ValueError("ExoFOP returned a different target identifier")
            if not all(isinstance(parsed.get(k), list) for k in ("tois", "ctois", "files")):
                raise ValueError("Changed ExoFOP target JSON schema")
            return parsed

        target = fetch(f"tic{tic}_exofop.json",
                       f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json",
                       exofop_validator)
        record = dict(tic=tic, checked_utc=utc(), exofop_valid=target is not None,
                      scope="Finite public target JSON, 30-arcsec SIMBAD cone and exact TIC-text arXiv check; nulls do not establish novelty.")
        if target is not None:
            record.update(aliases=target["basic_info"].get("star_names"),
                          tois=target["tois"], ctois=target["ctois"], files=target["files"],
                          planet_parameters=target.get("planet_parameters", []))
            ra, dec = numeric(target["coordinates"].get("ra")), numeric(target["coordinates"].get("dec"))
            if ra is not None and dec is not None:
                query = ("select main_id,ra,dec,otype,nbref from basic where "
                         f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{30/3600}))")
                url = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync?" + urlencode(
                    dict(request="doQuery", lang="ADQL", format="csv", query=query))
                rows = fetch(f"tic{tic}_simbad.csv", url,
                             lambda body: parse_csv(body, required_columns=("main_id", "ra", "dec", "otype", "nbref"), min_rows=0)[0])
                record.update(simbad_valid=rows is not None, simbad_rows=rows)
                vsx_url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({
                    "-source": "B/vsx/vsx", "-c": f"{ra} {dec}", "-c.rs": "30",
                    "-out": "OID,Name,RAJ2000,DEJ2000,Type,Period", "-out.max": "unlimited"})
                variables = fetch(f"tic{tic}_vsx.tsv", vsx_url, parse_vsx_tsv)
                record.update(vsx_valid=variables is not None, vsx_rows=variables)

        def arxiv_validator(body):
            feed = ET.fromstring(body)
            if feed.tag != "{http://www.w3.org/2005/Atom}feed":
                raise ValueError("Response is not an arXiv Atom feed")
            entries = feed.findall("{http://www.w3.org/2005/Atom}entry")
            parsed = []
            for entry in entries:
                get = lambda name: entry.findtext("{http://www.w3.org/2005/Atom}" + name, "")
                if "api/errors" in get("id"):
                    raise ValueError("arXiv returned an API error entry")
                parsed.append(dict(id=get("id"), title=get("title"), summary=get("summary")))
            return parsed

        arxiv_url = "https://export.arxiv.org/api/query?" + urlencode(
            dict(search_query=f'all:"TIC {tic}" OR all:"{tic}"', max_results=100))
        papers = fetch(f"tic{tic}_arxiv.xml", arxiv_url, arxiv_validator)
        record.update(arxiv_valid=papers is not None, arxiv_entries=papers)
        checks.append(record)
    write_json(directory / "summary.json", dict(checks=checks, created_utc=utc()))
    return directory


def summarize_live_context():
    """Derive angular separations without changing immutable HTTP responses."""
    records, sources = [], {}
    for path in sorted((BASE / "live").glob("*/summary.json")):
        body = path.read_bytes()
        sources[str(path.relative_to(ROOT))] = sha(body)
        for original in json.loads(body)["checks"]:
            record = dict(original)
            target_path = path.parent / f"tic{record['tic']}_exofop.json"
            if record["exofop_valid"]:
                target_body = target_path.read_bytes()
                manifest = json.loads(target_path.with_suffix(".json.manifest.json").read_text())
                if sha(target_body) != manifest["sha256"]:
                    raise ValueError("Changed live ExoFOP input")
                coordinates = json.loads(target_body)["coordinates"]
                ra, dec = float(coordinates["ra"]), float(coordinates["dec"])
                record.update(ra=ra, dec=dec)
                for row in record.get("vsx_rows", []) or []:
                    r1, d1, r2, d2 = map(math.radians, [ra, dec, float(row["RAJ2000"]), float(row["DEJ2000"])])
                    hav = math.sin((d2-d1)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((r2-r1)/2)**2
                    sep = math.degrees(2*math.asin(min(1, math.sqrt(hav))))*3600
                    row["separation_arcsec"] = sep
                    row["association_scope"] = "same_position_within_2arcsec" if sep <= 2 else "nearby_possible_contaminant_not_assigned_to_target"
            records.append(record)
    result = dict(created_utc=utc(), sources=sources, checks=records,
                  limitations="Only exact-TIC arXiv text, finite VSX/SIMBAD cones, and ExoFOP target page checked. No exhaustive alias/literature search, proper-motion propagation, or source localization is implied.")
    write_json(RESULTS / "live_context_summary.json", result)
    return result


def capture_deep_context(tics):
    """Lead-only stellar/neighbor properties and actual alias bibliography.

    Gaia's variability classification is prior context, not proof of the source
    of a TESS event. The ASAS-SN table is its published 2021 catalog, and the
    SIMBAD bibliography can omit papers and targets/aliases. Retain this scope.
    """
    directory = BASE / "deep" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    directory.mkdir(parents=True, exist_ok=False)
    records, jobs = {}, []

    def tap_url(service, query):
        base = {"gaia": "https://gea.esac.esa.int/tap-server/tap/sync",
                "simbad": "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"}[service]
        return base + "?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))

    for tic in dict.fromkeys(int(t) for t in tics):
        paths = sorted((BASE / "live").glob(f"*/tic{tic}_exofop.json"), reverse=True)
        if not paths:
            raise FileNotFoundError(f"Run --live-tic {tic} before its deeper checks")
        target_path = paths[0]
        body = target_path.read_bytes()
        expected = json.loads(target_path.with_suffix(".json.manifest.json").read_text())
        if not expected["valid"] or sha(body) != expected["sha256"]:
            raise ValueError("Invalid or changed ExoFOP source for deeper checks")
        target = json.loads(body)
        if int(target["basic_info"]["tic_id"]) != tic:
            raise ValueError("Wrong deep-check target identifier")
        ra, dec = float(target["coordinates"]["ra"]), float(target["coordinates"]["dec"])
        aliases = target["basic_info"].get("star_names", "")
        source_match = re.search(r"Gaia DR3 (\d+)", aliases)
        source_id = source_match.group(1) if source_match else None
        records[tic] = dict(tic=tic, ra=ra, dec=dec, aliases=aliases, gaia_dr3=source_id,
                            exofop_source=str(target_path.relative_to(ROOT)), exofop_sha256=sha(body))
        cone = f"1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{ra},{dec},{60/3600}))"
        query = ("select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.phot_bp_mean_mag,g.phot_rp_mean_mag,"
                 "g.bp_rp,g.parallax,g.parallax_error,g.pmra,g.pmdec,g.ruwe,g.non_single_star,"
                 "g.phot_variable_flag,g.radial_velocity,g.rv_amplitude_robust,g.rv_chisq_pvalue,"
                 "g.teff_gspphot,g.ipd_frac_multi_peak,g.astrometric_excess_noise,"
                 "v.best_class_name,v.best_class_score from gaiadr3.gaia_source g "
                 "left outer join gaiadr3.vari_classifier_result v on g.source_id=v.source_id "
                 f"where {cone} order by g.phot_g_mean_mag")
        jobs.append((tic, "gaia_neighbors", tap_url("gaia", query),
                     lambda b: parse_csv(b, required_columns=("source_id", "ra", "dec", "best_class_name"), min_rows=0)[0]))
        if source_id:
            query = "select * from gaiadr3.vari_eclipsing_binary where source_id=" + source_id
            jobs.append((tic, "gaia_eb", tap_url("gaia", query),
                         lambda b: parse_csv(b, required_columns=("source_id",), min_rows=0)[0]))
            query = "select * from gaiadr3.vari_short_timescale where source_id=" + source_id
            jobs.append((tic, "gaia_short_timescale", tap_url("gaia", query),
                         lambda b: parse_csv(b, required_columns=("source_id",), min_rows=0)[0]))
            for key, table in [("gaia_nss_orbit", "nss_two_body_orbit"),
                               ("gaia_nss_acceleration", "nss_acceleration_astro")]:
                query = "select * from gaiadr3." + table + " where source_id=" + source_id
                jobs.append((tic, key, tap_url("gaia", query),
                             lambda b: parse_csv(b, required_columns=("source_id",), min_rows=0)[0]))
        columns = ("ASASSN-V", "RAJ2000", "DEJ2000", "Vmag", "Amp", "Per", "Type", "Prob", "TIC", "GaiaDR3")
        url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({
            "-source": "II/366/catv2021", "-c": f"{ra} {dec}", "-c.rs": "60",
            "-out": ",".join(columns), "-out.max": "unlimited"})
        jobs.append((tic, "asassn", url, lambda b, cols=columns: parse_asu_tsv(b, "II/366/catv2021", cols)))
        sim_cone = f"1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',{ra},{dec},{60/3600}))"
        query = ("select b.main_id,b.ra,b.dec,r.bibcode,r.doi,r.title,r.abstract from basic b "
                 "join has_ref h on h.oidref=b.oid join ref r on h.oidbibref=r.oidbib where " + sim_cone)
        jobs.append((tic, "simbad_bibliography", tap_url("simbad", query),
                     lambda b: parse_csv(b, required_columns=("main_id", "bibcode", "title"), min_rows=0)[0]))
        query = "select b.main_id,b.ra,b.dec,i.id from basic b join ident i on i.oidref=b.oid where " + sim_cone
        jobs.append((tic, "simbad_aliases", tap_url("simbad", query),
                     lambda b: parse_csv(b, required_columns=("main_id", "id"), min_rows=0)[0]))
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [(tic, key, pool.submit(retain_public_response, directory,
                                         f"tic{tic}_{key}." + ("tsv" if key == "asassn" else "csv"), url, validator))
                   for tic, key, url, validator in jobs]
        for tic, key, future in futures:
            value = future.result()
            records[tic][key] = value
            records[tic][key + "_valid"] = value is not None
    result = dict(created_utc=utc(), checks=list(records.values()),
                  limitations="Gaia DR3 classifiers/properties and 60arcsec neighbors, actual target NSS tables, published ASAS-SN 2021 cone, SIMBAD 60arcsec bibliography/aliases only. References retain source names/coordinates: nearby sources must not be assigned to the target. No catalog nonmatch is novelty; nearby variables are not localized TESS hosts.")
    write_json(directory / "summary.json", result)
    write_json(RESULTS / "deep_context_summary.json", {**result, "directory": str(directory.relative_to(ROOT))})
    return directory
class FreshCatalogIndex(CatalogIndex):
    """CatalogIndex-compatible fresh recognition; previous code is unchanged."""

    def __init__(self, root=ROOT, snapshot=None, **kwargs):
        root = Path(root)
        if snapshot is None:
            pointer = json.loads((root / "data/multisector/catalogs/latest.json").read_text())
            snapshot = root / pointer["snapshot"]
        self.snapshot = Path(snapshot)
        if not self.snapshot.is_absolute():
            self.snapshot = root / self.snapshot
        meta = json.loads((self.snapshot / "snapshot.json").read_text())
        if not meta["complete"]:
            raise ValueError("An incomplete catalog snapshot cannot be used")
        super().__init__(root, **kwargs)

    def _source(self, name, path, required=True):
        if name in ("confirmed", "toi", "ctoi"):
            rows, metadata = verified_csv(self.snapshot / path.name)
            metadata["snapshot_kind"] = "fresh_core"
            self.provenance[name] = metadata
            return rows
        rows = super()._source(name, path, required)
        if name in self.provenance:
            self.provenance[name]["snapshot_kind"] = "reused_published_support_2026-10-13"
        return rows

    def _load(self):
        super()._load()
        # Preserve the community report's author/group/notes in match evidence;
        # dispositions alone do not explain a prior RAVEN or other report.
        ctoi, _ = verified_csv(self.snapshot / "ctoi.csv")
        details = {"CTOI " + r["CTOI"]: r for r in ctoi}
        for row in self.rows:
            if row["catalog"] == "ctoi":
                raw = details[row["name"]]
                row.update(submitter=raw.get("User", ""), group=raw.get("Group", ""),
                           tag=raw.get("Tag", ""), notes=raw.get("Notes", ""),
                           last_modified=raw.get("CTOI lastmod", ""))
        for name, rows, metadata in raven_sources(self.root):
            self.provenance["raven_" + name] = metadata
            for raw in rows:
                duration = numeric(raw.get("bls_tdur"))
                self._add("raven", "RAVEN TIC " + raw["ticid"] + "." + raw.get("candidate", ""),
                          raw["ticid"], raw.get("ra"), raw.get("dec"), raw.get("bls_per"),
                          raw.get("bls_t0"), duration / 24 if duration is not None else None,
                          transiting=True, disposition="previous public RAVEN report: " + name,
                          reference="https://zenodo.org/records/19661443; arXiv:2603.22597",
                          release_table=name)
        # The user-supplied manuscript is a known prior report even though the
        # target is absent from the bulk online catalogs. Never rediscover the
        # recovery control by treating a catalog nonmatch as novelty.
        source = self.root / "candidates/tic4206066/recovery.json"
        manuscript = json.loads(source.read_text())["manuscript"]
        self.provenance["user_supplied_manuscript_tic4206066"] = dict(
            file=str(source), sha256=sha(source.read_bytes()), manuscript=manuscript,
            scope="User-supplied 2026-09-25 draft; reported candidate and explicitly tentative second signal, not confirmation.")
        for period, epoch, duration, label in [
                (3.182785, 3998.736, 1.87 / 24, "planet candidate"),
                (11.13274, 4000.865, None, "tentative post hoc signal")]:
            self._add("supplied_report", "Rabtsevich TIC 4206066 " + str(period), 4206066,
                      None, None, period, epoch, duration, transiting=True,
                      disposition=label, reference=manuscript["title"] + "; draft 2026-09-25",
                      manuscript_sha256=manuscript["sha256"])

    def match_target(self, tic, ra=None, dec=None):
        result = super().match_target(tic, ra, dec)
        result["is_raven_host"] = any(r["host_match"] and r["catalog"] == "raven" for r in result["matches"])
        result["is_supplied_report_host"] = any(r["host_match"] and r["catalog"] == "supplied_report" for r in result["matches"])
        return result

    def match_signal(self, tic, period_days, ra=None, dec=None, **kwargs):
        result = super().match_signal(tic, period_days, ra=ra, dec=dec, **kwargs)
        host = [r for r in result["matches"] if r["host_match"]]
        result["recognition_status"] = (
            "previously_reported_period" if result["known_period_match"] else
            "harmonic_of_previously_reported_period" if result["known_harmonic_match"] else
            "known_host_different_period_pending_signal_and_literature_checks" if host else
            "nearby_reported_source_pending_localization" if result["possible_nearby_period_contamination"] else
            "no_report_in_checked_snapshots_pending_literature")
        result["catalog_snapshot"] = str(self.snapshot.relative_to(self.root))
        result["novelty_status"] = "not_established; finite_catalog_scope_only"
        return result

    def metadata(self):
        result = super().metadata()
        result["fresh_snapshot"] = str(self.snapshot.relative_to(self.root))
        result["limitations"] += " RAVEN release membership includes unvalidated detections. No catalog absence is a discovery or exhaustive-priority result."
        return result


def rank_unknown_signals(frame, limit=10):
    """Operational follow-up queue, not statistical significance or discovery.

    Require positive, shallow P>1-day events and recurrence across two median
    widths in the SAME underlying light curve. Width agreement is a robustness
    check; it is not independent data or an empirical false-alarm probability.
    Known hosts with a different period remain eligible for additional planets.
    """
    import pandas as pd
    eligible = frame[
        (frame["period_days"] > 1) & (frame["local_depth_ppm"] > 0) &
        (frame["local_depth_ppm"] < 5000) & (frame["n_positive_events"] >= 3) &
        (frame["local_snr"] >= 6) & (frame["event_snr_max_fraction"] < .6) &
        ~frame["known_period_match"].astype(bool) & ~frame["known_harmonic_match"].astype(bool)
    ].copy()
    records = []
    for (_, _, _), group in eligible.groupby(["tic", "sector", "label"]):
        for _, row in group.iterrows():
            matching = group[abs(group["period_days"] / row["period_days"] - 1) < .01]
            count = matching["width_days"].nunique()
            if count < 2:
                continue
            records.append({**row.to_dict(), "matching_detrend_widths": int(count),
                            "minimum_width_local_snr": float(matching.groupby("width_days")["local_snr"].max().min()),
                            "priority_scope": "Numerical follow-up priority only; no calibrated significance or novelty."})
    if not records:
        return pd.DataFrame(columns=list(frame.columns) + ["matching_detrend_widths", "minimum_width_local_snr", "priority_scope"])
    return pd.DataFrame(records).sort_values(
        ["matching_detrend_widths", "minimum_width_local_snr", "local_snr"], ascending=False
    ).drop_duplicates("tic").head(limit)


def run_controls(index):
    """Actual prior reports plus contamination, mask, and service-error checks."""
    import numpy as np
    checks = []

    def check(name, condition, detail):
        if not condition:
            raise AssertionError(name)
        checks.append(dict(check=name, passed=True, detail=detail))

    r = index.match_signal(55525572, 83.8911)
    check("confirmed_planet_direct_period", r["known_period_match"] and r["is_confirmed_host"], r["recognition_status"])
    r = index.match_signal(71316387, 2*.6518369465110034)
    check("ctoi_and_raven_harmonic", r["known_harmonic_match"] and r["is_ctoi_host"] and r["is_raven_host"], r["recognition_status"])
    r = index.match_signal(261136679, 124.64)
    check("rv_planet_included", r["known_period_match"] and r["known_rv_period_match"], r["recognition_status"])
    for period in [3.182785, 11.13274]:
        r = index.match_signal(4206066, period)
        check("supplied_manuscript_" + str(period), r["known_period_match"] and r["is_supplied_report_host"], r["recognition_status"])
    known = next(row for row in index.rows if row["catalog"] == "confirmed" and row["name"] == "TOI-813 b")
    r = index.match_signal(99999999999, 83.8911, known["ra"], known["dec"] + 10/3600)
    check("nearby_period_not_assigned_to_target", r["possible_nearby_period_contamination"] and not r["known_period_match"], r["recognition_status"])
    r = index.match_signal(99999999999, 3.31)
    check("absence_not_novelty", r["recognition_status"] == "no_report_in_checked_snapshots_pending_literature" and "not_established" in r["novelty_status"], r["recognition_status"])
    r = index.match_signal(55525572, 7.321)
    check("new_period_on_known_host_retained", "known_host_different_period" in r["recognition_status"], r["recognition_status"])
    keep = index.transit_mask(4206066, np.array([3998.736, 3998.736 + 3.182785/2]))
    check("transit_mask_actual_epoch", not bool(keep[0]) and bool(keep[1]), keep.tolist())
    for label, body in [("html_error", b"<html>error</html>"),
                        ("schema_error", b"wrong,column\n1,2\n"), ("short_error", b"tic,period\n")]:
        try:
            parse_csv(body, expected_columns=["tic", "period"])
        except ValueError:
            check(label, True, "Invalid response rejected; not treated as an empty catalog")
        else:
            raise AssertionError(label)
    rows = parse_vsx_tsv((ROOT / "data/new_stars/rv_recognition/4690783166942834816_vsx.tsv").read_bytes())
    check("vsx_prior_rotator_recovered", any(row["Name"] == "ASASSN-V J005926.01-712933.7" and abs(float(row["Period"]) - 20.766) < 1e-5 for row in rows), "Public named ROT source recovered at its published 20.766-day period")
    check("valid_vsx_no_row_response", parse_vsx_tsv((ROOT / "candidates/tic71316387/catalog/vsx_cone.tsv").read_bytes()) == [], "Complete table/column metadata required before accepting zero rows")
    try:
        parse_vsx_tsv(b"#INFO\tError=Table or Catalog not found\n")
    except ValueError:
        check("failed_vsx_not_an_absence", True, "ASU error rejected")
    else:
        raise AssertionError("failed_vsx_not_an_absence")
    result = dict(checked_utc=utc(), checks=checks, index=index.metadata(),
                  method_sha256=sha(Path(__file__).read_bytes()),
                  base_method_sha256=sha((ROOT / "src/broad_crossmatch.py").read_bytes()))
    write_json(RESULTS / "positive_controls.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retrieve", action="store_true")
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--tic", type=int)
    parser.add_argument("--period", type=float)
    parser.add_argument("--live-tic", type=int, action="append",
                        help="Lead-only live ExoFOP/alias and finite literature check; repeat for multiple TICs")
    parser.add_argument("--rank-unknown", type=Path, help="Write numerical top-10 follow-up queue from an annotated input")
    parser.add_argument("--controls", action="store_true", help="Run actual recognition and error-response controls")
    parser.add_argument("--deep-tic", type=int, action="append", help="Deeper Gaia/ASAS-SN/SIMBAD checks after --live-tic")
    args = parser.parse_args()
    snapshot = retrieve_snapshot() if args.retrieve else args.snapshot
    index = FreshCatalogIndex(snapshot=snapshot)
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_json(RESULTS / "catalog_index.json", index.metadata())
    if args.controls:
        print(json.dumps(dict(passed_controls=len(run_controls(index)["checks"])) ))
    if args.live_tic:
        print(json.dumps(dict(live_context=str(capture_live_context(args.live_tic)))))
        summarize_live_context()
    if args.deep_tic:
        print(json.dumps(dict(deep_context=str(capture_deep_context(args.deep_tic)))))
    if args.input:
        if not args.output or args.input.resolve() == args.output.resolve():
            parser.error("A distinct --output is required")
        import pandas as pd
        input_body = args.input.read_bytes()
        input_table = gzip.decompress(input_body) if args.input.suffix == ".gz" else input_body
        frame = pd.read_csv(io.BytesIO(input_table))
        result = index.annotate_dataframe(frame)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, index=False)
        saved_input = RESULTS / "inputs" / (args.input.stem + "_" + sha(input_table)[:16] + ".csv.gz")
        saved_input.parent.mkdir(parents=True, exist_ok=True)
        if not saved_input.exists():
            saved_input.write_bytes(gzip.compress(input_table, mtime=0))
        elif sha(gzip.decompress(saved_input.read_bytes())) != sha(input_table):
            raise ValueError("Changed retained exact input table")
        write_json(args.output.with_suffix(args.output.suffix + ".manifest.json"), dict(
            created_utc=utc(), input_path=str(args.input), input_file_sha256=sha(input_body),
            input_table_sha256=sha(input_table), retained_exact_input=str(saved_input.relative_to(ROOT)),
            output_sha256=sha(args.output.read_bytes()), row_count=len(result),
            method_sha256=sha(Path(__file__).read_bytes()),
            base_method_sha256=sha((ROOT / "src/broad_crossmatch.py").read_bytes()),
            catalogs=index.metadata()))
        if args.rank_unknown:
            args.rank_unknown.parent.mkdir(parents=True, exist_ok=True)
            rank_unknown_signals(result).to_csv(args.rank_unknown, index=False)
        print(json.dumps(dict(rows=len(result), output=str(args.output))))
    elif args.tic:
        print(json.dumps(index.match_signal(args.tic, args.period), indent=2))
    else:
        print(json.dumps(dict(snapshot=str(index.snapshot), entries=len(index.rows),
                              known_tics=len(index.known_tics))))


if __name__ == "__main__":
    main()
