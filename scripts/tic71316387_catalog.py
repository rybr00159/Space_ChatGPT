"""Retrieve immutable public provenance for the TIC 71316387 claim audit.

Existing snapshots are checksum-verified and reused unless --refresh is given.
No authenticated services or paid API calls are used. Large paper files remain
in TESS_DATA_DIR (default /workspace/tess-research-data), outside ordinary Git.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import os
import subprocess
from urllib.parse import urlencode

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "candidates/tic71316387/catalog"
CACHE = Path(os.environ.get("TESS_DATA_DIR", "/workspace/tess-research-data")) / "tic71316387/catalog"
SID = "6571343840265900544"
RA, DEC = 330.887940458872, -42.8355109237273


def utc():
    return datetime.now(timezone.utc).isoformat()


def fetch(name, url, refresh=False, large=False):
    directory = CACHE if large else OUT
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    meta_path = OUT / (name + ".manifest.json")
    if not refresh:
        # A prior HTTP failure may have been followed by a successful immutable
        # retry snapshot. Reuse it instead of repeating a completed download.
        for candidate in sorted(OUT.glob(Path(name).stem + "*.manifest.json"), reverse=True):
            old = json.loads(candidate.read_text())
            local = Path(old.get("local_path", str(directory / candidate.name.removesuffix(".manifest.json"))))
            if old.get("url") == url and old.get("http_status") == 200 and local.exists():
                if hashlib.sha256(local.read_bytes()).hexdigest() != old["sha256"]:
                    raise RuntimeError(f"Snapshot checksum mismatch: {local}")
                return local
    if path.exists() and meta_path.exists() and not refresh:
        old = json.loads(meta_path.read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest() != old["sha256"]:
            raise RuntimeError(f"Snapshot checksum mismatch: {path}")
        if old["http_status"] == 200:
            return path
        refresh = True  # Preserve unsuccessful responses before a retry.
    if refresh or path.exists():
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        path = path.with_name(path.stem + "_" + stamp + path.suffix)
        meta_path = OUT / (path.name + ".manifest.json")
    response = requests.get(url, timeout=(20, 120), headers={"User-Agent": "Space_ChatGPT-public-research/0.5"})
    body = response.content
    path.write_bytes(body)
    meta = dict(url=url, retrieved_utc=utc(), http_status=response.status_code,
                bytes=len(body), sha256=hashlib.sha256(body).hexdigest(),
                md5=hashlib.md5(body).hexdigest(), local_path=str(path),
                content_type=response.headers.get("Content-Type"))
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    response.raise_for_status()
    return path


def tap(name, query, refresh, service="gaia"):
    base = {"gaia": "https://gea.esac.esa.int/tap-server/tap/sync",
            "simbad": "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"}[service]
    url = base + "?" + urlencode(dict(request="doQuery", lang="ADQL", format="csv", query=query))
    return fetch(name, url, refresh)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--large-papers", action="store_true", help="Retain ~55MB paper/source files outside Git.")
    args = parser.parse_args()
    urls = {
        "exofop.json": "https://exofop.ipac.caltech.edu/tess/target.php?id=71316387&json",
        "target_page.html": "https://exofop.ipac.caltech.edu/tess/target.php?id=71316387",
        "zenodo.json": "https://zenodo.org/api/records/19661443",
        "raven_summary.pdf": "https://exofop.ipac.caltech.edu/tess/get_file.php?id=1113786",
        "arxiv_raven.xml": "https://export.arxiv.org/api/query?" + urlencode(dict(search_query="au:Lafarga AND all:RAVEN", max_results=10)),
        "raven_pipeline.html": "https://arxiv.org/html/2509.17645v1",
        "tesseb.html": "https://tessebs.villanova.edu/0071316387",
        "exofop_nearby.csv": "https://exofop.ipac.caltech.edu/tess/download_nearbytarget.php?id=71316387&output=csv",
        "vsx_cone.tsv": "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode({"-source": "B/vsx/vsx", "-c": f"{RA} {DEC}", "-c.rs": 120, "-out.max": "unlimited", "-out.all": "", "-out.form": "tsv"}),
        "arxiv_target_aliases.xml": "https://export.arxiv.org/api/query?" + urlencode(dict(search_query='all:"TIC 71316387" OR all:"71316387" OR all:"2MASS J22033305-4250076" OR all:"TYC 8001-154-1" OR all:"6571343840265900544"', max_results=100)),
    }
    paths, failures = {}, {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [(name, pool.submit(fetch, name, url, args.refresh)) for name, url in urls.items()]
        for name, future in futures:
            try:
                paths[name] = str(future.result())
            except Exception as exc:
                failures[name] = str(exc)
    if "zenodo.json" in paths:
        record = json.loads(Path(paths["zenodo.json"]).read_text())
        files = [x for x in record["files"] if x["size"] < 5_000_000]
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [(item, pool.submit(fetch, item["key"], item["links"]["self"], args.refresh)) for item in files]
            for item, future in futures:
                try:
                    path = future.result()
                    checksum = "md5:" + hashlib.md5(path.read_bytes()).hexdigest()
                    if checksum != item["checksum"]:
                        raise RuntimeError("Zenodo MD5 differs from record")
                    paths[item["key"]] = str(path)
                except Exception as exc:
                    failures[item["key"]] = str(exc)
    queries = [
        ("gaia_target.csv", f"select source_id,ra,dec,phot_g_mean_mag,bp_rp,ruwe,non_single_star,phot_variable_flag,parallax,parallax_error,astrometric_excess_noise,ipd_frac_multi_peak,radial_velocity,radial_velocity_error,rv_amplitude_robust,rv_chisq_pvalue from gaiadr3.gaia_source where source_id={SID}", "gaia"),
        ("gaia_eb.csv", f"select * from gaiadr3.vari_eclipsing_binary where source_id={SID}", "gaia"),
        ("gaia_class.csv", f"select * from gaiadr3.vari_classifier_result where source_id={SID}", "gaia"),
        ("gaia_neighbors.csv", f"select g.source_id,g.ra,g.dec,g.phot_g_mean_mag,g.bp_rp,g.ruwe,g.non_single_star,g.ipd_frac_multi_peak,v.in_vari_eclipsing_binary from gaiadr3.gaia_source g left outer join gaiadr3.vari_summary v on g.source_id=v.source_id where 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{RA},{DEC},0.0333333333333333)) order by g.phot_g_mean_mag", "gaia"),
        ("simbad_cone.csv", f"select oid,main_id,ra,dec,otype,nbref from basic where 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.0333333333333333))", "simbad"),
    ]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [(name, pool.submit(tap, name, query, args.refresh, service)) for name, query, service in queries]
        for name, future in futures:
            try:
                paths[name] = str(future.result())
            except Exception as exc:
                failures[name] = str(exc)
    if args.large_papers:
        for name, url in [("raven_paper.pdf", "https://arxiv.org/pdf/2603.22597v2"),
                          ("paper_source.tar", "https://arxiv.org/src/2603.22597v2")]:
            try:
                paths[name] = str(fetch(name, url, args.refresh, large=True))
            except Exception as exc:
                failures[name] = str(exc)
    for name in ["raven_summary.pdf", "raven_paper.pdf"]:
        if name in paths:
            dest = OUT / (Path(paths[name]).stem + ".txt")
            subprocess.run(["pdftotext", "-layout", paths[name], str(dest)], check=True)
            paths[dest.name] = str(dest)
    target_rows = {}
    for name in ["nsfp09_table.csv", "vet_table.csv", "val_table.csv", "val_rp8_table.csv"]:
        if name in paths:
            target_rows[name] = [row for row in csv.DictReader(Path(paths[name]).open()) if row.get("ticid") == "71316387"]
    summary = dict(checked_utc=utc(), tic=71316387, gaia_source_id=SID,
                   paths=paths, failures=failures, raven_target_rows=target_rows,
                   scope="Target catalog provenance and finite public checks; catalog absence does not establish novelty or rule out binarity.")
    summary_path = OUT / "catalog_check.json"
    if summary_path.exists():
        summary_path = OUT / ("catalog_check_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + ".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(dict(summary=str(summary_path), failures=failures), indent=2))


if __name__ == "__main__":
    main()
