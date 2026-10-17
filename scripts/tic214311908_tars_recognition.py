"""Retrieve independent light curves and check TARS' published default catalog.

This is a catalog-recognition audit, not an independent exoplanet validation.
Raw files are cached outside Git. Public URLs and checksums are retained.
"""
from __future__ import annotations

import concurrent.futures
import csv
import datetime as dt
import hashlib
import io
import json
from pathlib import Path
import zipfile

import requests
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/tic214311908/tars_recognition"
CACHE = Path("/workspace/tess-research-data/tic214311908/independent_reductions")
TIC = 214311908


def fetch(url: str, destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        with requests.get(url, timeout=(30, 240), stream=True) as response:
            response.raise_for_status()
            temporary = destination.with_suffix(destination.suffix + ".part")
            with temporary.open("wb") as handle:
                for block in response.iter_content(1024 * 1024):
                    handle.write(block)
            temporary.replace(destination)
    digest = hashlib.sha256()
    md5 = hashlib.md5()
    with destination.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
            md5.update(block)
    return {"url": url, "path": str(destination), "bytes": destination.stat().st_size,
            "sha256": digest.hexdigest(), "md5": md5.hexdigest(),
            "retrieved_or_verified_utc": dt.datetime.now(dt.timezone.utc).isoformat()}


def download_lightcurve(record: dict) -> dict:
    uri = record["dataURL"]
    url = requests.Request("GET", "https://mast.stsci.edu/api/v0.1/Download/file",
                           params={"uri": uri}).prepare().url
    destination = CACHE / uri.rsplit("/", 1)[-1]
    metadata = fetch(url, destination)
    metadata.update({"mast_uri": uri, "provenance_name": record["provenance_name"],
                     "sector": record["sequence_number"]})
    headers = []
    with fits.open(destination) as hdus:
        for i, hdu in enumerate(hdus):
            headers.append({"hdu": i, "name": hdu.name,
                            "header": hdu.header.tostring(sep="\n", endcard=True, padding=False),
                            "columns": list(hdu.columns.names) if hasattr(hdu, "columns") else [],
                            "rows": len(hdu.data) if hdu.data is not None else 0})
        metadata["hdus"] = [{k: v for k, v in h.items() if k != "header"} for h in headers]
        metadata["reported_period_metadata"] = {
            str(i): {k: hdus[i].header[k] for k in hdus[i].header
                     if k.startswith(("PER", "POW", "SNR", "P_SYS", "P_ALIAS", "PSYS"))
                     or k in ("APERTURE", "AMPFIT", "P2PRMS", "REDCHI2", "BIC1T", "BIC2T", "BIC2TS")}
            for i in range(len(hdus))
        }
    (OUT / (destination.name + ".headers.json")).write_text(json.dumps(headers, indent=2))
    (OUT / (destination.name + ".manifest.json")).write_text(json.dumps(metadata, indent=2))
    return metadata


def default_catalog() -> dict:
    record = json.loads((OUT / "zenodo18342590.json").read_text())
    item = next(x for x in record["files"] if x["key"] == "tars_table_2.csv.zip")
    manifest = fetch(item["links"]["self"], CACHE / item["key"])
    manifest["zenodo_checksum"] = item.get("checksum")
    if item.get("checksum", "").startswith("md5:"):
        assert manifest["md5"] == item["checksum"].split(":", 1)[1]
    # TARS headers store Gaia identifiers as floats, which cannot retain the
    # precision of 19-digit source IDs. Use the exact integer from raw Gaia CSV.
    gaia_source = ROOT / "data/tic214311908/catalogs/gaia_target.csv"
    with gaia_source.open() as handle:
        gaia_ids = {row["source_id"] for row in csv.DictReader(handle)}
    # Also check the rounded value displayed in the TARS preview, separately
    # documented as an identifier precision issue rather than a true Gaia ID.
    gaia_ids.add("5925373272850784256")
    result = {"tic": TIC, "gaia_identifiers_checked": sorted(gaia_ids),
              "gaia_identifier_source": {"path": str(gaia_source),
                  "sha256": hashlib.sha256(gaia_source.read_bytes()).hexdigest(),
                  "note": "Exact archive source ID plus rounded TARS preview value; TARS FITS Gaia IDs are stored as imprecise floats."},
              "matches": [], "scanned_rows": 0, "source": manifest,
              "note": "The default TARS rotation/variability catalog is not a planet candidate catalog."}
    with zipfile.ZipFile(manifest["path"]) as archive:
        result["archive_members"] = archive.namelist()
        for name in archive.namelist():
            if not name.endswith(".csv"):
                continue
            with archive.open(name) as binary:
                reader = csv.DictReader(io.TextIOWrapper(binary))
                result["fieldnames"] = reader.fieldnames
                tic_columns = [k for k in reader.fieldnames if "tic" in k.lower()]
                gaia_columns = [k for k in reader.fieldnames if k.lower() in ("dr2_source_id", "dr3_source_id")]
                assert tic_columns, reader.fieldnames
                for row in reader:
                    result["scanned_rows"] += 1
                    tic_match = any(row[k].strip() in (str(TIC), str(TIC) + ".0") for k in tic_columns)
                    gaia_match = any(row[k].strip() in gaia_ids for k in gaia_columns)
                    if tic_match or gaia_match:
                        result["matches"].append(row)
    (OUT / "default_catalog_check.json").write_text(json.dumps(result, indent=2))
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inventory = json.loads(Path("/workspace/tess-research-data/tic214311908/queries/mast_all_collections.json").read_text())
    records = [r for r in inventory["data"] if r["provenance_name"] in ("TARS", "GSFC-ELEANOR-LITE")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        results = list(executor.map(download_lightcurve, records))
    (OUT / "independent_reduction_manifest.json").write_text(json.dumps(results, indent=2))
    print(json.dumps([{k: x[k] for k in ("sector", "provenance_name", "path", "hdus", "reported_period_metadata")} for x in results], indent=2))
    result = default_catalog()
    print("Default catalog", result["scanned_rows"], "rows; target matches:", json.dumps(result["matches"]))


if __name__ == "__main__":
    main()
