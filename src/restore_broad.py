"""Verify or restore exact QLP products recorded in a gzip JSONL manifest.

Default mode uses no network. Downloads require both --download and a positive
--max-downloads cap. Repeated runs verify and skip retained products, never
overwrite a hash mismatch, and leave failed downloads available for inspection.
"""
from __future__ import annotations

from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import threading
import time
import uuid
from urllib.parse import urlparse

import numpy as np
import requests
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
PRODUCT = re.compile(r"hlsp_qlp_tess_ffi_s(\d{4})-(\d{16})_tess_v\d+_llc\.fits")
_local = threading.local()


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ones_complement(block):
    """Raw FITS sum includes padding and remains valid without DATASUM cards."""
    if len(block) % 4:
        raise ValueError("FITS HDU is not aligned to four-byte words")
    total = int(np.frombuffer(block, dtype=">u4").sum(dtype=np.uint64))
    while total >> 32:
        total = (total & 0xffffffff) + (total >> 32)
    return total


def inspect_product(path, record):
    """Accept a product only after expected byte hash and raw HDU validation."""
    size = path.stat().st_size
    if size != record["raw_bytes"]:
        raise ValueError(f"Expected {record['raw_bytes']} bytes, found {size}")
    digest = digest_file(path)
    if digest != record["sha256"]:
        raise ValueError(f"SHA256 differs from recorded observation: {digest}")
    raw = path.read_bytes()
    checks = []
    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        for index, hdu in enumerate(hdus):
            info = hdu.fileinfo()
            start, data = info["hdrLoc"], info["datLoc"]
            stop = data + info["datSpan"]
            item = {"hdu": index, "checksum": "absent", "datasum": "absent"}
            if "CHECKSUM" in hdu.header:
                if ones_complement(raw[start:stop]) != 0xffffffff:
                    raise ValueError(f"Raw CHECKSUM failed in HDU {index}")
                item["checksum"] = "verified_raw_ones_complement"
            if "DATASUM" in hdu.header:
                if ones_complement(raw[data:stop]) != int(hdu.header["DATASUM"]):
                    raise ValueError(f"Raw DATASUM failed in HDU {index}")
                item["datasum"] = "verified_raw_ones_complement"
            checks.append(item)
        header = hdus[1].header
        for key, expected in (("TICID", record["tic"]), ("SECTOR", record["sector"])):
            actual = header.get(key, hdus[0].header.get(key))
            if actual is None or int(actual) != expected:
                raise ValueError(f"FITS {key} differs from recorded product")
        if header.get("TIMESYS", hdus[0].header.get("TIMESYS")) != "TDB":
            raise ValueError("FITS time system is not TDB")
        reference = header.get("BJDREFI", 0) + header.get("BJDREFF", 0)
        if abs(reference - 2457000) > 1e-8:
            raise ValueError("FITS timing is not the expected BTJD reference")
    return {"sha256": digest, "raw_bytes": size, "fits_checksums": checks}


def normalize_record(record, cache):
    """Map recorded absolute paths to a portable, deterministic cache layout."""
    tic = int(record["tic"])
    sha = record["sha256"]
    size = int(record["raw_bytes"])
    if not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{64}", sha):
        raise ValueError("Expected SHA256 must be a lowercase 64-character digest")
    if size <= 0:
        raise ValueError("Expected raw byte count must be positive")
    url = urlparse(record["url"])
    if (url.scheme != "https" or url.netloc != "archive.stsci.edu" or
            not url.path.startswith("/hlsps/qlp/") or url.query or url.fragment):
        raise ValueError("Only recorded public HTTPS MAST QLP product URLs are supported")
    name = Path(url.path).name
    match = PRODUCT.fullmatch(name)
    if match is None or int(match[2]) != tic:
        raise ValueError("URL product name does not match TIC")
    sector = int(match[1])
    metadata = record.get("metadata", {})
    if int(metadata.get("SECTOR", sector)) != sector:
        raise ValueError("URL sector differs from manifest metadata")
    if int(metadata.get("TICID", tic)) != tic:
        raise ValueError("TIC differs from manifest metadata")
    path = cache / "broad" / f"qlp_s{sector:04d}" / "raw" / f"{tic % 1000:03d}" / name
    return {**record, "tic": tic, "sector": sector, "raw_bytes": size, "destination": path}


class DownloadBudget:
    """Bound attempts and concurrent disk reservations before opening files."""
    def __init__(self, cache, maximum, reserve):
        self.cache = cache
        self.maximum = maximum
        self.reserve = reserve
        self.attempts = 0
        self.in_flight_bytes = 0
        self.lock = threading.Lock()

    def claim(self, size):
        with self.lock:
            if self.attempts >= self.maximum:
                return "missing_download_limit"
            available = shutil.disk_usage(self.cache).free
            if available - self.in_flight_bytes - size < self.reserve:
                return "missing_storage_reserve"
            self.attempts += 1
            self.in_flight_bytes += size
            return None

    def release(self, size):
        with self.lock:
            self.in_flight_bytes -= size


def restore_one(record, download, budget, run_id):
    path = record["destination"]
    result = {"tic": record["tic"], "sector": record["sector"],
              "url": record["url"], "expected_sha256": record["sha256"],
              "recorded_path": record.get("path"), "path": str(path), "checked_utc": utc()}
    if path.exists():
        try:
            result.update(inspect_product(path, record), status="verified_cached")
        except Exception as error:
            result.update(status="existing_integrity_failure", error=f"{type(error).__name__}: {error}")
        return result
    if not download:
        return {**result, "status": "missing"}
    blocked = budget.claim(record["raw_bytes"])
    if blocked:
        return {**result, "status": blocked}
    temporary = path.with_name(path.name + f".restore-{run_id}.part")
    result["temporary_path"] = str(temporary)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        if not hasattr(_local, "session"):
            _local.session = requests.Session()
            _local.session.headers.update({"User-Agent": "Space_ChatGPT-public-TESS-restoration/1"})
        # Individual retries are explicit in the evidence log; no partial file is
        # silently removed. Re-running resumes at the level of verified files.
        with _local.session.get(record["url"], stream=True, timeout=(15, 60)) as response:
            response.raise_for_status()
            length = response.headers.get("Content-Length")
            if length is not None and int(length) != record["raw_bytes"]:
                raise ValueError("Archive Content-Length differs from recorded observation")
            written = 0
            with temporary.open("xb") as stream:
                for block in response.iter_content(65536):
                    if not block:
                        continue
                    written += len(block)
                    if written > record["raw_bytes"]:
                        raise ValueError("Archive response exceeds recorded product size")
                    stream.write(block)
                stream.flush()
                os.fsync(stream.fileno())
        integrity = inspect_product(temporary, record)
        # Atomic hard linking will fail if another worker/run has created the
        # destination. Unlike replace(), it cannot overwrite retained evidence.
        try:
            os.link(temporary, path)
            temporary.unlink()  # Verified bytes remain at their final location.
        except FileExistsError:
            inspect_product(path, record)
            result["duplicate_verified_download_retained"] = True
        result.update(integrity, status="downloaded_verified")
    except Exception as error:
        result.update(status="download_failure", error=f"{type(error).__name__}: {error}")
    finally:
        budget.release(record["raw_bytes"])
    return result


def manifest_rows(stream, tics, limit):
    seen = set()
    selected = 0
    with gzip.GzipFile(fileobj=stream, mode="rb") as zipped:
        for line_number, line in enumerate(zipped, 1):
            record = json.loads(line)
            tic = int(record["tic"])
            if tics and tic not in tics:
                continue
            # Multiple sectors are supported; repeated identical product URLs
            # indicate a malformed export rather than an extra observation.
            key = record["url"]
            if key in seen:
                raise ValueError(f"Duplicate product URL at manifest line {line_number}")
            seen.add(key)
            yield line_number, record
            selected += 1
            if limit and selected >= limit:
                break


def run(args):
    cache = args.cache.resolve()
    if args.download:
        cache.mkdir(parents=True, exist_ok=True)
    if not cache.exists():
        # Verification does not require constructing a new raw-data cache.
        if args.output_dir is None:
            raise ValueError("Missing cache: provide --output-dir for verification logs")
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + uuid.uuid4().hex[:8]
    output = args.output_dir or cache / "restoration"
    output.mkdir(parents=True, exist_ok=True)
    evidence = output / f"{run_id}.jsonl"
    summary_path = output / f"{run_id}.summary.json"
    counts = Counter()
    processed_tics = set()
    failures = []
    started = time.monotonic()
    budget = DownloadBudget(cache, args.max_downloads or 0, int(args.reserve_gib * 1024**3))
    summary = {"started_utc": utc(), "mode": "limited_download" if args.download else "verify",
               "manifest": str(args.manifest.resolve()), "cache": str(cache),
               "evidence_log": str(evidence), "max_downloads": args.max_downloads,
               "workers": args.workers, "reserve_gib": args.reserve_gib,
               "selection_limit": args.limit, "selected_tics": sorted(args.tic or []),
               "source_sha256": digest_file(Path(__file__)), "credentials_required": False}
    error = None
    try:
        # Hash and read the same open descriptor: a concurrent atomic survey
        # export cannot change the manifest snapshot underlying this run.
        with args.manifest.open("rb") as raw, evidence.open("x") as log:
            digest = hashlib.sha256()
            for block in iter(lambda: raw.read(1024 * 1024), b""):
                digest.update(block)
            summary["manifest_sha256"] = digest.hexdigest()
            raw.seek(0)
            rows = iter(manifest_rows(raw, set(args.tic or []), args.limit))
            exhausted = False
            jobs = {}
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                while not exhausted or jobs:
                    while not exhausted and len(jobs) < args.workers * 2:
                        try:
                            line_number, record = next(rows)
                        except StopIteration:
                            exhausted = True
                            break
                        normalized = normalize_record(record, cache)
                        future = executor.submit(restore_one, normalized, args.download, budget, run_id)
                        jobs[future] = line_number
                    if not jobs:
                        break
                    ready, _ = wait(jobs, return_when=FIRST_COMPLETED)
                    for future in ready:
                        item = future.result()
                        item["manifest_line"] = jobs.pop(future)
                        processed_tics.add(item["tic"])
                        counts[item["status"]] += 1
                        if item.get("error") and len(failures) < 20:
                            failures.append(item)
                        log.write(json.dumps(item, allow_nan=False) + "\n")
                    log.flush()
    except (Exception, KeyboardInterrupt) as caught:
        error = f"{type(caught).__name__}: {caught}"
    summary.update(finished_utc=utc(), elapsed_seconds=round(time.monotonic() - started, 3),
                   status_counts=dict(counts), processed_records=sum(counts.values()),
                   download_attempts=budget.attempts, error=error, first_failures=failures,
                   requested_tics_not_processed=sorted(set(args.tic or []) - processed_tics))
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"summary": str(summary_path), "counts": dict(counts),
                      "download_attempts": budget.attempts, "error": error}, allow_nan=False))
    if error or failures:
        return 1
    if (not counts or summary["requested_tics_not_processed"] or
            any(name.startswith("missing") for name in counts)):
        return 2
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "results/broad/observations.jsonl.gz")
    parser.add_argument("--cache", type=Path,
                        default=Path(os.environ.get("TESS_DATA_DIR", "/workspace/tess-research-data")))
    parser.add_argument("--download", action="store_true", help="Allow capped missing-product downloads")
    parser.add_argument("--max-downloads", type=int, help="Required positive cap with --download; counts attempts")
    parser.add_argument("--workers", type=int, default=4, help="Bounded parallel workers, 1 to 16")
    parser.add_argument("--reserve-gib", type=float, default=5.0, help="Free space retained before each download")
    parser.add_argument("--limit", type=int, default=0, help="Maximum selected manifest records; 0 means all")
    parser.add_argument("--tic", type=int, nargs="+", help="Restrict to specified TIC identifiers")
    parser.add_argument("--output-dir", type=Path, help="Directory for append-only per-run verification evidence")
    args = parser.parse_args()
    if not 1 <= args.workers <= 16:
        parser.error("--workers must be between 1 and 16")
    if args.limit < 0 or not math.isfinite(args.reserve_gib) or args.reserve_gib < 0:
        parser.error("--limit and --reserve-gib must be finite nonnegative values")
    if args.download and (args.max_downloads is None or args.max_downloads <= 0):
        parser.error("--download requires an explicit positive --max-downloads cap")
    if not args.download and args.max_downloads is not None:
        parser.error("--max-downloads requires --download")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
