"""Fast, provenance-aware matching of broad transit screens to public catalogs.

Catalog absence is never evidence of discovery. Position matches within 30 arcsec
are contamination warnings; only exact TIC matches or <=2 arcsec matches count
as host associations. Confirmed RV planets are retained even when no ephemeris
exists. No input table is modified by this module.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]

SUPPORT_URLS = {
    "diamante_v2.csv": "https://archive.stsci.edu/hlsps/diamante/hlsp_diamante_tess_lightcurve_catalog_tess_v2_cat.csv",
    "diamante_v2_readme.txt": "https://archive.stsci.edu/hlsps/diamante/hlsp_diamante_tess_lightcurve_all_tess_v2_readme-catalog.txt",
    **{name: "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + urlencode(
        {"-source": source, "-out.max": "unlimited", "-out.all": "", "-out.form": "tsv"})
       for name, source in [("tesseb_prsa2022.tsv", "J/ApJS/258/16"),
                            ("olmschenk2021.tsv", "J/AJ/161/273"),
                            ("yu2019.tsv", "J/AJ/158/25")]},
}


def _vizier_rows(path):
    lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    if len(lines) < 4 or "TIC" not in lines[0].split("\t"):
        raise ValueError(f"Missing expected TIC table in {path}")
    # ASU-TSV has one units row and one hyphen separator after the header.
    rows = list(csv.DictReader([lines[0]] + lines[3:], delimiter="\t"))
    return [{k: v.strip() if isinstance(v, str) else v for k, v in r.items()}
            for r in rows if str(r.get("TIC", "")).strip().isdigit()]


def _write_derived(path, rows, columns, sources, scope):
    if path.exists():
        # Preserve an existing snapshot, and validate before reusing it.
        _read_verified(path)
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    metadata = {"derived_utc": datetime.now(timezone.utc).isoformat(), "row_count": len(rows),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "scope": scope,
                "sources": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    path.with_suffix(path.suffix + ".manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")


def retrieve_support(root=ROOT):
    """Restore exact cached snapshots from manifests; preserve existing raw files.

    TLS verification remains enabled. Successful retrieval hashes are recorded;
    an existing hash is checked before reuse. Original NASA/ExoFOP tables are
    handled separately by data/catalogs/retrieve_catalogs.py.
    """
    destination = Path(root) / "data" / "broad_catalogs"
    destination.mkdir(parents=True, exist_ok=True)
    for name, url in SUPPORT_URLS.items():
        path = destination / name
        manifest_path = path.with_suffix(path.suffix + ".manifest.json")
        if path.exists():
            metadata = json.loads(manifest_path.read_text())
            if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["sha256"]:
                raise ValueError(f"Checksum mismatch in cached {path}; preserve and investigate")
            continue
        metadata = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        request_url = metadata.get("url", url)
        with urlopen(Request(request_url, headers={"User-Agent": "Space_ChatGPT-public-research/0.2"}), timeout=120) as response:
            body = response.read()
            status = response.status
        digest = hashlib.sha256(body).hexdigest()
        if metadata.get("sha256") and metadata["sha256"] != digest:
            raise ValueError(f"Remote catalog changed since snapshot {name}; do not silently replace provenance")
        path.write_bytes(body)
        if not metadata:
            manifest_path.write_text(json.dumps({"retrieved_utc": datetime.now(timezone.utc).isoformat(),
                                                "url": request_url, "http_status": status,
                                                "bytes": len(body), "sha256": digest}, indent=2) + "\n")
    eb_source = destination / "tesseb_prsa2022.tsv"
    eb_rows = [{"tic": int(r["TIC"]), "in_catalog": "True", "ra": r["RAJ2000"],
                "dec": r["DEJ2000"], "epoch_btjd": r["BJD0"], "period_days": r["Per"],
                "morphology": r["Morph"], "signal_id": r["m_TIC"],
                "reference": "2022ApJS..258...16P; CDS J/ApJS/258/16; sectors 1-26"}
               for r in _vizier_rows(eb_source)]
    _write_derived(destination / "tesseb_members.csv", eb_rows,
                   ["tic", "in_catalog", "ra", "dec", "epoch_btjd", "period_days", "morphology", "signal_id", "reference"],
                   [eb_source], "Published EB membership table for sectors 1-26; not the complete current live catalog.")
    published = []
    for name, label, period_column, epoch_column, reference in [
            ("olmschenk2021.tsv", "Olmschenk 2021", "Porb", "BJD", "2021AJ....161..273O; CDS J/AJ/161/273"),
            ("yu2019.tsv", "Yu 2019", "Per", "T0", "2019AJ....158...25Y; CDS J/AJ/158/25")]:
        for r in _vizier_rows(destination / name):
            duration = numeric(r.get("Dur"))
            published.append({"tic": int(r["TIC"]), "name": label + " TIC " + r["TIC"],
                              "ra": r.get("_RA", ""), "dec": r.get("_DE", ""),
                              "period_days": r[period_column], "epoch_btjd": r[epoch_column],
                              "duration_days": duration / 24 if duration is not None else "",
                              "reference": reference})
    _write_derived(destination / "published_candidates.csv", published,
                   ["tic", "name", "ra", "dec", "period_days", "epoch_btjd", "duration_days", "reference"],
                   [destination / "olmschenk2021.tsv", destination / "yu2019.tsv"],
                   "Published vetted candidates only; no complete all-sector/all-TCE claim.")


def numeric(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def tic_number(value):
    if value is None:
        return None
    value = re.sub(r"^TIC\s*", "", str(value).strip(), flags=re.I)
    number = numeric(value)
    return int(number) if number is not None and number > 0 and number == int(number) else None


def _xyz(ra, dec):
    ra, dec = np.deg2rad([ra, dec])
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _read_verified(path, required=True):
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return [], None
    body = path.read_bytes()
    digest = hashlib.sha256(body).hexdigest()
    metadata_path = path.with_suffix(path.suffix + ".manifest.json")
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    if metadata.get("sha256") and metadata["sha256"] != digest:
        raise ValueError(f"Catalog checksum mismatch: {path}; preserve and investigate")
    metadata = {**metadata, "file": str(path), "sha256": digest}
    rows = list(csv.DictReader(body.decode("utf-8-sig").splitlines()))
    return rows, metadata


class CatalogIndex:
    """Load catalogs once, index TIC identifiers and unit-sphere coordinates.

    ``match_target`` and ``match_signal`` return JSON-serializable dictionaries.
    ``transit_mask`` returns True at samples that may be searched after masking
    complete known transit ephemerides. It does not mask RV or nearby-object
    entries, and never guesses an absent epoch.
    """

    def __init__(self, root=ROOT, position_radius_arcsec=30.0, host_radius_arcsec=2.0):
        self.root = Path(root)
        self.position_radius_arcsec = float(position_radius_arcsec)
        self.host_radius_arcsec = float(host_radius_arcsec)
        self.rows = []
        self.provenance = {}
        self.missing_optional = []
        self._load()
        self.by_tic = defaultdict(list)
        sky_xyz, sky_indices = [], []
        for index, row in enumerate(self.rows):
            if row["tic"]:
                self.by_tic[row["tic"]].append(index)
            if row["ra"] is not None and row["dec"] is not None:
                sky_xyz.append(_xyz(row["ra"], row["dec"]))
                sky_indices.append(index)
        self._sky_indices = np.asarray(sky_indices, dtype=int)
        self._tree = cKDTree(sky_xyz) if sky_xyz else None
        self.known_tics = set(self.by_tic)
        self.eb_tics = {r["tic"] for r in self.rows if r["catalog"] == "tesseb" and r["tic"]}
        self.confirmed_tics = {r["tic"] for r in self.rows if r["catalog"] == "confirmed" and r["tic"]}

    def _source(self, name, path, required=True):
        rows, metadata = _read_verified(path, required)
        if metadata is None:
            self.missing_optional.append(name)
        else:
            self.provenance[name] = metadata
        return rows

    def _add(self, catalog, name, tic, ra, dec, period, epoch, duration,
             transiting=False, rv=False, disposition="", reference="", **extra):
        self.rows.append({"catalog": catalog, "name": str(name), "tic": tic_number(tic),
                          "ra": numeric(ra), "dec": numeric(dec), "period_days": numeric(period),
                          "epoch_btjd": numeric(epoch), "duration_days": numeric(duration),
                          "transiting": bool(transiting), "rv": bool(rv),
                          "disposition": str(disposition), "reference": str(reference), **extra})

    def _load(self):
        p = self.root / "data" / "catalogs"
        b = self.root / "data" / "broad_catalogs"
        for r in self._source("confirmed", p / "confirmed_all.csv"):
            epoch = numeric(r.get("pl_tranmid"))
            duration = numeric(r.get("pl_trandur"))
            self._add("confirmed", r.get("pl_name"), r.get("tic_id"), r.get("ra"), r.get("dec"),
                      r.get("pl_orbper"), epoch - 2457000 if epoch is not None else None,
                      duration / 24 if duration is not None else None,
                      transiting=r.get("tran_flag") == "1", rv=r.get("rv_flag") == "1",
                      disposition=r.get("discoverymethod", ""), reference=r.get("pl_orbper_reflink", ""))
        for r in self._source("toi", p / "toi.csv"):
            epoch = numeric(r.get("pl_tranmid"))
            duration = numeric(r.get("pl_trandurh"))
            self._add("toi", r.get("toidisplay"), r.get("tid"), r.get("ra"), r.get("dec"),
                      r.get("pl_orbper"), epoch - 2457000 if epoch is not None else None,
                      duration / 24 if duration is not None else None,
                      transiting=True, disposition=r.get("tfopwg_disp", ""))
        for r in self._source("ctoi", p / "ctoi.csv"):
            epoch = numeric(r.get("Transit Epoch (BJD)"))
            duration = numeric(r.get("Duration (hrs)"))
            self._add("ctoi", "CTOI " + r.get("CTOI", ""), r.get("TIC ID"), r.get("RA"), r.get("Dec"),
                      r.get("Period (days)"), epoch - 2457000 if epoch is not None else None,
                      duration / 24 if duration is not None else None,
                      transiting=True, disposition=r.get("User Disposition", ""), reference=r.get("Paper", ""))
        for r in self._source("diamante_v2", b / "diamante_v2.csv", False):
            period = numeric(r.get("periodLin"))
            epoch = numeric(r.get("t0Lin"))
            if period is None or period <= 0 or epoch is None:
                period, epoch = numeric(r.get("periodBLS")), numeric(r.get("t0Fit"))
            self._add("diamante", "DIAmante v2 TIC " + r.get("ticID", ""), r.get("ticID"),
                      r.get("ra"), r.get("dec"), period, epoch, r.get("duration"),
                      transiting=True, disposition="published candidate; binflag=" + r.get("binflag", ""),
                      reference="2023MNRAS.518L..31M; DOI 10.17909/t9-p7k6-4b32")
        for r in self._source("tesseb", b / "tesseb_members.csv", False):
            # Membership is explicit, never inferred from an object detail link.
            if str(r.get("in_catalog", "")).strip().lower() not in ("true", "1", "yes"):
                continue
            self._add("tesseb", "TESS EB " + r.get("tic", ""), r.get("tic"), r.get("ra"),
                      r.get("dec"), r.get("period_days"), r.get("epoch_btjd"), None,
                      disposition="explicit catalog member", transiting=False,
                      reference=r.get("reference", "https://tessebs.villanova.edu/"))
        for r in self._source("published_candidates", b / "published_candidates.csv", False):
            self._add("published", r.get("name"), r.get("tic"), r.get("ra"), r.get("dec"),
                      r.get("period_days"), r.get("epoch_btjd"), r.get("duration_days"),
                      transiting=True, disposition="previously published candidate",
                      reference=r.get("reference", ""))

    def _matches(self, tic, ra=None, dec=None):
        tic = tic_number(tic)
        indices = set(self.by_tic.get(tic, []))
        ra, dec = numeric(ra), numeric(dec)
        separations = {}
        if ra is not None and dec is not None and self._tree is not None:
            point = _xyz(ra, dec)
            radius = 2 * np.sin(np.deg2rad(self.position_radius_arcsec / 3600) / 2)
            for local_index in self._tree.query_ball_point(point, radius):
                index = int(self._sky_indices[local_index])
                row = self.rows[index]
                chord = np.linalg.norm(point - _xyz(row["ra"], row["dec"]))
                separations[index] = float(np.rad2deg(2 * np.arcsin(min(chord / 2, 1))) * 3600)
                indices.add(index)
        matches = []
        for index in sorted(indices):
            row = self.rows[index]
            identifier_match = bool(tic and row["tic"] == tic)
            distance = separations.get(index)
            host_match = identifier_match or (distance is not None and distance <= self.host_radius_arcsec)
            matches.append({**row, "identifier_match": identifier_match,
                            "separation_arcsec": distance, "host_match": host_match})
        return matches

    def match_target(self, tic, ra=None, dec=None):
        matches = self._matches(tic, ra, dec)
        host_catalogs = {r["catalog"] for r in matches if r["host_match"]}
        return {"tic": tic_number(tic), "is_confirmed_host": "confirmed" in host_catalogs,
                "is_toi_host": "toi" in host_catalogs, "is_ctoi_host": "ctoi" in host_catalogs,
                "is_eb_host": "tesseb" in host_catalogs,
                "is_diamante_host": "diamante" in host_catalogs,
                "is_published_host": "published" in host_catalogs,
                "nearby_catalog_object": any(not r["host_match"] for r in matches),
                "matches": matches,
                "novelty_status": "unchecked_literature_and_all_aliases",
                "missing_optional_catalogs": list(self.missing_optional)}

    def match_signal(self, tic, period_days, ra=None, dec=None, tolerance=0.01, max_harmonic=5):
        result = self.match_target(tic, ra, dec)
        period = numeric(period_days)
        known_matches = []
        if period is not None and period > 0:
            factors = [1.0] + list(range(2, max_harmonic + 1)) + [1 / n for n in range(2, max_harmonic + 1)]
            for row in result["matches"]:
                known = row["period_days"]
                if known is None or known <= 0:
                    continue
                factor = min(factors, key=lambda f: abs(period / (known * f) - 1))
                mismatch = abs(period / (known * factor) - 1)
                if mismatch <= tolerance:
                    known_matches.append({**row, "harmonic_factor": float(factor),
                                          "fractional_period_mismatch": float(mismatch)})
        host = [r for r in known_matches if r["host_match"]]
        result.update({"period_days": period, "known_matches": known_matches,
                       "known_period_match": any(r["harmonic_factor"] == 1 for r in host),
                       "known_harmonic_match": any(r["harmonic_factor"] != 1 for r in host),
                       "known_rv_period_match": any(r["rv"] for r in host),
                       "possible_nearby_period_contamination": any(not r["host_match"] for r in known_matches),
                       "period_match_tolerance_fraction": tolerance})
        result["period_status"] = ("known_period_match" if result["known_period_match"] else
                                   "known_harmonic_match" if result["known_harmonic_match"] else
                                   "no_cached_period_match")
        return result

    def transit_mask(self, tic, t_btjd, ra=None, dec=None, width_factor=1.7):
        t = np.asarray(t_btjd, dtype=float)
        keep = np.ones(t.shape, dtype=bool)
        for row in self._matches(tic, ra, dec):
            p, t0, duration = row["period_days"], row["epoch_btjd"], row["duration_days"]
            if not row["host_match"] or not row["transiting"] or p is None or p <= 0 or t0 is None:
                continue
            # Missing duration alone can use an explicit conservative 3-hour mask.
            duration = duration if duration is not None and duration > 0 else 3 / 24
            halfwidth = max(width_factor * duration / 2, .05)
            # Avoid nonsensical ephemerides masking most of the time series.
            if 2 * halfwidth >= .3 * p:
                continue
            phase = (t - t0 + p / 2) % p - p / 2
            keep &= np.abs(phase) > halfwidth
        return keep

    def annotate_dataframe(self, dataframe):
        """Return a copy with scalar crossmatch flags and compact match evidence."""
        out = dataframe.copy()
        records = []
        for r in out.to_dict("records"):
            check = self.match_signal(r.get("tic"), r.get("period_days", r.get("period")),
                                      r.get("ra"), r.get("dec"))
            records.append({k: v for k, v in check.items()
                            if k not in ("tic", "period_days", "matches", "known_matches", "missing_optional_catalogs")})
            records[-1]["catalog_matches_json"] = json.dumps(check["matches"], separators=(",", ":"))
            records[-1]["known_period_matches_json"] = json.dumps(check["known_matches"], separators=(",", ":"))
        import pandas as pd
        extra = pd.DataFrame(records, index=out.index)
        for column in extra:
            out[column] = extra[column]
        return out

    def metadata(self):
        return {"checked_utc": datetime.now(timezone.utc).isoformat(),
                "position_radius_arcsec": self.position_radius_arcsec,
                "host_radius_arcsec": self.host_radius_arcsec,
                "row_count": len(self.rows), "tic_count": len(self.known_tics),
                "catalogs": self.provenance, "missing_optional_catalogs": self.missing_optional,
                "limitations": "Cached catalog matching only. Nonmatch is not novelty. Sky positions are not propagated for proper motion; a 30 arcsec match may identify a contaminant, not the target. Harmonics are ranking flags, not signal identifications. Per-target literature, aliases, TCE/DV and pixel checks remain required."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, help="CSV with tic, optional ra/dec and period or period_days")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--tic", type=int)
    parser.add_argument("--period", type=float)
    parser.add_argument("--retrieve-support", action="store_true", help="Restore missing supportive catalog files and normalized tables")
    args = parser.parse_args()
    if args.retrieve_support:
        retrieve_support()
    index = CatalogIndex()
    if args.input:
        if not args.output or args.output.resolve() == args.input.resolve():
            parser.error("A distinct --output is required; source tables are never overwritten")
        import pandas as pd
        result = index.annotate_dataframe(pd.read_csv(args.input))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, index=False)
        print(f"Annotated {len(result)} rows -> {args.output}")
    elif args.tic:
        result = index.match_signal(args.tic, args.period) if args.period else index.match_target(args.tic)
        print(json.dumps(result, indent=2))
    else:
        print(json.dumps(index.metadata(), indent=2))
    if args.metadata:
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(index.metadata(), indent=2) + "\n")


if __name__ == "__main__":
    main()
