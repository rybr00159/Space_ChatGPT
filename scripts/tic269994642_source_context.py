"""Derive TIC 269994642 neighborhood evidence from retained public responses.

No catalog downloads or original response files are modified. Request URLs,
ADQL, timestamps and hashes in the adjacent manifests enable retrieval replay.
"""
import csv
import hashlib
import io
import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "broad_catalogs"
PREFIX = "tic269994642_"


def separation(ra, dec, target_ra, target_dec):
    a, d, ta, td = map(math.radians, (ra, dec, target_ra, target_dec))
    v = math.sin((d - td) / 2) ** 2
    v += math.cos(d) * math.cos(td) * math.sin((a - ta) / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, v)))) * 3600


def main():
    gaia_path = BASE / (PREFIX + "gaia_nearby.csv")
    tic_path = BASE / (PREFIX + "tic_cone120.json")
    rows = list(csv.DictReader(gaia_path.open()))
    target = next(r for r in rows if r["source_id"] == "6806240691223695744")
    tr, td, tg = (float(target[k]) for k in ("ra", "dec", "phot_g_mean_mag"))
    tic_rows = json.loads(tic_path.read_text())["data"]
    # TIC GAIA is DR2; common IDs map directly. Do not force changed DR3 IDs.
    tic_by_gaia = {str(r["GAIA"]): r for r in tic_rows if r.get("GAIA")}
    for row in rows:
        row["separation_arcsec_from_gaia_target"] = separation(
            float(row["ra"]), float(row["dec"]), tr, td)
        row["full_g_flux_ratio_to_target"] = (
            10 ** (-0.4 * (float(row["phot_g_mean_mag"]) - tg))
            if row["phot_g_mean_mag"] else "")
        tic = tic_by_gaia.get(row["source_id"], {})
        row["tic_dr2_id_match"] = tic.get("ID", "")
        row["tic_tmag"] = tic.get("Tmag", "")
        row["full_tess_flux_ratio_to_target"] = (
            10 ** (-0.4 * (float(tic["Tmag"]) - 12.9234))
            if tic.get("Tmag") is not None else "")
        row["is_gaia_target"] = row["source_id"] == target["source_id"]
    rows.sort(key=lambda r: r["separation_arcsec_from_gaia_target"])
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    path = BASE / (PREFIX + "gaia_neighbors_derived.csv")
    content = output.getvalue().encode()
    if path.exists() and path.read_bytes() != content:
        raise RuntimeError("Preserve previous derived evidence before replacing it")
    path.write_bytes(content)
    manifest = {
        "derived_utc": datetime.now(timezone.utc).isoformat(),
        "sha256": hashlib.sha256(content).hexdigest(),
        "bytes": len(content),
        "row_count": len(rows),
        "method_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "120 arcsec Gaia cone; separations use Gaia target coordinates. "
                 "Full catalog flux ratios are not aperture fractions. TIC GAIA "
                 "is DR2, so changed source IDs remain unmapped here.",
        "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in (gaia_path, tic_path)},
    }
    manifest_path = Path(str(path) + ".manifest.json")
    if not manifest_path.exists():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    inventory = []
    for p in sorted(BASE.glob(PREFIX + "*.manifest.json")):
        source = Path(str(p)[:-len(".manifest.json")])
        metadata = json.loads(p.read_text())
        actual = hashlib.sha256(source.read_bytes()).hexdigest()
        if actual != metadata["sha256"]:
            raise RuntimeError("Source hash mismatch: " + str(source))
        inventory.append({"file": str(source.relative_to(ROOT)),
                          "manifest": str(p.relative_to(ROOT)),
                          "bytes": source.stat().st_size, "sha256": actual})
    inventory_path = ROOT / "reports" / "tic269994642_prior_check_sources.json"
    inventory_path.write_text(json.dumps({
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "verified_source_count": len(inventory), "sources": inventory}, indent=2) + "\n")
    print("Verified", len(inventory), "source snapshots; derived", len(rows), "Gaia rows")


if __name__ == "__main__":
    main()
