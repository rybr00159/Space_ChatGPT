"""Restore and independently vet two further public TCEs after the pilot.

Exact existing raw-product hashes are immutable. The two-star sample and eight
sectors were selected from metadata before flux inspection. No paid API.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import urllib.request

import pandas as pd
from astropy.io import fits

from archive import verify_fits_checksums
from broad_crossmatch import CatalogIndex
from vet_tce107 import ROOT, OUT, safe, vet

TICS = [49617264, 89187590]
SELECTED_PUBLIC_TCES = {49617264: "00049617264-01", 89187590: "00089187590-01"}
MANIFEST = ROOT / "data/tce107_extra_observations.json"


def verify_restore(records, retrieve):
    cache = Path(os.environ.get("TESS_DATA_DIR", "/workspace/tess-research-data"))
    if sum(r["bytes"] for r in records) > 20_000_000 or len(records) > 12:
        raise ValueError("This frozen additional-light-curve experiment exceeds its 20 MB / 12-product scope")
    for record in records:
        path = cache / "lightcurves" / Path(record["local_path"]).name
        if not path.exists():
            if not retrieve:
                raise FileNotFoundError(f"{path}; --retrieve restores only exact pinned products")
            path.parent.mkdir(parents=True, exist_ok=True)
            part = path.with_suffix(path.suffix + ".part")
            with urllib.request.urlopen(record["retrieval_url"], timeout=120) as response, part.open("wb") as output:
                while block := response.read(1024 * 1024): output.write(block)
            if hashlib.sha256(part.read_bytes()).hexdigest() != record["sha256"]:
                raise ValueError(f"Retrieved bytes mismatch pinned SHA; retain {part} for diagnosis")
            part.replace(path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Raw input changed; preserve and investigate: {path}")
        with fits.open(path) as h:
            h.verify("exception")
            verify_fits_checksums(path, h)
            if int(h[0].header["TICID"]) != record["tic"] or int(h[0].header["SECTOR"]) != record["sector"]:
                raise ValueError("Target identity or sector mismatch")
        record["local_path"] = str(path)
    return records


def assessment(results):
    text = ["# Further independent follow-up of sector-107 public TCEs", "",
        "Two additional stars were selected from the frozen public-TCE metadata ranking before observing their flux: TIC 49617264 (MES 17.84, four reported events, no severe initial DV flag) and TIC 89187590 (MES 14.83, three events, geometric-model S/N 7.79). Eight native SPOC observations totaling **18,095,040 bytes** were retrieved and integrity checked. This extends time-series follow-up from three to five stars and from nine to seventeen curves. The 846 public TCE metadata rows are separate from independent light-curve analysis counts.", "",
        "**TIC 49617264 recovers a previously reported neighboring planet candidate, not a new object.** The approximately 4.832-day dips repeat in sectors 97 and 105–107, persist in raw SAP and corrected PDC, and survive quadratic local baselines and alternative gap-aware detrending. The restricted sector-107 fine search gives P≈4.832023 days, epoch BTJD≈4270.47597 and box duration≈2.29 hours; these are selected follow-up estimates, not a physical orbital posterior. A wider prior search then found original bright neighbor TIC 49617263 = TOI 2417.01 (FP) only 43.18 arcseconds away, outside the default 30-arcsecond match radius, at nearly the same period. This neighboring period alone does not identify an eclipsing binary: the decisive public 2024 follow-up notes already call the actual faint source TIC 49617264.11 / TOI 2417.11 and explicitly describe it as an NPC to the original bright candidate. The source was already reported and followed as a planetary candidate before this analysis.", "",
        "The [Brierfield 2024-10-05 report](https://exofop.ipac.caltech.edu/tess/get_file.php?id=1004492) describes an on-time, on-target approximately 13.1 ppt event in a claimed uncontaminated 6.6-arcsecond aperture, with a reported 2h46m transit duration, and cites TRAPPIST 2024-09-15 zp follow-up at approximately 12.8 ppt and revised period 4.83261 days. It quotes '[NPC to original candidate TIC 49617263 (TOI-2417.01)] SAFFIs constrain event to this star relative to Gaia DR2 stars bright enough to cause the TESS detection.' These are prior observers' claims, not new confirmation by this project; the notes also contain transcription/date/telescope inconsistencies. The raw note, URL/hash manifest and context are preserved in `data/broad_catalogs/tic49617264_followup20241005_notes.txt` and `reports/tic49617264_prior_check.md`. An empty TOI field on the actual target and the original bright target's FP disposition must not be converted into either a novelty claim or an EB classification of the relocated candidate.", "",
        "Crowding is substantial (CROWDSAP approximately 0.46–0.61), corrected depth changes between sectors, and independent pixel localization is useful as a check on the prior relocation. PDC linear-local depths are approximately 0.64%, 0.86%, 1.03%, 1.15% for sectors 97,105,106,107, while SAP depths are approximately 0.30%, 0.54%, 0.52%, 0.70%. Formal odd/even comparisons are weak in the three recent sectors; an older-sector formal difference is sensitive to one questionable event and must not be treated as a calibrated binary classification. A physical radius inferred by assuming the catalog target and its dilution remains conditional. The known candidate is a useful positive recovery and catalog-recognition lesson, without a new planetary designation or discovery claim.", "",
        "**TIC 89187590 fails the proposed approximately 8.726-day planetary interpretation.** Every sector contains a strong approximately 1.09-day broad waveform. Quadratic native local baselines reduce sector-107's selected depth from approximately 1518 ppm to approximately 5 ppm, and sector-106's from approximately 278 ppm to approximately −616 ppm. The restricted fine BLS selects a different phase about 3.28 days from the public epoch at only diagonal S/N≈8.6. This is evidence of baseline curvature/variability rather than stable local transits. Earlier-sector fixed-phase negative depths alone are not used as decisive falsification because the public uncertainty columns are unreliable over seven years; recent native baselines already undermine the hypothesis. This does not exclude other planets around the star.", "",
        "| TIC | Sector | Flux | Local linear depth ppm | Formal S/N | Event-scatter S/N | Events | Local quadratic depth ppm |", "|---|---|---|---|---|---|---|---|"]
    for result in results:
        for sector in result["sector_results"]:
            for m in sector["measurements"]:
                a, b = m["fixed_linear"], m["fixed_quadratic"]
                def f(v): return "unavailable" if v is None else f"{v:.2f}"
                d, q = a.get("depth"), b.get("depth")
                text.append(f"|{result['tic']}|{sector['sector']}|{m['column']}|{f(d * 1e6 if d is not None else None)}|{f(a.get('formal_snr'))}|{f(a.get('event_scatter_snr'))}|{a['n_events']}|{f(q * 1e6 if q is not None else None)}|")
    text += ["", "Reproduce with `python src/tce107_followup_extra.py`. If the cache is lost, use `--retrieve`; exact URLs, MAST identifiers, sectors, processing versions, SHA-256 and byte counts are frozen in `data/tce107_extra_observations.json`. The original nine-product manifest remains unchanged. Raw inputs are retained outside Git; diagnostic plots and native event-level measurements are in this directory. The wrapper checks hashes and raw FITS checksums before each analysis. These follow-up ranking statistics and simple covariance errors are not global false-alarm probabilities.", "",
        "Next: preserve the known-candidate recovery and independent pixel comparison; focus discovery effort on other unrecognized signals. Future shortlisted targets need a match radius appropriate to their pixel aperture, nearby-star period checks, aliases and ground-follow-up uploads, even when the target's own catalog row is empty. This investigation found three variability/baseline aliases, one previously published complex periodic variable, and one previously reported neighboring planet candidate among the five independently followed public-TCE stars. No credible new planet emerged from this track.", ""]
    (OUT / "EXTRA_FOLLOWUP.md").write_text("\n".join(text))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, default=MANIFEST)
    p.add_argument("--retrieve", action="store_true")
    p.add_argument("--reuse-measurements", action="store_true", help="Reuse matching preserved native measurements after verifying raw hashes")
    args = p.parse_args()
    records = verify_restore(json.loads(args.manifest.read_text()), args.retrieve)
    triage = pd.read_csv(ROOT / "results/tce107_triage.csv")
    index = CatalogIndex()
    results = []
    for tic in TICS:
        # A later recognition check can legitimately mark a frozen selected
        # signal ineligible; that must not prevent reproducing its old data.
        selected = triage[triage.tceid == SELECTED_PUBLIC_TCES[tic]]
        if len(selected) != 1:
            raise ValueError("Pinned public TCE identifier missing or duplicated")
        lead = selected.iloc[0]
        rows = [r for r in records if r["tic"] == tic]
        saved = OUT / f"tic{tic}_measurements.json"
        if args.reuse_measurements and saved.exists():
            result = json.loads(saved.read_text())
            if {r["data_uri"]: r["sha256"] for r in result["observations"]} != {r["data_uri"]: r["sha256"] for r in rows}:
                raise ValueError("Existing measurement source set changed; do not reuse")
            if abs(result["ephemeris"]["period_days"] - float(lead.period_days_followup)) > 1e-12 or abs(result["ephemeris"]["epoch_btjd"] - float(lead.epoch_btjd_followup)) > 1e-10:
                raise ValueError("Existing measurement selection ephemeris changed; do not reuse")
        else:
            result = vet(tic, rows, lead, index)
        result["additional_experiment"] = {"sample_selected_from_metadata_before_flux": True,
            "manifest": str(args.manifest.relative_to(ROOT)), "previous_baseline": "three stars / nine curves"}
        result["scientific_claim"] = ("recovery of previously reported neighboring planet candidate TIC 49617264.11 / TOI 2417.11; no novel claim" if tic == 49617264 else
                                     "public 8.726-day interpretation undermined by native local curvature and short-period variability")
        if tic == 49617264:
            notes = ROOT / "data/broad_catalogs/tic49617264_followup20241005_notes.txt"
            provenance = json.loads(notes.with_suffix(notes.suffix + ".manifest.json").read_text())
            if hashlib.sha256(notes.read_bytes()).hexdigest() != provenance["sha256"]:
                raise ValueError("Positive prior-recognition source changed; preserve and investigate")
            result["positive_prior_recognition"] = {"status": "previously reported neighboring planet candidate",
                "name": "TIC 49617264.11 / TOI 2417.11", "original_bright_host": "TIC 49617263 / TOI 2417.01 (FP)",
                "neighbor_separation_arcsec": 43.18, "source_date": "2024-10-05",
                "source_manifest": provenance, "ground_claimed_depth_ppt": 13.1,
                "ground_claimed_uncontaminated_aperture_arcsec": 6.6,
                "independent_confirmation_by_this_project": False,
                "scope": "authoritative public upload establishes prior recognition, not planetary confirmation"}
        (OUT / f"tic{tic}_measurements.json").write_text(json.dumps(safe(result), indent=2, allow_nan=False))
        results.append(result)
    assessment(results)


if __name__ == "__main__":
    main()
