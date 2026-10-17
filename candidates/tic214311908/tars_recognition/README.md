# TIC 214311908: TARS recognition audit

The presence of this target in the TESS All-Sky Rotation Survey (TARS) is **not a prior planet-candidate designation**. The MAST HLSP includes light curves for approximately 39 million star-sector combinations of 7.5 million targets, selected by magnitude and parallax. Its principal published product is a stellar rotation/variability-period catalog. It does not classify individual transit candidates.

This audit was performed independently, without the other model's findings. Public source material and exact retrieval/checksum records are in `data/tic214311908/tars_recognition/`. Original FITS light curves and preview images are retained outside ordinary Git under `/workspace/tess-research-data/tic214311908/independent_reductions/`.

## Reported per-sector periods

The following are the **published TARS header values**, not fitted transit periods in this project:

| Sector | Highest Lomb–Scargle period (d) | Other relevant listed period (d) | Its rank |
|---|---:|---:|---:|
| 12 | 8.5715718 | 5.7368405 | 3 |
| 39 | 8.7464183 | none near 5.69 among the five listed peaks | — |
| 66 | 8.0918201 | 5.7542783 | 5 |
| 93 | 4.5891302 | 5.8904204 | 3 |

The near-5.7-day values are broad, variable lower-ranked sinusoidal-periodogram peaks. The TARS headers do not report a coherent approximately 5.6915-day transit ephemeris, transit epoch, or planet classification. This finite check does not prove that a transit signal is previously unpublished: a generic LS peak can contain power from transits and is itself a prior public measurement.

The headers' `SNR1`–`SNR5` mean peak power divided by median periodogram power, according to the TARS README. These are **not Gaussian transit-significance measurements**. The preview classifier scores are rotation/variability and alias diagnostics, not exoplanet probabilities.

## Independent reduction availability and limitations

Four TARS original FITS files (Sectors 12, 39, 66, 93) and one GSFC-ELEANOR-LITE original FITS file (Sector 12) were retrieved and parsed. They use the same TESS photons as the QLP observations; agreement across reductions tests processing sensitivity but does not add independent observing seasons.

TARS light curves contain `TIME` and normalized `FLUX`, with `APERTURE=1` in these four files. They do not provide per-cadence uncertainties or quality flags. Removing data or estimating uncertainties therefore requires an explicit policy; these values should not silently inherit QLP uncertainties. The TARS previews display a TIC contamination ratio of **1.8755**, and a crowded DSS image. The target's raw aperture flux cannot be assumed to consist of target light alone.

GSFC-ELEANOR-LITE Sector 12 provides `TIME`, `RAW_FLUX`, `CORR_FLUX`, `PCA_FLUX`, `FLUX_ERR`, `QUALITY`, background, and centroid columns (1,289 rows). Its aperture is `L_180`; HDU 2 contains 15×15 aperture/reference-image arrays, not a time-resolved target-pixel cube. Full header and schema snapshots accompany each retrieval manifest.

## Published default catalog

The complete published default rotation catalog, Zenodo release **19917941**, contains **1,046,317 rows**. A full row scan finds **no TIC 214311908 entry**. The audit also checks exact Gaia source ID **5925373272850784768** and the rounded ID printed in TARS previews. The original 136,182,994-byte CSV ZIP is retained outside Git; its MD5 matches the publisher's checksum (`49e0b7538840dea83ecbde3a48387fd6`).

The complete scan result, version, checksums, row count, identifier policy, and source URL are stored in `data/tic214311908/tars_recognition/default_catalog_check.json`. The script is `scripts/tic214311908_tars_recognition.py`. The 13.37 GB full sector-measurement table is not required for the four archived per-sector headers and is not downloaded by this script.

**Identifier precision caution:** TARS stores `GAIADR2`/`GAIADR3` as floating-point values in these FITS headers (`5.92537327285078E+18`). This loses integer precision. The previews display `5925373272850784256`; the actual Gaia archive returns source ID `5925373272850784768` at the same coordinates. The exact integer from the original Gaia CSV governs identifier matching. Neither a FITS floating-point conversion nor the rounded preview should replace it.

Thus TARS establishes that this star has pre-existing public light curves and generic sinusoidal period measurements. It does **not** supply a default adopted rotation-period entry or a planet-candidate designation for this TIC. This conclusion is limited to these specific public products and this pinned catalog release.

## Sources

- [MAST TARS overview and data products](https://archive.stsci.edu/hlsp/tars), DOI [10.17909/2sff-fn29](https://doi.org/10.17909/2sff-fn29).
- [Official TARS README](https://archive.stsci.edu/hlsps/tars/hlsp_tars_readme.txt).
- Boyle, Bouma & Mann, *The TESS All-Sky Rotation Survey*, [arXiv:2603.05586](https://arxiv.org/abs/2603.05586). The downloaded HTML is explicitly version 1; its default-catalog size differs from the subsequently released catalog, so the actual pinned Zenodo files govern the row-level audit.
- [Current pinned Zenodo release 19917941](https://zenodo.org/records/19917941), reached via concept record 18342590 on 2026-10-17 UTC. Download URLs and server-published MD5 are recorded and checked.

No contact with the authors or external submission was made. No paid API was used.
