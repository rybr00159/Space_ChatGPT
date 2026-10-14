# TIC49617264: independent native-pixel source test

**The 4.83194-day signal localizes to faint TIC49617264, rather than the
nearby bright TIC49617263. It is an already-reported, ground-followed neighbor
planetary candidate, TOI2417.11/TIC49617264.11. This investigation did not
discover it and does not confirm its planetary nature.** Checked 2026-10-14 UTC.

The original bright target TIC49617263 carries TOI2417.01 with a nearly
identical listed period, but is marked false positive. A coordinate-only or
period-only catalog rejection would be misleading: existing follow-up notes
explicitly relocate the transit to TIC49617264. The preserved public ExoFOP
2024-10-05 report describes a 13.1-ppt on-target transit in an uncontaminated
6.6-arcsec Rc aperture and separately cites a 2024-09-15 TRAPPIST observation
with a 12.8-ppt event in a 4.5-arcsec aperture. These are prior observer reports,
not new measurements made here. See
`data/broad_catalogs/tic49617264_followup20241005_notes.txt`, ExoFOP fid1004492,
and their retrieval manifest for provenance.

## Pixel result

Two native 2-minute SPOC target-pixel products, from sectors 106 and 107,
provide different camera orientations and 9 complete local event windows.
All nonzero QUALITY flags and nonfinite flux/error/pointing/background rows
were removed. The target-pixel stamp contains both stars and the bright
neighbor's close Gaia companion. Pixels are already background-subtracted;
the pipeline FLUX_BKG is used as a temporal nuisance regressor, not subtracted
a second time.

| Sector | Camera / CCD | Usable samples | Event windows | Difference-source offset from target | Distance from bright neighbor | Target r=1 pixel depth |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 106 | 3 / 2 | 12,739 | 5 | 0.015 pixels | 2.059 pixels | 7.14 ± 0.28 ppt |
| 107 | 2 / 1 | 12,884 | 4 | 0.061 pixels | 2.067 pixels | 5.38 ± 0.34 ppt |

The independent raw SPOC-optimal pixel aperture yields fixed-trapezoid depths
6.65 ± 0.24 ppt (sector 106) and 8.91 ± 0.49 ppt (sector 107). These are
**diluted observed-aperture depths**, not intrinsic stellar eclipse depths or
planet-radius measurements. Increasing the target-centered aperture radius
adds the much brighter neighbor's light and decreases the observed dip: at
r=2 pixels the depths are 1.62 ± 0.14 ppt and 2.56 ± 0.24 ppt. The brightest
pixel's core aperture is consistent with no dip, while its broader apertures
include some light from the transiting faint target. A nonzero broader
neighbor-aperture dip therefore does not establish that neighbor as the source.

The positive transit difference lies at TIC49617264 in both camera
orientations. The Gaia DR3 target position is 0.004 pixels from the sector-106
difference centroid and 0.066 pixels from the sector-107 centroid. The bright
neighbor's Gaia DR3 source is 4952839741811566720; its close companion
4952839741810218880 is also ~2 pixels from the difference source. Their similar
positions are unresolved in these pixels, but both are separated from the
faint transit source.

## Method and limitations

Each event uses a +/-0.40-day local simultaneous fit of a quadratic time
baseline, temporal background, POS_CORR1/2 pointing terms, and a
cadence-integrated trapezoid. Its fixed period is the public TCE period
4.83194208145142 days. Initial T0=4265.644028635728 BTJD follows the native
sector-107 fine-BLS timing correction (~6 minutes earlier than the public
epoch); duration is 2.286 hours with 20% ingress. This fixed template describes
the source test and is not a physical planet model or updated global ephemeris.
TIME is TDB with BJD reference 2457000; its barycentric correction is already
applied.

Alternative quadratic-only, broad 2.77-hour box, and 1.56-hour core-box fits
all locate the source near the faint target. Their catalog-position offsets
span 0.011–0.027 pixels in sector 106 and 0.036–0.103 pixels in sector 107.
An elliptical pixel-integrated Gaussian plus affine image background fits
the combined difference image with diagonal chi-square/dof of 106.8/113 and
110.9/113. Fixing that model at the bright neighbor instead substantially
worsens its fit (delta chi-square ~749 and ~470). Those values are morphology
diagnostics, not calibrated source probabilities or search false-alarm rates.

Five hundred bootstrap draws per sector resample complete events plus
synchronized residual images in 10-cadence (~20-minute) blocks. Conditional
16th–84th-percentile target offset ranges are 0.022–0.077 pixels and
0.042–0.128 pixels. All draws are closer to the faint target than to either
bright-neighbor Gaia component. Gaussian/PRF systematics, WCS uncertainty,
proper-motion epoch differences, and unresolved blends are **not** covered
by these intervals. No official TESS PRF scene model was fitted. Source
localization alone cannot exclude an unresolved eclipsing companion or
establish a planetary mass.

The TPF extractions share photons with the SPOC light curves; they do not
provide independent statistical confirmation. Their value is the spatial
source check and aperture-dependence test. The archival ground photometry
is separate prior evidence and independently establishes that this source
was already investigated years before this project.

## Evidence and reproduction

`sector106/` and `sector107/` preserve JSON measurements, per-event aperture
CSVs, Gaia-to-pixel positions, alternative difference images with errors,
bootstrap coordinates, extracted light curves, native optimal-aperture masks,
and diagnostic plots. Original raw TPFs remain at
`/workspace/tess-research-data/tce107/pixels/`, totaling about 101 MB decimal
(97 MiB). They are excluded from ordinary Git and never deleted by the script.

`data/tce107_pixel_manifest.json` records MAST observation/product identifiers,
download URLs, native processing versions, camera/CCD, cadence and time
system, original SHA-256 hashes, and passing raw FITS one's-complement
checksum verification. The raw verifier includes full HDU data even when a
DATASUM card is absent; no verification was disabled. `products_s106.json`
and `products_s107.json` preserve original MAST product-list responses.
Cache persistence in a future environment is not guaranteed, but exact
public product URIs allow re-retrieval; retain divergent downloads for
diagnosis instead of replacing expected hashes silently.

From the repository root, use its installed Python environment:

```bash
OPENBLAS_NUM_THREADS=1 python src/pixel_vet_49617264.py --sectors 106 107 --bootstrap 500 --gaia data/broad_catalogs/tic49617264_gaia_neighbors_current.csv
```

The script caches missing native TPFs and verifies their original FITS
integrity before analysis. It does not spend paid API funds. Prior notes,
Gaia queries, and recognition checks were retrieved by the catalog work and
are preserved in `data/broad_catalogs/` with their own provenance.

**Disposition: known prior candidate; exclude from a new-discovery list.**
These tests provide a useful positive control for source localization in a
strongly contaminated TESS aperture. Recommended scientific follow-up would
belong to the existing candidate's validation effort, including high-resolution
imaging and mass constraints, rather than a novelty claim from this search.
