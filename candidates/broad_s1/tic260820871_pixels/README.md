# TIC260820871: source falsification with sector-28 pixels

**Both alternating eclipse sequences localize near the already-cataloged
faint eclipsing source Gaia DR3 6381543192438078336, rather than the bright,
hot TIC260820871. Blended eclipses strongly explain the shallow dips in the
bright target's aperture. Exclude this signal from a new planet list on
TIC260820871.** Checked 2026-10-14 UTC.

This conclusion comes from the spatial difference images and aperture
dependence. The neighbor's VSX catalog period, 2.692080 days, is only
approximately twice the tested 1.346665-day alternating-eclipse period.
The discrepancy is real; neither an exact period identity nor a revised
catalog orbital period is claimed here. Pixel localization supports the
known neighbor's vicinity independently of that approximate harmonic.

## Independent spatial result

A retained 15x15 TESSCut stamp from sector 28 contains both stars separated
by approximately 1.67 pixels (~33.7 arcsec). It uses original 10-minute SPOC
FFIs, with 3,449 usable QUALITY=0 frames of 3,464 total. The bright target is
on science silicon. Deep and shallow eclipse sequences were analyzed
separately at phase 0 and 0.5 of the doubled photometric period.

| Eclipse sequence | Local windows | Difference centroid x,y | Offset from known EB | Offset from hot TIC target |
| --- | ---: | --- | ---: | ---: |
| Phase 0, deeper | 18 | 6.077, 7.566 | 0.140 pixels | 1.659 pixels |
| Phase 0.5, shallower | 17 | 6.086, 7.384 | 0.064 pixels | 1.629 pixels |

The neighbor's catalog WCS position is (6.042,7.431), and the bright TIC
position is (7.712,7.284). The pooled difference center is (6.070,7.502),
about 0.077 pixels from the known EB position. Both sequences have their
strongest positive difference pixels near column6, while the bright target's
core lies near columns7–8.

Quadratic-only, background/pointing, broad-box, and core-box extractions
agree spatially. Across these choices, deeper-sequence positions are
0.125–0.156 pixels from the neighbor and 1.602–1.659 pixels from the host.
Shallower-sequence positions are 0.064–0.108 pixels from the neighbor and
1.600–1.728 pixels from the host.

Five hundred whole-event plus synchronized residual-block bootstrap draws
per sequence all lie closer to the catalog EB than to the bright TIC star.
Conditional 16th–84th-percentile EB-offset ranges are 0.105–0.199 pixels for
the deeper sequence and 0.070–0.307 pixels for the shallower sequence.
These are not source-identification probabilities or astrometric exclusion
limits: Gaussian/PRF and WCS systematics, proper motion, and close unresolved
sources are not represented.

Fixing the empirical Gaussian source at the hot target instead of allowing
a free position worsens diagonal chi-square by approximately 1,098 and 315.
Fixing it at the EB position worsens chi-square by approximately 20.2 and
0.9. These are descriptive morphology comparisons, not calibrated likelihood
ratios or search-wide false-alarm statistics. The small remaining neighbor
offset can reflect the empirical image model; no official TESS PRF scene
fit was performed.

## Aperture falsification

The smallest bright-target core aperture is consistent with little or no
dip, whereas the smallest neighbor-centered aperture detects both eclipse
sequences strongly:

| Aperture, r=1 pixel | Deeper sequence | Shallower sequence |
| --- | ---: | ---: |
| Bright TIC core | 103 ± 107 ppm | 86 ± 61 ppm |
| Known EB neighborhood | 3,946 ± 132 ppm | 2,189 ± 135 ppm |

Larger bright-target apertures incorporate light from the eclipsing neighbor
and yield shallow diluted dips: r=2-pixel target apertures measure
445 ± 58 ppm and 279 ± 48 ppm. These independent extractions reproduce the
mechanism needed to contaminate the original bright-target light curves.
The neighbor aperture also contains substantial bright-target wings, so its
measured depths are still diluted and do not describe the intrinsic EB
eclipse depth. No stellar or planet radius follows from these aperture depths.

## Method and limitations

The fixed trial period is 1.3466651951673554 days, with epoch
1325.7398150267488 BTJD and duration 2.131 hours. These values come from
the independent multi-sector QLP timing/alternating-depth investigations;
the phase is not selected to favor either sky position. The integrated
trapezoid uses ingress fraction0.408. Phase 0.5 uses the same duration and
a freely fitted separate depth.

Each +/-0.25-day event window simultaneously fits the eclipse amplitude,
quadratic time baseline, measured blank-sky background, and image-pointing
terms. Blank sky is the lowest25% of the median reference-image pixels.
Registration fits reference-image amplitudes and gradients outside the
bright target's 2.5-pixel region and the neighbor's 2-pixel region, protecting
both source cores from influencing the pointing solution. The quadratic-only
comparison does not use those measured pointing terms and reaches the same
source conclusion.

Pixel difference images combine individual events with inverse pixel-error
variances. A pixel-integrated elliptical Gaussian plus affine background
provides an empirical centroid. Bootstrap residual images are resampled
synchronously in blocks of10 usable cadences, nominally100 minutes, alongside
resampling complete events. Spatial covariance and some local temporal noise
are thus retained. Conditional aperture errors take the larger of propagated
local residual scatter and scatter among event depths divided by sqrt(N).
These approximations do not establish confirmation or a discovery.

The evidence excludes an eclipse centered on the resolved hot star as a
good interpretation of these shallow events. A hypothetical unresolved
source near the known faint EB cannot be eliminated by this simple image
model. Exact orbital-period reconciliation and detailed source validation
would require further catalog/epoch-photometry work; neither is necessary
to reject a new planet claim on the bright TIC star.

## Preservation and reproduction

The original sector-28 cutout is retained outside ordinary Git tracking at
`/workspace/tess-research-data/broad/pixels/tic260820871/`, about15.9 MB.
`data/broad_pixel_260820871.json` records the original byte and scientific-array
SHA-256 hashes, passing full raw FITS HDU checksum verification, contributing
FFI identifiers, archive request, processing version, data release41,
WCS, camera3/CCD2, detector position, and time metadata.

TIME is already barycentric TDB, with BJD reference2457000 and
TIMEREF=SOLARSYSTEM; TIMECORR is not added again. DATE-OBS is
2020-07-31T08:23:31.355 UTC. The first raw midpoint converts to
08:28:31.346 UTC, +299.991 seconds, half a 10-minute cadence. That first
frame has QUALITY32 and is excluded; the first usable midpoint is
08:38:31.346 UTC. The timestamp/cadence consistency check passes.

`phase0/` and `phase_half/` preserve all baseline/box variants, per-event
difference images and errors, bootstrap locations, individual aperture
measurements, extracted local light curves, and diagnostic images.
`pixel_nuisance_arrays.npz` retains usable timestamps, reference and sky
masks, background and pointing vectors. The combined difference image,
Gaia pixel positions, original retrieval information, and complete FFI
product identifiers are also retained.

The Gaia/VSX recognition responses and their provenance are in
`data/broad_catalogs/tic260820871_*current*`. Existing code and data for other
targets were not changed. No paid API requests were used.

Run from the repository root with its installed Python environment:

```bash
OPENBLAS_NUM_THREADS=1 python src/pixel_vet_260820871.py --bootstrap 500
```

The script caches a bounded100-MB public archive response when the raw cutout
is absent. Cache persistence across future tasks is unverified. Regenerated
Astrocut checksum timestamps may change byte hashes; original expected hashes
remain preserved, and any regenerated hash is recorded only when scientific
arrays and WCS match exactly and FITS integrity verification passes. A science
or WCS mismatch stops the analysis and retains the downloaded file for review.
