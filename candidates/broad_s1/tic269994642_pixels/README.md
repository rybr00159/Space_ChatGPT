# TIC269994642: independent FFI recurrence and source checks

**A ~1.7464-day eclipse-like signal recurs in independently extracted TESS
pixels from sectors 1 and 92, separated by almost seven years. Both observing
geometries favor the target's neighborhood over resolved catalog neighbors.
Planetary origin and novelty remain unestablished.** Checked 2026-10-14 UTC.

The resolved-source checks support further investigation. The native
200-second sector-92 dip is strongly V-shaped, and Gaia reports target
RUWE 2.08 with significant astrometric excess noise. These are material
reasons to investigate grazing stellar eclipses and unresolved multiplicity.
The source test cannot resolve those explanations.

| Sector | Camera / CCD | Cadence | Usable frames | Local event windows | Target r=2-pixel trapezoid depth | Approximate aperture SNR |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| 1 | 1 / 4 | 30 minutes | 1,066 | 12 | 0.02719 ± 0.00146 | 18.6 |
| 92 | 2 / 1 | 200 seconds | 9,953 | 13 | 0.02751 ± 0.00108 | 25.5 |

These fixed-template depths are diluted observed-aperture measurements, not
physical radius ratios or planet-radius estimates. The fixed trapezoid has
duration 1.727 hours and ingress fraction 0.414, initially inherited from
the sector-1 QLP analysis. Its large ingress fraction and template limitations
make physical modeling essential. Errors take the larger of residual-scaled
diagonal propagation and the scatter of event depths divided by sqrt(N).
They do not provide complete correlated-noise uncertainties.

## Recurrence without assuming an accurate seven-year extrapolation

The discovery-screening period and epoch have limited short-baseline
precision. The sector-92 source test therefore uses an independent search
over all phases within P=1.72649–1.76649 days on a raw target-centered circular
aperture. A 0.5-day running median removes smooth trends separately in
continuous time segments. Only upward 8-MAD excursions are clipped for the
period consistency test; downward events are retained. The original raw
aperture measurements also remain saved without this search-only mask.

The sector-92 regional BLS maximum is P=1.7463145200 days,
T0=3804.4520957294 BTJD, box duration 0.049 days, and diagonal SNR 57.32.
The sector-1 independent pixel-aperture maximum is P=1.7468445200 days,
T0=1326.3465103111 BTJD, with diagonal SNR 36.18. A selected-signal consistency
search is not a blind discovery or a search-wide false-alarm test. Local
period peaks are not period uncertainties and should not be directly
extrapolated between these epochs. Event-timing and cycle-count fits are
needed to establish a joint ephemeris.

The sector-1 local difference image retains its original QLP-derived
P=1.7464945200053 days and T0=1326.3488754694 BTJD. Sector 92 uses its
independently detected epoch/period. Source attribution does not depend on
assuming that the sector-1 epoch predicts the 2025 events accurately.

## Pixel source evidence

The public TESSCut cubes contain un-subtracted sky background. The lowest
30% of median reference-image pixels define a fixed blank-sky mask; their
per-frame median is subtracted from every pixel. All nonzero QUALITY frames
and nonfinite flux/error frames are removed. No original raw file is altered.

Each local +/-0.25-day event window fits a quadratic baseline and a
cadence-integrated transit template simultaneously. Alternative fits add
blank-sky background and empirically measured pointing shifts. Image shifts
use reference-image gradients outside the target's 2.5-pixel region, with
independent source-amplitude and constant terms. No pipeline POS_CORR
information is assumed when those columns are zero. Quadratic-only fits,
pointing fits, broad boxes, and core boxes all place the difference source
near the target.

The full difference-image Gaussian centroids are (7.995, 7.667) in sector 1
and (6.890, 7.727) in sector 92. Their offsets from the TIC/header WCS position
are 0.170 and 0.194 pixels, roughly 3–4 arcsec. **Those offsets cannot be
interpreted using the much smaller conditional bootstrap errors alone.**
Reference-image and difference-image Gaussian fits in matching 3-pixel
regions are separated by 0.032 and 0.094 pixels. Undersampled, asymmetric
TESS PRFs and background/ROI choices shift Gaussian centroids even on the
constant reference source. Fixing a Gaussian at the exact catalog position
worsens diagonal chi-square despite clear reference-source alignment; no
calibrated TESS PRF scene model has been fitted.

Five hundred bootstrap draws per usable sector resample complete event
windows and synchronized residual images in blocks of 10 usable cadences.
The conditional 16th–84th-percentile catalog-offset intervals are
0.149–0.191 pixels and 0.170–0.223 pixels. They omit PRF/WCS systematics,
proper-motion epoch differences, and unresolved-source uncertainty. They
are not astrometric exclusion limits or source-identity probabilities.

All resolved G<17 off-target stars are farther from every bootstrap source
location than the TIC neighborhood. The relatively bright TIC269994635,
Gaia DR3 6806287729705524352, lies ~50 arcsec away and 2.562/2.570 pixels
from the difference centroid. Fixing the empirical Gaussian there worsens
diagonal chi-square by about 2,431/3,869. These are morphology diagnostics,
not calibrated likelihood ratios. Its smallest aperture has no significant
matched dip: -0.7 ± 3.3 ppt in sector 1 and 4.2 ± 3.1 ppt in sector 92.
Larger neighbor apertures overlap the target's wings and show diluted dips;
that expected leakage does not identify the neighbor as their source.

The closest listed neighbor TIC269994643 is ~21 arcsec away and much fainter,
Tmag18.75 versus target12.92. Its full TESS flux is ~0.0047 of the target's,
insufficient to explain a ~2.7% event by a full eclipse unless aperture flux
fractions differ substantially. This is a catalog-flux constraint, not a
measured PRF dilution model. Close unresolved companions remain unexcluded.

The Gaia source cone, current target properties, live recognition checks,
and their provenance are in `data/broad_catalogs/tic269994642_*`. The saved
Gaia cone has 44 rows; an adjacent error JSON records a later display-parser
failure on a missing G value, not a failed archive response. No prior
recognition was identified in the catalog team's finite search, which
does not establish that the signal is previously undiscovered.

## Sector 27 is collateral geometry, not a stellar observation

TESSCut's coverage endpoint lists sectors 1, 27, and 92. The sector-27 cutout
reveals why a coverage listing alone is insufficient: its physical CCD
columns are RAWX20–34, and the projected target column is ~28. These are
prescan/collateral columns outside the illuminated science area. The median
image contains electronic column stripes with flux ~0.05 e-/s and no stars.

An initial search of that low-flux stamp returned an insignificant, misleading
maximum at P~1.7363 days (diagonal SNR4.46). It is an invalid stellar analysis
and contributes **zero** usable target observations or recurrence evidence.
The failed experiment and all data were retained under
`sector27/exploratory_unusable_overscan/` rather than discarded. The current
`sector27/pixel_vetting.json` explicitly records the exclusion and contains
no accepted transit search. The reusable script now checks the target's
physical detector coordinates before detrending or searching and skips this
collateral stamp. Sector 27 cannot bridge the long-baseline ephemeris.

## Preserved data, time validation, and reproduction

Original cutouts remain outside Git in
`/workspace/tess-research-data/broad/pixels/tic269994642/`, about 73 MB decimal
in total. Files, processing versions, data releases, original SHA-256 hashes,
scientific-array hashes, FITS checksums, complete contributing FFI identifiers,
and exact retrieval URLs are recorded in `data/broad_pixel_269994642.json`.
No paid services were used.

Both usable cubes declare TIME=TDB, TIMEREF=SOLARSYSTEM, and BJD reference
2457000. TIMECORR is not added again. The first sector-1 midpoint converts to
DATE-OBS +899.990 seconds, half the 30-minute cadence. Sector 92's
DATE-OBS is 2025-05-07T11:11:19.038 UTC and its first midpoint converts to
2025-05-07T11:12:59.027 UTC, +99.989 seconds, half its 200-second cadence.
The source is original sector-92 SPOC FFIs in archive data release126, not
the affected later QLP99–104 products. Full raw FITS HDU checksum verification
and byte/scientific hashes passed.

The sector directories preserve alternative difference images with errors,
per-event aperture measurements, Gaia pixel positions, bootstrap coordinates,
reference images, pointing/background vectors, folded diagnostics, complete
raw aperture photometry, errors, periodograms, and retrieval identifiers.
`independent_aperture_periodogram.npz` provides unflattened
`raw_aperture_flux_e_s`, `raw_aperture_flux_error_e_s`, and `time_btjd` for
physical fits. The raw errors include quadrature native FLUX_ERR plus an
approximate propagated median-sky uncertainty. Moving-median or local
decorrelation outputs should not be mistaken for uninfluenced raw data.

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 python src/pixel_vet_269994642.py --sectors 1 27 92 --bootstrap 500 --gaia data/broad_catalogs/tic269994642_gaia_nearby.csv
```

Missing cutouts are re-retrieved with a bounded 200-MB archive request.
Cached files are reused. Persistent cache availability in a future task is
not guaranteed. Astrocut regeneration can change header checksum timestamps;
the script preserves the original expected byte hash and separately records
a regenerated byte hash only when all recorded scientific arrays match
exactly and FITS checksums pass. Divergent science arrays cause a stop and
remain available for diagnosis.

Next required checks are joint event timing with cycle-count uncertainty,
physical grazing-transit versus stellar-eclipse fits, alternate-period and
odd/even comparisons, correlated-noise and selection calibration, calibrated
PRF/scene fitting, and fuller novelty/literature review. High-resolution
imaging and spectroscopic mass/binarity constraints would be valuable if
the signal survives those tests. This pixel report does not assign a planet.
