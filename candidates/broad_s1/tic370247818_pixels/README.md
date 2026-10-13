# Independent pixel source checks: TIC 370247818

Checked 2026-10-13. **The ~3.5222-day dip is independently recovered in FFI
apertures and localized to the target's neighborhood in three observing
geometries. The bright star CD-79 151, 86.6 arcsec away, is strongly disfavored
as its source. These tests do not establish a planetary origin or exclude an
unresolved blended eclipsing binary.** Absolute source attribution remains
limited by the empirical Gaussian image model.

The data are independently extracted public TESSCut FFI pixels, rather than
the QLP aperture columns. They are independent extractions of the same sky
observations, not statistically independent of the QLP detections in these
sectors. Sectors 5 and 39 additionally check different epochs and camera/CCD
geometries. There are 22 usable local event windows in these three cutouts.

| Sector | Camera / CCD | Local windows | Target r=2 px diluted depth | Conservative aperture SNR | Difference/WCS offset | Difference/bright-neighbor distance |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 3 / 3 | 8 | 0.00681 ± 0.00058 | 11.83 | 0.096 px | 4.187 px |
| 5 | 4 / 3 | 6 | 0.00802 ± 0.00098 | 8.22 | 0.079 px | 4.538 px |
| 39 | 3 / 4 | 8 | 0.00773 ± 0.00040 | 19.54 | 0.331 px | 4.746 px |

Depth errors take the larger of the propagated local residual-scaled diagonal
error and scatter of per-event depths divided by sqrt(number of windows).
They are approximate, not a full correlated-noise posterior. Aperture flux
includes diluted target and neighbor wings; these depths are not stellar
radius-ratio measurements. The integrated trapezoid has duration 2.254 hr and
fixed ingress fraction 0.331. It is not a physical transit model.

## Method and falsification

All nonzero `QUALITY` cadences and nonfinite pixels are excluded. TESSCut
`TIME` is already barycentrically corrected, in TDB at BJD reference 2457000;
`TIMECORR` is not added again. Thirty-minute cadences in sectors 1/5 and
ten-minute cadences in sector 39 are integrated when computing the template.

For each event, a +/-0.30-day local fit jointly estimates the positive
out-minus-in transit image and a quadratic time baseline for every pixel.
Two alternative fits add the contemporaneous blank-sky background and
background plus measured x/y image shifts. `POS_CORR1/2` are zero in these
cutouts, so image shifts are estimated from the median reference image and
its gradients with an independent source-amplitude parameter. The target
region is excluded from this registration. Simultaneous nuisance and transit
terms avoid explicitly fitting the target transits into an out-of-transit
trend. All three methods give a compact positive difference near the target.

Difference images combine events with pixel inverse variances. An elliptical
Gaussian plus affine background is fitted inside a target-centered 3-pixel
region for a conditional position diagnostic. Five hundred bootstrap draws
resample whole events and synchronously resample residual images in blocks of
three cadences. Thus spatial covariance is retained, and some local temporal
noise is represented. Every draw in all three sectors is closer to the target
than to CD-79 151. This is not a posterior probability of source identity.

Independent apertures at the target, CD-79 151, an intermediate position,
and the third bright star reproduce the dip at the target. At r=2 pixels,
the bright-neighbor event depths are 0.000178 ± 0.000873 (sector 1),
0.000303 ± 0.000200 (sector 5), and -0.000020 ± 0.000076 (sector 39).
This is inconsistent with an obvious matched deep eclipse on that resolved
bright source. The intermediate aperture has much smaller fractional dips.

The first full-sector BLS consistency test deliberately retained all usable
`QUALITY=0` frames. Sector 1 then produces a **spurious 4.105-day maximum** in
both target and bright-neighbor apertures, caused by common image disruptions.
QUALITY flags alone are insufficient. An ephemeris-independent 6-MAD clipping
criterion on the median-detrended bright comparison aperture removes 46 of
1,267 usable sector-1 frames, including five transit-window points belonging
to the particularly noisy cycle 6. The corrected target period maximum is
3.522625 days with diagonal BLS SNR 19.66. The unmasked original curves and
periodograms remain saved. The comparison-star mask must not be used to
search for eclipses on that comparison star itself.

Sector 5 and sector 39 recover the signal naturally with all usable frames:
regional BLS maxima are 3.521500 days (diagonal SNR 16.46) and 3.523000 days
(diagonal SNR 19.46), respectively. Comparison-star clipping changes neither
period; it removes five sector-5 frames and zero sector-39 frames. These
2.8–4.3-day searches test period consistency. Their scores are neither
survey-wide false-alarm probabilities nor independent confirmations.

## Localization uncertainty and unresolved blends

The 3-pixel-ROI bootstrap conditional position offsets from catalog/WCS are
0.063–0.145 pixels (sector 1), 0.053–0.136 pixels (sector 5), and
0.303–0.361 pixels (sector 39), at the 16th–84th percentiles. **These intervals
do not cover Gaussian/PRF-model or WCS systematics.** Sector 39's apparent
~7-arcsec catalog offset cannot be interpreted as an accurately localized
separate source.

The Gaussian result depends on aperture radius and background modeling.
Larger reference-image fitting regions include the much brighter neighbor's
halo and can become unphysical; width-boundary hits are recorded, rather
than treated as reliable centroids. With matching 2-pixel regions, the
reference-image and difference-image centroid separations are 0.015 pixels
(sector 1), 0.130 pixels (sector 5), and 0.103 pixels (sector 39). In sector
39 the observed reference star also has a substantial catalog/WCS offset.
The sector-39 2-pixel difference fit reaches a Gaussian width bound, so its
centroid agreement is particularly model dependent. This supports the
target neighborhood, but demonstrates why the smaller
conditional bootstrap errors cannot be used as absolute astrometric limits.
See each `localization_roi_sensitivity.csv` for every fitting radius.

Gaia positions are projected with each cutout's WCS. Proper motion is not
propagated here, and a Gaia source cone can miss unresolved stars. A calibrated
TESS PRF scene model, high-resolution imaging, or ground-based transit
photometry is needed to exclude close blending. Stellar eclipses on the
target itself also remain possible. Novelty catalogs, odd/even and secondary
checks, physical fitting, and statistical calibration are assessed elsewhere
in the candidate's main report; pixel localization alone answers none of them.

## Preserved evidence and reproduction

- Root directory: sector-1 JSON, per-event aperture CSV, Gaia pixel positions,
  FFI identifiers, difference-image arrays, extracted photometry, bootstrap
  arrays, BLS arrays, and diagnostic plots.
- `sector5/` and `sector39/`: the same products for independent geometries.
- `../../../data/broad_pixel_manifest.json`: original byte hashes, science-array
  hashes, archive API requests, processing versions, cadence/time metadata,
  and passing FITS CHECKSUM/DATASUM results.
- All three raw cutouts remain outside ordinary Git tracking in
  `/workspace/tess-research-data/broad/pixels/`. They are never deleted by the
  analysis script. Persistence of this local cache between future tasks is
  not guaranteed.

TESSCut retrieval produces a ZIP containing a FITS cutout. Future downloads
may contain newly generated FITS checksum timestamps, causing the original
byte hash to differ even if data are unchanged. The manifest also provides a
scientific array hash computed from labeled, little-endian TIME(float64),
QUALITY(int32), FLUX/FLUX_ERR(float32), and 38-byte FFI_FILE strings, in file
order. Preserve divergent downloads for diagnosis; do not silently replace
the original product or expected hashes.

Run with the repository's Python environment, from its root:

```bash
OPENBLAS_NUM_THREADS=1 python src/pixel_vet_370247818.py --bootstrap 500
OPENBLAS_NUM_THREADS=1 python src/pixel_vet_370247818.py --sector 5 --period 3.522194945921 --epoch 1326.9162961298 --bootstrap 500
OPENBLAS_NUM_THREADS=1 python src/pixel_vet_370247818.py --sector 39 --period 3.522194945921 --epoch 1326.9162961298 --bootstrap 500
```

The sector-1 analysis retains its initial short-baseline ephemeris
P=3.522500444983 days, T0=1326.915994193442 BTJD. Later independent epochs use
the refined joint ephemeris P=3.522194945921 days, T0=1326.9162961298 BTJD.
The tiny change within sector 1 does not invalidate its source test. The
script accepts `--raw` for a re-retrieved cutout path.

Next scientific checks: calibrated PRF localization relative to the target
and all neighboring Gaia sources, broader catalog/literature novelty checks,
physical fits with dilution constraints, and observation-based exclusion of
unresolved stellar eclipses. A transit on a star is not by itself a planet.
