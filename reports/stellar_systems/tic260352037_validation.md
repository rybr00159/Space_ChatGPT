# Falsification and characterization of the eclipse blend near TIC 260352037

Checked 2026-10-14. **The repeated eclipse signal survives substantially more
stringent tests, but the identity of the eclipsing star remains unresolved.**
This strengthens the case for a stellar eclipsing system in the giant's blended
TESS image. It does not establish a new component of the giant, a bound triple,
or first discovery of a particular stellar system.

## Independent later observations

Sixteen additional public QLP products were retained: sectors 61–66, 68–69,
88–90, 93, and 95–98. These have native 200-second exposures. Sector 67's
predicted product returned 404 and that failed retrieval is retained in the
manifest. The TessCut coverage query has no coverage for this position in
sectors 74–79 or 99–104, avoiding the known affected QLP v01 releases. No
SPOC observations were returned by the exact-TIC MAST target-name query; that
finite query is not proof that no other independent reduction exists.

The original 25 sectors plus the new 16 span 2018-07-26 to 2026-01-05, or
7.45 years. Before looking for a new optimum, the old ephemeris was held fixed:

- Recurrence period: 0.880788466 days.
- Primary epoch: BTJD 1326.593218, with BTJD = BJD(TDB) − 2457000.
- Primary template duration: 0.07718209 days, including each native exposure.

Both eclipse families persist at these predicted phases in the new data.
Raw SAP flux was fitted independently around each event with a simultaneous
polynomial baseline and a cadence-integrated trapezoid. These are **central
trapezoid depths**, so they should not be compared directly with the earlier
972/445 ppm means of all in-eclipse samples, which include shallow ingress.

| Raw baseline | Original 25 sectors, primary / secondary | New 16 sectors, primary / secondary | All 41 sectors, primary / secondary |
|---|---:|---:|---:|
| Linear | 1729 ± 35 / 801 ± 41 ppm | 1649 ± 40 / 817 ± 32 ppm | 1691 ± 26 / 809 ± 28 ppm |
| Quadratic | 1722 ± 37 / 764 ± 43 ppm | 1669 ± 40 / 801 ± 30 ppm | 1696 ± 27 / 784 ± 28 ppm |
| Cubic | 1719 ± 35 / 765 ± 43 ppm | 1674 ± 40 / 801 ± 30 ppm | 1697 ± 26 / 784 ± 28 ppm |

The quadratic results contain 962 covered primary and 967 covered secondary
events; 435 and 433 are genuinely later holdout events. Errors take the largest
of residual/diagonal and five-cadence correlation estimates per event, then
include observed event scatter and sector heterogeneity. They remain targeted,
conditional errors, not global discovery significances. Sector depth changes
can reflect aperture dilution, giant variability, and template mismatch.

## Attempts to falsify the period and eclipse interpretation

**Cadence and detrending.** The later 200-second data resolve the dip shapes;
the events are not produced solely by 30-minute integrations. Independent
linear, quadratic and cubic raw baselines give similar results. The two-day
median-normalized fold is a display and descriptive morphology check; the
raw local measurements above do not require that normalization.

**Control phases.** Unmasked quarter-phase control windows initially generated
apparent 266–311 ppm dips. Their baseline windows include the real neighboring
eclipses and distort the polynomial: this was a failed control construction,
not another discovery. After masking the actual primary and secondary events,
ten later-sector control phases give depths between −48 and +28 ppm with
30–39 ppm errors. The failed and corrected controls are both preserved.
They share data and are not an independent null population.

A second conditional noise check reverses signs in synchronized four-cycle
blocks, preserving primary/secondary pairing and event magnitudes. In 512
trials none approaches the observed primary or secondary mean; maximum
absolute control means are 513 and 301 ppm, respectively. The validity of a
sign-reversal check requires symmetric errors, and the signal was already
selected, so this does **not** supply a survey-wide false-alarm probability.

**Half and double periods.** Folding at 0.440394 days merges unequal dips.
The persistent primary/secondary difference makes that a poor single-eclipse
recurrence interpretation. At 1.761577 days, alternating primary events differ
by −49 ± 45 ppm and alternating secondary events by +8 ± 46 ppm after
sector offsets/scatter. There is no evidence that two 0.8808-day cycles differ.
An exactly repeating 2P configuration remains a mathematical photometric alias;
radial velocities, rather than a parity nondetection, determine a secure orbit.

**Timing stability.** Sector timing profiles include all covered events, rather
than selecting only individually significant events. A linear refit gives:

- P = **0.880788521 ± 0.000000330 days**.
- Primary epoch = **BTJD 1326.593182 ± 0.000543**.

These use fixed templates, a conservative three-minute sector correlation
floor and 1,000 whole-sector resamples. The original and later observations
separately give 0.880788193 ± 0.000001133 and
0.880788449 ± 0.000001148 days. A quadratic timing term corresponds to
+0.004 ± 0.030 seconds/year: no measurable secular drift is established.
All uncertainties are conditional; no eclipse timing dynamics claim is made.

**Shapes and stellar variability.** Refitting the two shapes independently in
later sectors gives primary duration about 1.87 hours and secondary about
1.68 hours. Joint central depths are 1713 ± 47 and 777 ± 35 ppm before
additional sector/model uncertainty. The secondary is near half phase; separate
sector estimates allow a few minutes' offset, and no significant eccentricity
is inferred. The short isolated dips remain coherent through 7.45 years while
the giant's broad stochastic variability changes. A spot or ordinary giant
oscillation interpretation is substantially less plausible than an eclipse
blend. Phenomenological template fits cannot exclude every nonlinear variable
star model or determine component radii/temperatures.

## Pixel evidence and an unsuccessful source attribution

Four raw native TessCut cubes (sectors 1, 28, 39, 95) were checked for actual
science-silicon coverage and reprocessed using the same frozen old ephemeris.
The new sector-95 cube contains 26 primary and 24 secondary covered events.
In all four sectors, the main deficit is in the central blended image, rather
than on the bright resolved star near the stamp edge.

A new relative astrometric test uses the observed giant-dominated reference
image as an empirical profile. It fits translated reference images to each
eclipse deficit, rather than assuming a symmetric Gaussian. Gaia positions
are propagated from their reference epoch using proper motions. The tests
retain background and pointing nuisance variants and 2.5, 3.5 and 4.5-pixel
regions. Whole-event resampling gives statistical centroid dispersion.

This experiment **does not identify the source**. The reference and deficit
profiles have systematic differences, including deficits narrower than the
reference and structured residuals in its wings. Subpixel interpolation,
neighbor light, source color and intra-pixel response remain uncalibrated.
Region choices change the fits. In sector 28 the primary relative horizontal
shift changes from about +0.090 to −0.167 pixel when pointing proxies are
included; in sector 95 it changes from +0.139 to −0.113 pixel. These proxies
are estimated from the whole image and may absorb genuine eclipse changes,
so neither variant is a calibrated correction. Some low-S/N sector-1 secondary
fits hit their position bounds. The nominal diagonal likelihood differences
between proposed sources must not be treated as decisive source odds.

The primary and secondary centroids also differ somewhat within sectors.
Their weakness, background dependence and template residuals prevent a secure
claim that they come from different stars. Equally, these data do not prove
that the two deficit images originate on precisely the same component.

Therefore the 4.21-arcsecond and 7.73-arcsecond neighbors, a closer unresolved
source, or a more complex blend remain live hypotheses. A simple centroid
comparison cannot identify a new inner binary in either catalogued neighbor.

## Physical context and dilution constraints

The independent context audit finds Gaia GSP-Spec Teff approximately 4834 K,
log g approximately 2.51, and a FLAME spectroscopic radius approximately
9.447 solar radii for the bright giant. Its stable Gaia radial velocities
(p-value for constancy 0.807, 21 transits) do not reveal a giant short-period
spectroscopic orbit. The 0.8808-day separation for illustrative 1–3 solar-mass
systems is only about 3.9–5.6 solar radii. Merely achieving a separation equal
to the giant radius would require about 14.6 solar masses; a Roche-safe orbit
is more restrictive. A normal detached eclipse on this giant is physically
implausible. The same data support a genuinely large visible star, so simply
assigning a dwarf radius to TIC 260352037 does not solve the blend.

Current TIC T magnitudes are 11.169 ± 0.006 for the giant,
16.5989 ± 0.0076 for TIC 737466840 and 18.1961 ± 0.0411 for TIC 260352041.
For equal aperture coupling, with just giant plus proposed source, their
fractions of total T-band light are 0.006685 and 0.001543. A 1710 ppm primary
central loss would require about **25.6%** intrinsic eclipse depth on
TIC 737466840, physically plausible for a stellar pair, or **110.8%** on
TIC 260352041, nominally impossible. The latter needs at least approximately
11% greater relative coupling or a changed flux estimate; neither is excluded
by these undersampled images. Additional contaminating light worsens the
required intrinsic depth. TIC errors do not include all color/crowding errors:
737466840's Gaia BP/RP photometry is contaminated, so accurate component
spectral models from those magnitudes would be unwarranted.

If two similarly sized occulting regions generate both eclipses, their
central depth ratio around 0.45–0.46 can reflect unequal surface brightness.
That is a conditional geometry interpretation, not a measured temperature
ratio: dilution, grazing geometry and limb darkening remain uncertain.

The possibility of a wide association with TIC 737466840 also deserves
skepticism. At the giant's nominal distance, 4.215 arcseconds corresponds to
7057 AU. The relative proper-motion vector (−0.13676, −0.00596) mas/year
implies nominal transverse speed 1.086 km/s. Bound motion would require total
mass at least 4.70 solar masses using only projected separation and transverse
speed, which is awkward for a straightforward giant plus faint dwarf picture.
Within-source covariance and astrometric errors are important: illustrative
Monte Carlo 16–84 percentile minimum masses span 2.21–16.75 solar masses.
Those transformed error draws are **not** posterior binding probabilities;
common distance, radial motion, true separation and cross-source Gaia
systematics remain unknown. Similar astrometry alone does not prove binding.

## Recognition, preservation and follow-up

Recognition/novelty results are covered by the [prior audit](prior_audit_supplement.md)
and [deep-context audit](deep_context.md). No new catalog absence is promoted to proof of first discovery here.
There is convincing evidence for an eclipsing system in this field, but its
host identifier is not secure enough to register a specific new binary solely
from these TESS data.

The decisive next observation is resolved differential photometry of the giant
and both neighbors across one predicted primary and one secondary event. The
intrinsic primary on 737466840 could be approximately 25% deep, vastly easier
to detect than the diluted 0.17% signal. That outcome can also directly falsify
737466840 as the source. Spectroscopy becomes useful after spatial source
identification. Calibrated TESS PRF scene fitting is another option, but the
empirical-image experiment here does not substitute for that calibration.

Code: `src/stellar_giant_validate.py`. Outputs:
`candidates/stellar_systems/tic260352037_validation/`. Exact later product URLs,
versions, byte hashes, FITS checksums and failed retrievals:
`data/stellar_systems/tic260352037_validation/later_observations.json`.
All new raw light curves and the sector-95 cube are retained outside Git in
`/workspace/tess-research-data/stellar_systems/tic260352037_validation/`.
Original raw cubes/light curves remain in their historical cache locations.
Each pixel directory records raw URLs, byte hashes, full relevant headers,
science detector position, FFI identifiers and derived per-event deficits.
A separate externally durable storage volume remains unverified.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/workspace/openrouter-check-venv/bin/python scripts/stellar_restore_lightcurves.py --download
/workspace/openrouter-check-venv/bin/python src/stellar_giant_validate.py --mode retrieve
/workspace/openrouter-check-venv/bin/python src/stellar_giant_validate.py --mode photometry
/workspace/openrouter-check-venv/bin/python src/stellar_giant_validate.py --mode morphology
/workspace/openrouter-check-venv/bin/python src/stellar_giant_validate.py --mode pixels
/workspace/openrouter-check-venv/bin/python data/stellar_systems/tic260352037_validation/verify_manifest.py
```

Retrieval modes retain existing observations and reject changed expected raw
hashes. They retrieve missing observations when the environment cache is lost;
original old-sector restoration follows its retained manifest. `morphology`
also writes timing-stability and wide-pair-dynamics checks. The verified checkpoint lists 45 raw products (41 light curves and four cubes),
checks three method hashes and 58 derived artifact hashes, and retains the
16 later light curves plus the new cube in 38.3 MB of additional raw cache.
A fresh-download
full reconstruction of this new dossier has not been independently exercised;
all present-cache inputs were integrity checked during the analyses.
