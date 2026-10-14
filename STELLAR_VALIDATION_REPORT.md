# Further falsification and characterization of two eclipsing-system leads

Public data checked 2026-10-14 UTC. **Both repeating eclipse signals survive
the new tests. The background source TIC 273879360 is strongly supported as an
eccentric eclipsing binary, but its eclipse-like variability was already
publicly machine-flagged. The giant-field signal is strongly supported as a
blended eclipsing system; the actual eclipsing star is still unidentified.**
Neither result establishes a newly uncatalogued star, measured component masses,
or a bound multiple system.

The [searchable measurements](results/stellar_systems/deep_validation_measurements.csv)
separate conditional period/model intervals, source-aperture depths and unresolved
source assignments. Current ranked assessments point to these validation dossiers;
the earlier phenomenological measurements remain in their original reports.

| Source | Additional evidence | Revised assessment |
|---|---|---|
| TIC 273879360 / Gaia 4692910966756908288 | Independent native-pixel physical fit; aperture/baseline/limb controls; archived ASAS-SN eclipses | Strong eccentric binary interpretation; independent characterization of a previously flagged source |
| TIC 260352037 giant field | 16 additional QLP sectors; four native cubes; phase, parity, timing and source-model controls; Gaia spectroscopy | Real coherent eclipse blend; ordinary detached orbit on bright giant is implausible; exact source and binding unresolved |

## TIC 273879360: physical model and attempts to break it

The source is 65.59 arcseconds from the foreground TIC 273879358. Earlier difference
images place both eclipse families on the background source in sectors 27/29/67/95.
The new extraction starts again from the original checksum-verified native
pixel cubes. A one-pixel-radius aperture centered on the Gaia source, a
time-independent blank-sky mask and QUALITY=0 selection produce 22 local eclipse
windows and 6234 measurements. A stricter gap/coverage rule excludes one window
used by the previous 23-event phenomenological analysis; all exclusions are
recorded. No retained raw data were removed.

Two separate limb-darkened occultations represent two luminous spheres with
a shared eccentric two-body orbit. Their secondary time is tied to Keplerian
geometry, using the near-edge-on conjunction approximation. The fit integrates
each native exposure, profiles a quadratic raw baseline for every event, and
allows a constant binary-light fraction in each sector. Limb-darkening
coefficients are approximate fixed choices, tested below. It omits calibrated
pixel-response/color effects, spots, reflection and tidal distortion.

| Conditional fitted quantity | Value |
|---|---:|
| Period | 7.586890837 days |
| Primary epoch, BJD(TDB)−2457000 | 1327.900292 |
| Secondary phase | 0.691727 |
| Eccentricity | 0.32145 |
| Argument of periastron under this component convention | 341.24 degrees |
| Inclination | 88.48 degrees |
| Radius ratio R2/R1 | 0.680 |
| Semimajor axis/R1 | 20.63 |
| TESS luminosity ratio L2/L1 | 0.145 |
| Surface brightness ratio | 0.314 |

Eight distinct starting geometries converge to the same chi-square ≈ 6037.
Dark-companion and circular-binary alternatives fail strongly in these targeted
windows, with conditional chi-square penalties ≈ 22788. These numbers are model
comparisons under empirical errors, not calibrated survey false-positive
probabilities or evidence for first discovery. A luminous secondary naturally
explains the roughly 7–9% source-aperture secondary eclipse.
The circular and dark-model penalties mostly express the same failure to
produce a secondary in these selected windows; they are not independent
pieces of statistical evidence. Exact minimum projected separation differs
from the model's conjunction approximation by about 0.68 minutes, negligible
for this eccentricity interpretation but relevant to precision timing.

Raw linear/cubic baselines, 0.8/1.2-pixel apertures and weaker/stronger limb
darkening give eccentricity 0.3204–0.3224. Radius ratio 0.663–0.696 and light
ratio 0.142–0.148 change more with these choices. Increasing exposure integration
from 15 to 61 samples changes the model by at most 39 ppm, far smaller than the
percent-level eclipses. Twenty-four 30-minute residual-block resamples yield
conditional 16th–84th percentile ranges: e 0.3197–0.3222, radius ratio 0.664–0.738,
and light ratio 0.141–0.166. The small resample count and omitted inter-event
calibration/stellar-model uncertainty preclude treating these as total errors.
The wider bootstrap radius range illustrates the dilution/geometry degeneracy.

The secondary timing alone requires e ≥ approximately 0.306 for an edge-on
two-body orbit. Duration asymmetry and the physical shape constrain a plausible
family near e 0.32. The combined-source catalog radius is not R1: assigning its
approximately solar size to one component would silently convert an unresolved
stellar model into a component measurement. No absolute radii, masses or
mass ratio are claimed.

The half-period hypothesis would require extra primary and secondary eclipses
at phases 0.5 and 0.1918 of the 7.5869-day cycle. Fourteen and twelve covered
control windows show no such flux loss: event-mean template depths are
0.067% and −0.095%, compared with approximately 22–27% and 7–9% predicted
losses. Event-resample intervals and signed individual measurements are
preserved. This rejects an ordinary single binary repeating at half the
period; it does not exhaust contrived multiple-source configurations.

A deliberately failed early-to-later prediction is also retained. The early
sample has only four primary events. Its original restricted period fit did
not converge and approached a bound; fixing that output as precise produced
a large later timing penalty. Widening the bounds and profiling the early
period reveals much poorer timing precision than the full multiyear sample.
The preferred period shifts by about seven seconds from the full-data value,
and that modest difference accumulates over years. A small residual
early/late timing or shape discrepancy remains; the failed nonconverged fit
cannot establish orbital inconsistency. See the converged likelihood-profile
results and initialization diagnostics rather than interpreting its original
9034 chi-square penalty as physical evidence.
The final profile has conditional delta-chi-square=1 crossings at
7.58676659–7.58685246 days and delta=3.84 crossings at
7.58672511–7.58689416 days. The full-data period lies within the latter
interval (training penalty 3.54); all final profile fits converge. Holding the
full-data period while refitting early shape reduces the later prediction
penalty to 50.5. This is a diagnosis using an already measured period, not a
new independent prediction, and the remaining penalty warrants caution.

See [physical fit and diagnostics](candidates/stellar_systems/tic273879360_validation/README.md).

## Independent earlier observations and the novelty correction

Archived ASAS-SN source AP36443078 matches the background coordinates within
0.158 arcseconds and has 552 V-band observations from 2014–2018. CSV and embedded
HTML measurements agree exactly. At the fixed TESS ephemeris, eleven primary
measurements represent nine separate events, eight before TESS. Every event
is fainter after alternative local/camera baselines. A descriptive primary
template amplitude is 35.2%, with conditional night-block 95% interval 28.4–39.7%;
leave-one-event-out and baseline alternatives retain the signal. Different
bandpass, source coupling and fixed shape prevent interpreting this as the
same physical depth as the TESS aperture.

The secondary is only supportive: its fitted 6.7±2.2% amplitude is
camera-dependent. Fifty resolved SkyMapper epochs miss both eclipse windows,
so they neither confirm nor rule out the signal. Archived ASAS-SN HJD values
use field centers and approximately 200-second timing accuracy. They cannot
extend a high-precision BJD(TDB) ephemeris without corrected time standards
and source-specific light-travel-time conversion. Public SkyPatrol maintenance
blocked its current light-curve route; the older release remained accessible.

**Our earlier novelty wording was too strong.** Kostov et al.2025 Table 2 is
a published collection of roughly 900000 targets with neural-network scores> 0.9
for eclipse-like features. It is not a generic observed-target inventory.
TIC 273879360 is in that collection. Its row supplies no vetted classification
or period, and the vetted Tables 3/4 contain no match, but prior public eclipse
recognition exists. This work adds independent source attribution, ground
recurrence and orbital characterization. It cannot claim first recognition
or prove first published binary classification. The dated correction is
recorded in the research log and all current ranking documents.

Gaia context places the combined background source at an inverse-parallax
scale 812±11 pc, with dwarf-like luminosity/gravity. TIC and Gaia single-star
temperature/radius models disagree enough to discourage component estimates.
Its clean Gaia photometric/astrometric diagnostics and lack of convincing
warm infrared excess do not exclude a close binary. The large Gaia RV error
is a combined-cross-correlation uncertainty, not epoch scatter or an SB
detection. Gaia MSC's two-star fits are supplied under an assumed binary
model and explicitly do not classify whether a source is binary.

Read the [ground evidence](candidates/stellar_systems/deep_context/ground_photometry/README.md)
and [catalog/stellar context audit](reports/stellar_systems/deep_context.md).

## TIC 260352037: signal strengthened, source assignment weakened

Sixteen genuinely later 200-second QLP sectors bring the total to 41 across
2018–2026,7.45 years. The old 0.880788466-day ephemeris was fixed before fitting
these data and correctly predicts both eclipse families. Independent raw
linear/quadratic/cubic baselines preserve them. In the quadratic extraction,
962 primary and 967 secondary events give central template depths 1696±27 and
784±28 ppm. The new sectors alone give 1669±40/801±30 ppm. These central depths
are distinct from the old 972/445 ppm in-eclipse window averages.

Sector timing profiles give conditional P=0.880788521±0.000000330 days,
primary epochBTJD 1326.593182±0.000543. Later shapes last approximately 1.87 and
1.68 hours. Alternating primary and secondary contrasts at twice the period
are−49±45 and+8±46 ppm; no extra alternation is measured. A quadratic timing
term gives period change+0.004±0.030 seconds/year, with no established drift.

An initial unmasked quarter-phase control was flawed: its baseline included
real nearby eclipses and manufactured 266–311 ppm dips. After masking both
families, ten control phases give−48..+28 ppm with 30–39 ppm errors. The failed
control is retained. Synchronized four-cycle sign-reversal blocks also fail
to reproduce the observed means in 512 trials; this targeted check assumes
symmetric errors and supplies no global discovery probability.

Four native cubes, including new sector 95, put the deficit in the central
blend. Empirical translated reference profiles, event resampling, ROI changes,
background and pointing variants **do not securely identify the eclipsing
star**. Structured profile residuals and position changes ≈ 0.2–0.25 pixel exceed
the separation of the proposed close source. Nominal source likelihoods or
bootstrap centroids would overstate confidence without calibrated pixel
response and systematic errors. The primary and secondary deficits cannot
be securely declared coincident or on different sources from these fits.

Independent Gaia RVS spectroscopy finds T_eff ≈ 4834 K, logg ≈ 2.51 and radius
approximately 9.45 solar radii for the visible star. Its catalog RV-constancy
p ≈ 0.807 shows no reported significant variation of that bright giant. At the
measured period, illustrative 1–3 solar-mass semimajor axes ≈ 3.9–5.6 solar radii
are smaller than the giant. An ordinary detached eclipse orbit on the giant
is physically implausible; assigning it a dwarf radius cannot rescue it.

TIC 737466840 is 4.215 arcseconds away. Equal aperture coupling gives a faint-light
fraction 0.006685, so the observed central primary requires approximately 25%
intrinsic eclipse loss, plausible for a faint stellar binary. Its contaminated
Gaia colors prevent precise stellar component estimates. TIC 260352041 at
7.73 arcseconds requires approximately 111% intrinsic loss under equal coupling,
but an 11% coupling change removes that contradiction; it is disfavored rather
than securely excluded. Neither neighbor has been positively identified as
the eclipse source.

Compatible parallax/proper motion does not prove that 737466840 is bound to
the giant. At the nominal 7057 AU projected separation, the central transverse
speed 1.086 km/s implies a minimum binding mass 4.70 solar masses, awkward for
a simple giant-plus-faint-dwarf picture. Uncertain proper motion, unmeasured
radial motion and Gaia cross-source systematics prevent a binding probability
or an unbound conclusion. No new triple system is established.

Read the [full giant-field validation](reports/stellar_systems/tic260352037_validation.md),
including data products, uncertainties and failed source controls.

## Reproduction, preservation and next observations

All original raw observations remain retained. New giant observations and
cube add 38.3 MB outside Git;45 old/new raw products and 58 outputs have a
106-entry method/data/output verification manifest. Background raw inputs
are verified during each extraction. Compact fit arrays, plots, measured
events, catalogs, public retrieval manifests and failed responses are tracked.
The cloud overlay has no independently verified durable external volume;
exact URLs, identifiers and hashes permit restoration. No paid models or
external submissions were used; research OpenRouter expenditure remains$0.

Use the [resumption checkpoint](PROGRESS.md) and source-specific dossiers for
commands. The decisive next observations are spatially resolved primary and
secondary photometry of the giant's neighbors, and multiband eclipses plus
mass-sensitive spectra of TIC 273879360. Public data have supplied strong
repeating-signal evidence but cannot replace these observations. The earlier
reserve queue remains a separate route to finding a less-recognized system.
