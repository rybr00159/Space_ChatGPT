# Multi-sector lead falsification

This is an evolving follow-up checkpoint, not a planet catalogue. The parent survey searches fixed samples independently of these outcomes. All raw observations, rejected signals and conditional diagnostics are retained. Warnings remain available for joint searches: a warning alone does not establish a false positive or exclude an overlapping planet.

## Methods and scope

`src/multisector_followup.py` verifies raw SHA-256, FITS integrity, explicit BTJD/TDB timing and the TIC/sector identity. It supports the newer QLP `DET_FLUX` and older `KSPSAP_FLUX` columns without modifying frozen broad-survey code. Source snapshots under `candidates/multisector/method_snapshots/` preserve the exact code hashes attached to outputs.

Selected periodic signals are locally refined and checked using corrected flux, raw SAP flux and both extra apertures. Diagnostics include 0.75/1.5-day gap-safe median detrending, raw local quadratic baselines, individual event depths and scatter, alternating events at twice the period, half-phase and searched additional-phase events, orbital variability, centroid changes, and proximity to flagged cadences. The local polynomial errors include residual floors and an event-scatter floor. These errors do not calibrate the survey false-alarm probability.

Independent QLP retrieval uses geometric TESSCut coverage and exact public archive URLs. The distinction between geometric coverage and an available extracted light curve is preserved. Sectors 74–79 use the corrected `v02` release; sectors 99–104 are excluded because of the documented QLP processing warning. Successes, unavailable products and error response bytes remain cached. A targeted common-epoch period search across available sectors is follow-up optimization, not an independent blind detection or a calibrated significance.

Large signals were examined as controls for false alarms, while modest-depth and longer-period leads were also tested. This queue is scientifically selected follow-up, not an additional unbiased sample or a completeness measurement. `candidates/multisector/followup_summary.csv` enumerates the detailed reports.

## Candidate and false-alarm investigations

### TIC 368248545: coherent recurrence with a concrete known-EB alternative

The sector-86 signal initially has period 1.124589 d and corrected/raw local S/N 10.47/9.92. It recurs in independent sectors 73, 59 and 19. A common-epoch refinement gives 1.1241347 d and corrected/raw local S/N 10.34/9.47, 12.64/12.39, 10.09/9.71 and 6.44/6.44 in sectors 86, 73, 59 and 19, respectively.

This recurrence does not make it a planet. The dominant out-of-event frequency is approximately 6.22896 cycles/day (0.1605405 d), around seven times the putative orbital frequency. Multiple other troughs occur per folded cycle. A Fourier baseline fitted outside the proposed events removes the purported primary; the eight-harmonic residual is a brightening rather than a transit. A flexible baseline alone cannot exclude a coexisting planet.

Fresh live context identifies the nearby known EW binary ZTF J043919.48+501122.7 / ATO J069.8311+50.1896, Gaia DR3 260154557196225536, 24.25 arcsec away, with a published VSX period 0.321179 d. Its half-period agrees with the dominant TESS modulation to approximately 0.03%. The QLP centroid warnings recur in independent sectors. Native 9×9 TESSCut pixels now localize the dominant modulation approximately 1.194 pixels from the catalog target and 0.019 pixels from the known EW source. An alternative background/baseline treatment gives 1.193 and 0.010 pixels, respectively. The Gaussian source fit is stable across these treatments and the fixed-neighbor position fits much better than the target position (conditional diagonal χ² 63 versus 2,179 with daily baselines). This identifies the known neighbor as the source of the dominant modulation, subject to PRF/WCS and correlated-noise limits; the χ² difference is not a calibrated significance. It does not exclude a weaker overlapping transit.

A simpler three-parameter single-frequency sinusoid fitted outside the putative events leaves only corrected/raw local S/N 1.32/1.52. This provides a less flexible alternative to the eight-harmonic demonstration.

The earlier pulsation hypothesis remains recorded in messages and diagnostics, but this known eclipsing source is a stronger alternative. The catalog target itself is Gaia DR3 260154552898141312. See the catalog review's preserved raw queries and `results/multisector_recognition/live_context_summary.json` for exact dates, coordinates and recognition scope.

### TIC 302550769: a basic-check survivor that did not survive broader scrutiny

The sector-86 lead at 1.350214 d has corrected/raw local S/N 8.24/8.09 and no initial automatic warning. Visual inspection revealed repeated waves across the folded cycle. The strongest out-of-event frequency is 4.45276 cycles/day, or 0.224580 d, close to six times the proposed orbital frequency; the original two-harmonic variability check did not capture this higher-frequency alternative.

Actual QLP observations were retrieved for sectors 85, 79, 78 and 59. At the targeted common-epoch optimum 1.3475406 d, the original sector-86 raw local S/N is −0.08; the other sectors give only 1.93–4.42. Individually optimized attractive folds therefore do not establish coherent recurrence of the proposed primary. A planet is not excluded by these limited tests, but this is not a credible new candidate.

The exact catalog target, Gaia DR3 521607117926177408, already has a Gaia non-single-star classification. It has RUWE 2.284, published RV variability (amplitude 21.82 km/s; constant-RV p=0.00324), and a published NSS astrometric `Orbital` solution at 161.673704 ± 0.333587 d (eccentricity 0.14239 ± 0.05391). The solution flags are 64, efficiency 0 and goodness of fit 6.34; the flags have not been decoded or the solution independently validated. This is a prior orbital report, not proof of a confirmed binary or any companion mass. These properties do not assign the 0.22458-day photometric signal or prove that no planet exists; they are important prior context that a core exoplanet-catalog nonmatch missed.

### TIC 261083172: weak initial signal fails to gain independent support

A sector-83 peak near 5.725 d has initial screening local S/N 9.95, but independent local corrected/raw estimates are 3.27/2.26. Actual sectors 56 and 16 were retrieved; sector 15's unavailable-product response is retained. The targeted common-epoch optimum 5.7298253 d yields corrected/raw S/N 3.71/2.79, 2.18/2.40 and 0.98/1.68 for sectors 83, 56 and 16. No convincing recurring transit is established. A catalogued nearby pulsator 16.32 arcsec away is an additional contamination concern, not a proven assignment.

### TIC 329083521: one strong local fluctuation is not periodic support

The sector-83 signal near 9.37 d is weak under local raw baselines. Independent sectors 77, 84 and 76 were retrieved, including corrected `v02` products. An independently optimized sector-77 fold initially yields local S/N 12.9, but this comes from only one bracketed dip beside large positive excursions. At a common epoch across all four observations, corrected/raw local S/N falls to 0.45/0.19, 2.26/2.34, 0.99/1.04 and 0.59/0.74. The source is the already catalogued emission-line star G 262-38. Activity is a plausible explanation of the isolated excursions; no planetary recurrence is established.

### TIC 182517223: real five-sector eclipses, already reported SB1

The producer-error column is nonfinite in sector 95, although 7,795 quality-zero cadences have valid positive flux. The separate empirical-error search recovered a strong 2.20768-day signal. Follow-up explicitly labels that alternative uncertainty assumption, retains the original producer-error array, and stores its reports in `sector95_empirical_...` folders. It does not replace the frozen producer-error results.

The signal is real: independent QLP sectors 94, 68, 67 and 28 all support a common period approximately 2.2076282 d. Local raw primary S/N ranges from 9.43 to 65.77 in those observations and 60.98 in sector 95 at the common ephemeris. In sector 95, the locally refined raw primary depth is 7,636 ± 103 ppm, while the half-phase secondary is 938 ± 109 ppm (conditional S/N 8.65). The small/large aperture depths are 5,113/8,815 ppm. A secondary/primary-ratio warning threshold of 0.15 missed this significant secondary because the ratio is approximately 0.12; the quantitative measurements were reviewed rather than relying on that flag.

Fresh recognition is decisive: exact Gaia DR3 4684552891681165440 has a published NSS SB1 solution with period 2.2076461525 ± 0.0000230023 d, velocity semi-amplitude 63.9837 ± 0.5519 km/s, eccentricity 0.01358 ± 0.00905, flags 0 and goodness of fit 1.032. The TESS period matches within about 3.25 seconds. This is a known spectroscopic-binary recovery, not a new planet or newly recognized star. Pixel work on this case is lower priority than validating genuinely unrecognized weak leads.

### Other tested leads

The strongest short-period/deep sector-83/86 signals typically exhibit multiple secondary troughs, large orbital variability, aperture dependence or centroid changes. TIC 366489851 has a nearby known RRab star 21.34 arcsec away with published period 0.492685 d, close to half its approximately 0.986-day BLS peak. That is a plausible variable-source alias, not pixel-confirmed attribution.

Longer-period selected peaks, including TIC 364119536, 470310305, 243336717 and 262001813, were tested rather than discarded because they have only a few events. Their local raw support is weak (approximately 1.61, 0.01, 3.35 and −0.26, respectively). Other modest-depth tests and every measured event remain in the machine-readable reports.

## The TIC 4206066 recovery control and a selection limitation

The known user-supplied candidate is explicitly marked as a prior report. Applying this follow-up to actual QLP sector-98 observations recovers 3.182880 d with corrected/raw local S/N 6.91/6.90. It also raises a raw half-period-variability warning. This is a direct example of why an orbital-variability flag must not automatically remove a weak candidate or its host star from multi-sector searches. The control is not a validated planet, so its recovery is a sensitivity check, not a complete false-positive calibration.

The warning thresholds and a fixed box model are imperfect. Higher-frequency stellar variation can evade a two-harmonic check; active stars can have real transits; a selected local baseline can make a fluctuation appear stronger; and independently chosen sector epochs can look persuasive without one common ephemeris. Joint searching, injection/recovery, source localization and recognition remain necessary.

## Preservation and interpretation

Raw QLP and TESSCut files live outside ordinary Git tracking under `/workspace/tess-research-data/multisector/followup/`. SHA-256 manifests, public retrieval URLs, source identifiers, failed response hashes, scientific diagnostic arrays and reasonably sized figures are versioned under `candidates/multisector/`. Retention is intentional; no evidence was deleted. The cache's survival outside this environment is not assumed.

No lead in this report is a new planet or a new star. Strong repeatability of an aperture signal can reveal a known contaminant. Nonmatches against finite catalogs do not establish priority. The blind pooled search of all eligible repeated targets proceeds independently of these follow-up warnings.
