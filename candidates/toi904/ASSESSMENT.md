# TOI-904 / TIC 261257684: additional-period residual fails screening

A weak 13.46-day residual appears after known-planet masking in sectors 93/94. Its depth depends strongly on detrending width and local baseline curvature. The frozen-ephemeris comparison in sector 88 provides no significant additional support. These observations do not justify a credible new exoplanet candidate or physical transit fit.

Selection: P=13.45728090 d, T0=3831.18733888 BTJD, duration=5.508 h, formal BLS S/N=7.93. The selected duration sits at the longest searched box. The alternate 1-day detrending prefers P=13.20616 d, demonstrating a broad/inconsistent period preference. Sectors 93/94 supplied the selected ephemeris; sector 88 is observationally separate.

## Independent local baseline measurements

The same selected period/epoch/duration are frozen for all comparisons. Each event has a separate local linear or quadratic baseline using points 1–3 durations from center. The candidate box uses points within half a duration. Known TOI-904 b/c transits and 0.25-day neighborhoods around raw-data boundaries/gaps are excluded. Formal uncertainties are diagonal-noise estimates and are not false-alarm probabilities.

|Sectors|Reduction|Width (d)|Baseline|Depth (ppm)|Formal error (ppm)|Formal S/N|Events|
|---|---|---|---|---|---|---|---|
|[93, 94]|PDCSAP_FLUX|1.0|linear|33|78|0.42|4|
|[93, 94]|PDCSAP_FLUX|1.0|quadratic|27|99|0.27|4|
|[88]|PDCSAP_FLUX|1.0|linear|187|144|1.30|2|
|[88]|PDCSAP_FLUX|1.0|quadratic|290|181|1.60|2|
|[93, 94]|PDCSAP_FLUX|2.0|linear|285|79|3.63|4|
|[93, 94]|PDCSAP_FLUX|2.0|quadratic|117|99|1.18|4|
|[88]|PDCSAP_FLUX|2.0|linear|222|145|1.53|2|
|[88]|PDCSAP_FLUX|2.0|quadratic|216|182|1.19|2|
|[93, 94]|SAP_FLUX|1.0|linear|18|79|0.23|4|
|[93, 94]|SAP_FLUX|1.0|quadratic|32|100|0.32|4|
|[88]|SAP_FLUX|1.0|linear|257|143|1.80|2|
|[88]|SAP_FLUX|1.0|quadratic|360|180|2.00|2|
|[93, 94]|SAP_FLUX|2.0|linear|161|82|1.96|4|
|[93, 94]|SAP_FLUX|2.0|quadratic|88|103|0.85|4|
|[88]|SAP_FLUX|2.0|linear|315|144|2.18|2|
|[88]|SAP_FLUX|2.0|quadratic|340|182|1.87|2|
|[93, 94]|Raw SAP|—|linear|205|75|2.74|4|
|[93, 94]|Raw SAP|—|quadratic|48|94|0.51|4|
|[88]|Raw SAP|—|linear|175|139|1.26|2|
|[88]|Raw SAP|—|quadratic|281|175|1.60|2|

## Duration consistency

The saved confirmed-planet catalog gives Mstar=0.557 Msun, Rstar=0.527 Rsun, density=3.81 solar. At the selected period, the maximum central, circular, small-planet transit duration is approximately 2.77 hours. The BLS box duration is 1.99 times this value. A central transit near apastron would need eccentricity approximately 0.60 under these assumptions; other possibilities include incorrect stellar density or a biased box duration. This is a consistency warning, not a measured eccentricity or proof against a planet. A larger impact parameter shortens the circular duration.

## Assessment and limitations

- The local linear depth falls from roughly 285 ppm at 2-day detrending to 33 ppm at 1-day detrending; a quadratic baseline reduces the 2-day value further. This sensitivity is consistent with slow stellar/instrumental curvature entering the long box statistic. No repeated isolated transit-shaped events are established.
- The 13.46-day trial is distinct from catalog periods 10.8772 and 83.9997 days, but that alone establishes neither novelty nor astrophysical origin. Current known-transit masks are used, with ephemeris uncertainties not fully propagated.
- Sector 88 tests the frozen nominal ephemeris without a new search. Its weak positive measurements do not independently validate the signal. The selected broad peak has no rigorous timing/period posterior, so failure at the nominal phase would not alone exclude a true planet over the extrapolation interval.
- SAP and PDCSAP are reductions of the same photons. Their comparison tests processing sensitivity, not independent astrophysical confirmation. Raw SAP baseline fits use empirical adjacent-point noise floors, while correlated stellar/instrumental variability remains incompletely modeled.
- BLS explored 110,061 period values in the selected two-sector screen, plus epochs/durations and multiple stars, sectors and detrenders. An attractive formal maximum is insufficient. Local off-phase controls and binned noise ratios in the JSON are descriptive and cannot calibrate this global trial factor.
- The source is not localized using pixels; background binaries and other blends remain untested. No exhaustive additional-candidate literature novelty claim is made. A physical fit, planet radius, impact parameter or false-positive probability would overinterpret the weak evidence.
- This screening rejection is not evidence that TOI-904 lacks further planets. Longer baseline data and variability models with injection-tested short-transit preservation could revisit the system without presuming the selected period.

Reproduce with `python src/vet_toi904.py`; inputs and product checksums are preserved in `data/observations.json`. The script does not download or rerun the period search.
