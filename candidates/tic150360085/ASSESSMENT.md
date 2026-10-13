# TIC 150360085: weak residual rejected as a credible candidate

A joint four-sector BLS screen selected a roughly 13.123-day signal. Independent local baseline measurements weaken it, and early predicted events lack significant individual depth. Their uncertainties are large enough that a shallow planet cannot be ruled out. This remains an exploratory residual, with no evidence sufficient to designate a credible new exoplanet candidate.

Frozen selection ephemeris: P=13.12303747 d, T0=3831.21574333 BTJD, duration=5.508 h. The 5.5-hour duration lies near the longest searched duration. Both detrenders select this period, but share observations and are correlated analyses. All sectors 93–96 contributed to selection; none are independent confirmation.

## Independent local-baseline measurements

Each predicted event uses a separate linear baseline from points 1–3 durations away; a box depth measures the central half-duration interval. The table uses the same frozen ephemeris for every reduction. Quadratic baselines are a diagnostic for local curvature. Formal errors assume independent points; event scatter is an additional warning, not a calibrated red-noise model.

|Reduction|Baseline|Depth (ppm)|Formal error (ppm)|Formal S/N|Event-scatter S/N|Events|
|---|---|---|---|---|---|---|
|PDCSAP_FLUX width 1.0 d|linear|262|72|3.62|2.42|8|
|PDCSAP_FLUX width 1.0 d|quadratic|184|90|2.04|2.04|8|
|PDCSAP_FLUX width 2.0 d|linear|275|73|3.77|2.94|8|
|PDCSAP_FLUX width 2.0 d|quadratic|160|91|1.76|1.76|8|
|SAP_FLUX width 1.0 d|linear|247|77|3.21|2.27|8|
|SAP_FLUX width 1.0 d|quadratic|165|96|1.72|1.72|8|
|SAP_FLUX width 2.0 d|linear|268|78|3.43|3.04|8|
|SAP_FLUX width 2.0 d|quadratic|140|98|1.44|1.44|8|
|Raw SAP, no median detrending|linear|241|72|3.32|2.77|8|
|Raw SAP, no median detrending|quadratic|118|90|1.31|1.31|8|

## Per-sector support

|Sector|Reduction|Depth (ppm)|Formal error (ppm)|Formal S/N|Events|
|---|---|---|---|---|---|
|93|PDCSAP_FLUX width 1 d|-37|200|-0.19|2|
|94|PDCSAP_FLUX width 1 d|374|209|1.79|2|
|95|PDCSAP_FLUX width 1 d|157|153|1.02|2|
|96|PDCSAP_FLUX width 1 d|355|100|3.56|2|
|93|SAP_FLUX width 1 d|-27|203|-0.13|2|
|94|SAP_FLUX width 1 d|278|224|1.24|2|
|95|SAP_FLUX width 1 d|180|162|1.11|2|
|96|SAP_FLUX width 1 d|346|107|3.23|2|

## Predicted event consistency

|Cycle|Center (BTJD)|Depth (ppm)|Formal error (ppm)|Formal S/N|
|---|---|---|---|---|
|0|3831.21574|10|286|0.03|
|1|3844.33878|-83|281|-0.29|
|2|3857.46182|139|304|0.46|
|3|3870.58486|583|287|2.03|
|4|3883.70789|960|292|3.29|
|5|3896.83093|-149|180|-0.83|
|6|3909.95397|322|132|2.45|
|7|3923.07701|399|152|2.62|

## Interpretation and limitations

- Counted events mean sufficient sampling, not significant repeated transits. Local depths include zero/negative events. A constant-depth interpretation has excess event scatter; sector-by-sector support varies strongly.
- Alternative detrending widths and local baseline curvature change the result. SAP comparison shares the same photons and cannot independently verify astrophysical origin. Raw SAP local-linear and quadratic measurements are included in the JSON; their diagonal errors can understate stellar/instrumental covariance.
- Noise ratios on duration-scale bins and off-phase fixed-ephemeris control statistics are descriptive. They are neither a period-search false-alarm probability nor a proven error inflation rule. The ephemeris was selected among hundreds of thousands of period trials and multiple stars/sectors/reductions.
- Odd/even or secondary-eclipse nondetections have little diagnostic power at this signal strength. Favorable aperture crowding fractions alone cannot rule out background blends, nor establish a transit on the target.
- A physical transit/trapezoid fit is not justified: consistent transit-shaped events have not survived screening. No radius, impact parameter, or planetary interpretation is claimed. Pixel-level source localization and exhaustive novelty checks would be required if stronger independent evidence emerges.
- No archived sectors 87–90 were used in this evaluation. A fixed-ephemeris comparison with independently retrieved older data could falsify or revive the hypothesis, provided ephemeris error from the selected broad peak is propagated. Lack of a nominal predicted event alone would be inconclusive over a long extrapolation.
- This negative screening judgment does not prove the star lacks planets; it means the selected residual lacks sufficient evidence to survive this pilot's vetting.

All inputs, observation identifiers, processing versions, and retrieval checksums are linked through `data/observations.json`. `src/vet_tic150360085.py` regenerates this assessment and diagnostic plots without a new period search.
