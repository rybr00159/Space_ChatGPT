# Pipeline validation on actual TESS data

These tests assess basic detection and filtering behavior; they do not establish novelty or calibrated false-alarm probabilities.

## Catalog ephemeris controls

All known transiting L 98-59 planets are measured at fixed catalog ephemerides, with other known planets excluded. Each event uses a separate weighted linear local baseline. Errors assume diagonal noise; event scatter and off-transit phase controls expose limitations.

|Sector|Filter (d)|Planet|Depth (ppm)|Formal S/N|Event scatter S/N|Events|Phase-control S/N std|
|---|---|---|---|---|---|---|---|
|2|1.0|L 98-59 b|484|9.4|8.9|11|0.87|
|2|1.0|L 98-59 c|1524|25.3|25.3|6|1.02|
|2|1.0|L 98-59 d|1387|15.5|12.3|4|1.01|
|2|2.0|L 98-59 b|495|9.5|9.4|11|1.16|
|2|2.0|L 98-59 c|1519|24.8|24.8|6|0.87|
|2|2.0|L 98-59 d|1395|15.3|12.8|4|1.00|
|5|1.0|L 98-59 b|636|13.0|10.7|12|1.12|
|5|1.0|L 98-59 c|1495|23.0|23.0|6|1.15|
|5|1.0|L 98-59 d|1148|3.7|3.7|2|1.03|
|5|2.0|L 98-59 b|637|12.9|10.6|12|1.27|
|5|2.0|L 98-59 c|1499|22.9|22.9|6|1.08|
|5|2.0|L 98-59 d|1124|3.6|3.6|2|0.94|
|8|1.0|L 98-59 b|505|8.5|8.5|9|1.02|
|8|1.0|L 98-59 c|1543|21.6|21.6|5|0.76|
|8|1.0|L 98-59 d|1302|11.7|11.7|3|1.00|
|8|2.0|L 98-59 b|502|8.4|8.4|9|1.07|
|8|2.0|L 98-59 c|1534|21.3|21.3|5|1.09|
|8|2.0|L 98-59 d|1301|11.6|11.6|3|1.06|

## Pre-detrending injections and unchanged-data control

Box-shaped signals are injected multiplicatively into quality-selected raw PDCSAP flux before median filtering and 10-minute binning. Errors are preserved. A fixed epoch is chosen using coverage only, then the full 0.7–13-day production BLS search runs without the true ephemeris. This favorable epoch has 100% of injected raw cadences outside known-planet masks; it is not a random-phase completeness test. Matching requires period within 1% and phase within half the injected duration. Zero-depth rows are unchanged real residual data, not simulations from a validated noise null.

|Filter (d)|Injected (ppm)|Top-five rank|Measured added depth (ppm)|Depth retained|Best blind P (d)|Best formal S/N|
|---|---|---|---|---|---|---|
|1.0|0|—|0|—|1.05026|5.67|
|1.0|500|1|487|0.974|4.73424|8.24|
|1.0|1000|1|973|0.973|4.73424|16.22|
|1.0|2000|1|1948|0.974|4.73322|33.06|
|2.0|0|—|0|—|1.05026|5.47|
|2.0|500|1|490|0.980|4.73424|8.05|
|2.0|1000|1|981|0.981|4.73424|16.22|
|2.0|2000|1|1958|0.979|4.73322|33.23|

## Interpretation and limitations

- The shallow known planet is a more meaningful end-to-end sanity check than only retrieving the deepest peak. Highest peaks include harmonics; top-three output alone cannot recover every known planet.
- Injections use one star, one sector, one period, one epoch, and ideal boxes. They cannot quantify pilot-wide completeness, detectability of longer periods, ingress/egress effects, or global false-positive frequency. There are no repeated noise realizations.
- Broad known-planet masks preferentially remove certain trial phases and reduce observing duty cycle. They should be expanded using ephemeris uncertainties for recent sectors, and verified against predicted events. Flattening occurs before masking, so known dips can affect neighboring trends.
- Median-filter width, incomplete events near data gaps, and 10-minute binning can alter depths and durations. Event scatter and 60/90/120-minute noise ratios in results/validation.json are descriptive diagnostics, not calibrated uncertainty models.
- BLS S/N and robust SDE are ranking scores, not Gaussian significance or discovery probabilities. Multiple stars, periods, epochs, durations and detrending choices create a substantial trial factor. Off-phase tests here cannot calibrate the maximum over that search.
- The shortest searched box is 1.08 hours; genuine shorter events such as L 98-59 d (catalog 0.84 h) are diluted. Single/double transits and periods beyond half-sector baseline are outside the blind periodic screen.
- A positive in-transit count is not evidence of repeated events: per-event depths, independent sectors, raw/SAP comparisons, spacecraft artifacts, source localization and novelty checks remain essential.

## Empirical full-grid null appendix

See [null_calibration.md](null_calibration.md) and `results/null_calibration.json` for half-day block permutation comparisons. These illustrative maxima calibrate neither the full multi-target search nor discovery probabilities.
