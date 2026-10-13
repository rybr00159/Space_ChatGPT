# TOI-4342: follow-up of the edge-trimmed two-event alignment

After boundary rejection, sector 102 gives P=6.67121422 d, epoch BTJD 4132.986646, duration 0.180 d. The observed events are separated by 13.342428 d. The intervening 6.671d prediction lies in the large quality gap, so 13.342d is equally compatible with these two event times. This specific signal remains unpromoted.

| Sector | Flux | Width (d) | Fixed 13.342d diagonal S/N |
| --- | --- | --- | --- |
| 101 | PDCSAP_FLUX | 1 | 0.38 |
| 101 | PDCSAP_FLUX | 2 | -0.66 |
| 101 | SAP_FLUX | 1 | 0.95 |
| 101 | SAP_FLUX | 2 | -0.63 |
| 102 | PDCSAP_FLUX | 1 | 8.53 |
| 102 | PDCSAP_FLUX | 2 | 9.93 |
| 102 | SAP_FLUX | 1 | -1.32 |
| 102 | SAP_FLUX | 2 | -1.20 |
| 103 | PDCSAP_FLUX | 1 | -1.02 |
| 103 | PDCSAP_FLUX | 2 | -0.28 |
| 103 | SAP_FLUX | 1 | -1.40 |
| 103 | SAP_FLUX | 2 | -0.15 |

In the source sector, the 1d-detrended PDC local depths are approximately 1055 and 1433 ppm; corresponding SAP depths are approximately −228 and −401 ppm. The same ephemeris has no significant PDC repetition at fixed times in neighboring sectors. Exact event measurements and a ±0.3d timing-sensitivity search are saved in `toi4342_trimmed_alias_check.json`. Extra timing searches introduce extra trials; their maxima are not discovery significance.

[Event-shape comparison](toi4342_trimmed_alias_shapes.png) includes both original raw flux columns, independently normalized, and both independently flattened columns. Opposite signs can be caused by additive instrumental backgrounds, so this is evidence of reduction dependence rather than a theorem excluding all transits. Pixel photometry and longer-baseline alias fitting would be required if future observations revive the signal. No current evidence justifies calling it a credible new candidate.

The two-event period uncertainty was not fitted. Fixed-time null results alone do not exclude every nearby period or timing variation. No unpublished/discovered status, formal false-positive probability, or planetary radius is assigned.
