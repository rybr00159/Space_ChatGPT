# TOI-1064: separate vetting of revised gap-trimmed peaks

The 8.87447 d sector-94 and 14.85263 d joint sector-104/105 peaks were re-examined independently after the original boundary peaks were rejected. Neither revised peak is currently a credible new planetary candidate. This finding applies to these signals and does not exclude all additional planets in the system.

Every test includes the known b/c transit masks. Raw SAP and PDC local models omit masked cadences before fitting. Both 1 d and 2 d rolling-median reductions, 0.25 d and 0.5 d good-segment edge exclusions, individual event depths, local polynomial baselines, duration response, and timing-tolerance searches are recorded in `toi1064_revised_trimmed_vetting.json`.

| Source period (d) | Flux | Median width (d) | Source-group fixed S/N, edge 0.25 d | Source-group fixed S/N, edge 0.5 d |
| --- | --- | --- | --- | --- |
| 8.874475 | PDCSAP_FLUX | 1 | 4.49 | 3.90 |
| 8.874475 | PDCSAP_FLUX | 2 | 10.08 | 8.76 |
| 8.874475 | SAP_FLUX | 1 | 1.69 | 1.91 |
| 8.874475 | SAP_FLUX | 2 | 5.06 | 3.45 |
| 14.852630 | PDCSAP_FLUX | 1 | 4.48 | 5.24 |
| 14.852630 | PDCSAP_FLUX | 2 | 9.63 | 9.15 |
| 14.852630 | SAP_FLUX | 1 | -1.72 | -2.28 |
| 14.852630 | SAP_FLUX | 2 | -1.74 | -2.70 |

## 8.87447 d signal

The source sector has three nominal supported windows, but detrending changes the PDC ranking statistic from 10.08 (2 d) to 4.49 (1 d). Under 1 d detrending, local PDC depths are approximately 450, 115, and 308 ppm, while SAP gives approximately −38, +186, and −516 ppm. The third window is only 0.36 d before a quality gap. The original formal odd/even discrepancy is 4.13 sigma, with three events and correlated-error limitations. These tests establish instability rather than a robust shallow transit sequence.

At the exact source ephemeris, other seasons give no coherent depth: the 2 d PDC fixed scores are −2.12 in sector 104 and −0.53 in sector 105 with the 0.25 d edge cut. That is a retrospective exact-ephemeris check. The several-hundred-day separation amplifies an unfit period uncertainty, so it does not eliminate every neighboring alias. Source-event ±0.3 d raw baseline scans are saved rather than interpreting fixed predictions as a complete posterior.

## 14.85263 d signal

The first nominal event at BTJD 4193.86913 occurs only 0.27 d after a large quality gap; a 0.5 d exclusion removes its usable transit window. The remaining two sector-105 event windows are away from the edges. Their local 1 d PDC depths are approximately 357 and 349 ppm, compared with SAP depths of +149 and −579 ppm. Sector 105 alone gives PDC fixed scores of 5.08 (1 d) or 8.99 (2 d), while SAP gives −2.81 or −3.65. The joint score thus relies on broad, reduction-dependent variability, with one source event adjacent to a gap.

The strongest reduction disagreement persists with direct fits to raw 2 min flux after known-transit masking: the event at BTJD 4223.574 has SAP depths from −458 to −1095 ppm across all four local linear/quadratic baselines (−3.1 to −8.5 diagonal sigma), while PDC gives +307 to +1300 ppm. A polynomial-baseline change does not restore a shared downward event. At the first event, raw PDC significance is only 0.1–1.6 diagonal sigma across these models. Source PDC red-noise binning factors at the proposed duration are approximately 1.7–2.1 even for the 1 d detrending.

Both revised searches choose the longest original duration grid point, about 5.5 h. An eccentric planet or uncertain host density could give a long transit; duration alone is not a rejection. The raw polynomial baseline and duration-response measurements test whether a compact event survives alternative local trend choices. They provide no justified physical planet fit for these unstable residuals.

## Scientific scope

The earlier rejection of single boundary excursions was not reused as a blanket decision. The revised peaks fail their own SAP/PDC and detrending-consistency checks. If a future independent pixel reduction restores coherent, similarly shaped events, revisit them with a multi-season ephemeris posterior, injection recovery, and complete false-positive diagnostics. No global significance, planetary radius, discovery, or confirmation is claimed here. All raw FITS remain unchanged and retained.

Reproduce with `python src/vet_toi1064_trimmed.py` after retaining the revised screen JSON and FITS products named by SHA256 in the evidence file.
