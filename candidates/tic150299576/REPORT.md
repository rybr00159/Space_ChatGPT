# TIC 150299576: sector 13 weak long-period alignment

The two central sector 13 events are investigated as a weak lead. They do not currently support a credible recurring planet. Original and independently normalized scientific curves, exact FITS hashes, sector tests and conditional null samples are retained here.

Original screen: P = 11.26922899 d, epoch BTJD = 1663.907703, box duration = 5.51 h. The maximum duration grid point was selected, making baseline choice especially relevant. Both SAP and PDC were analyzed, with 1 d and 2 d rolling-median detrending and independent local linear/quadratic baseline fits to raw 2 min data.

| Sector | PDC fixed S/N, 1 d | PDC fixed S/N, 2 d | SAP fixed S/N, 1 d | SAP fixed S/N, 2 d |
| --- | --- | --- | --- | --- |
| 10 | 0.82 | 0.95 | 0.48 | -0.61 |
| 11 | -1.00 | -0.69 | -0.57 | -0.41 |
| 12 | 1.14 | 1.38 | 1.53 | 2.88 |
| 13 | 6.67 | 8.81 | 2.52 | 2.77 |

The joint 10.5–12 d search uses 1201 periods and five durations. Its best PDC scores are 5.27 at 11.70375 d (width 1 d), 6.28 at 11.29875 d (width 2 d). They shift away from the sector 13 winning period. This does not establish a coherent extension of the original two events.

An empirical moving-block screen used 64 realizations (0.3 d and 0.75 d block sizes) and maxima over the regional grid and both PDC detrenders. 30 null maxima equaled/exceeded the observed joint maximum; the plus-one empirical tail is 0.477. It omits the earlier many-star/full-period selection and is not a global false-alarm probability. The small ensemble does not support precise tail claims.

Independent local raw SAP baselines give only 52–332 ppm for the first source event (0.2–1.9 diagonal sigma), compared with 568–1366 ppm in PDC. The second event has 544–655 ppm in SAP versus 1012–1106 ppm in PDC. These values are conditional on the fixed box and baseline window; their discrepancy argues against promoting the PDC interpretation. Sector 13 PDC red-noise binning factors at the proposed duration are 1.31 and 1.57 for the 1 d and 2 d detrenders.

Each original ephemeris prediction, local raw-baseline depth/error, and red-noise binning factor is in `vetting.json`. The observations in sectors 10–12 were already available for analysis, so this is retrospective robustness. Two-event period uncertainty, period aliases and transit timing variations have not been eliminated by a formal posterior; fixed-time null results alone cannot exclude every nearby ephemeris.

The inferred circular central-transit duration for a roughly 0.62 solar-radius/0.60 solar-mass TIC dwarf at 11.27 d is about 3 h. The selected 5.5 h box is long compared with that scale, but eccentricity and stellar-parameter uncertainty can lengthen it. No radius or physical planet parameters are claimed for a reduction-dependent, nonrepeating lead. No confirmation or formal validation is inferred.

[Diagnostic figure](diagnostics.png). All normalized 10 min scientific curves are retained as CSV.GZ; raw FITS retrieval products and SHA256 values remain in `vetting.json`. Reproduce with `python src/vet_tic150299576.py --null-draws 64`. Investigate pixels or additional sectors only if a transit-preserving independent reduction produces consistent events at a jointly constrained ephemeris.
