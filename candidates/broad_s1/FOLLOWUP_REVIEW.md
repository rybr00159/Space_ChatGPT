# Full strict-queue review and independent follow-up

The completed broad survey searched 159,505 sector-1 QLP light curves. This dossier phase examined every **389 distinct signal rows on 357 stars** in the frozen current-stellar-property uncapped queue. The main batch reviewed 357 distinct stars; the separate alternate-peak batch reviewed the other 32 signal rows. Both completed with zero errors. `all_signal_evidence.csv` joins their actual refined measurements to current stellar/catalog provenance. `all_signal_review_summary.json` verifies exact queue membership and the SHA-256 of all 357 original FITS files.

These are stringent, prioritized periodic signals, not a completeness census of the survey. The current selection excludes many already known or physically inappropriate contexts and does not vet every single event or low-scoring periodic feature. No credible new planet candidate is accepted from this phase. Eight first-sector signals passed the automated warning thresholds; later recurrence, source, catalog and physical-model checks change their interpretation.

The numeric method remains frozen at `raw-aperture-all-orbital-phase-raw-variability-v3`, source SHA-256 `23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333`. Queue/checkpoint/ranking histories and original raw observations are retained. Supplementary scripts below do not mutate that method. The erroneous CD-79 151 contaminant sentence inherited by unrelated multi-sector reports was removed with the text-only renderer; it applies only to TIC 370247818, and no scientific arrays changed.

## Strong recurring signals

| TIC | Follow-up observations | Important evidence | Current assessment |
|---|---:|---|---|
| 238202207 | 17 QLP sectors, 140 conditional event times | Timing period 2.321825237 ± 0.000001923 d; integrated trapezoid depth 17011 ± 310 ppm, duration 1.989 h. Ordinary half-phase secondary 33 ± 187 ppm in raw SAP. | Real recurring eclipse; planet and binary remain physically degenerate. Alternating depths and severe Gaia unresolved-source warnings strengthen the binary alternative. Source/physical checks are required. |
| 269994642 | Sector-1 QLP plus independent native-pixel work reported separately | Original QLP period 1.74649452 d, raw local S/N 18.7, ten events; integrated depth 24126 ± 2716 ppm, duration 1.727 h. | Grazing transit/equal-light binary degeneracy; new alternating-depth evidence motivates an unequal binary at twice the recurrence period. Later-sector QLP products were absent, which is not a photometric nondetection. |
| 300814564 | 23 usable QLP sectors, 224 event times | Period 1.860044776 ± 0.000000356 d; raw half-phase secondary **2179.7 ± 62.1 ppm**, conditional S/N 35.1. | Stellar eclipse strongly favored. Sector 37 was retained but excluded from this analysis because no quality-zero rows had valid finite positive flux and error. |
| 273879358 | Four QLP sectors | Raw additional-phase secondary at phase 0.6937: **2006 ± 220 ppm**, S/N 9.1; independent Gaia RV variability. | Stellar binary/source alternative strongly favored. |
| 270312735 | Six QLP sectors | Raw additional-phase secondary at phase 0.5126: **5584 ± 619 ppm**, S/N 9.0, also present in larger/smaller apertures. | Stellar binary/source alternative strongly favored. |

The last two illustrate why examining the whole orbital phase is essential: an eccentric eclipsing binary can avoid an exactly half-phase test. Aperture/centroid warnings by themselves were treated as unresolved source evidence, not as proof against a planet. Multi-sector secondaries supply much stronger astrophysical counterevidence.

## Alternating-depth tests

`tic238202207_multisector/parity_phase_qa.py` verifies the original FITS hashes, repeats local quadratic baselines at the fitted timing ephemeris and fits separate mean depths for each sector plus an even-minus-odd term. Exposure-integrated fixed-shape event fits give **1983 ± 399 ppm** in raw SAP, conditional S/N 4.97. Smaller and larger apertures give 1592 ± 385 and 1878 ± 509 ppm. The conservative error is the larger of chi-square-inflated covariance and sector-cluster covariance. Sixteen of 17 sector differences have the same sign. Earlier 30-minute sectors give 1770 ± 626 ppm; later 10-minute sectors give 2129 ± 518 ppm. Leave-one-sector S/N spans 4.34–5.10. Alternating event times are insignificant, 0.000302 ± 0.000585 d. Sector-dependent raw harmonics vary substantially and do not supply a clean ellipsoidal-mass measurement.

`tic269994642_p1.746821/parity_qa.py` uses only the original independent QLP sector-1 source. Raw integrated event fits give **5818 ± 1630 ppm**, conditional S/N 3.57; smaller/larger apertures give 4683 ± 1580 and 5970 ± 2266 ppm. Error includes empirical even/odd event dispersion. This agrees qualitatively with independent native-pixel binary-model tests, but is only one sector and does not establish a calibrated exclusion probability.

These tests concern alternating eclipses, not an ordinary off-phase secondary. At twice a planet's genuine recurrence period, the same planetary transit appears at two phases; that fact alone cannot reject the planet. A statistically supported difference between those two eclipse populations is more informative. Errors remain conditional on the estimated ephemeris, template and baselines. Search trials, correlated noise, source blending and possible time-dependent stellar activity preclude converting these values directly into discovery or exclusion probabilities.

The independent physical-model dossiers `tic238202207_physical/` and `tic269994642_physical/` support the same unequal-binary alternative. Their equal-light binaries fit the grazing profiles as well as dark-companion transits, while unequal-brightness binary fits at twice the recurrence periods improve conditional data chi-square by approximately 27.6 and 33.1. Independent native-pixel event checks for TIC 269994642 agree in sign across sectors 1 and 92; those photons are not independent of their source pixels, so the two analyses must not be multiplied as independent probabilities. Binarity is favored, with mass and exact source still unmeasured.

`tic238202207_multisector/parity_cadence_controls.py` tests three constant-depth circular limb-darkened BATMAN grazing planets at every observed timestamp, quality mask and native exposure. Despite using a different physical profile from the fitted trapezoid, noiseless controls create only 2.3–2.6 ppm artificial alternation. Three seeded white-noise controls give 131, -317 and 367 ppm, all below 0.91 conditional sigma. This rules out the specific native-cadence/template bias for these examples, not correlated stellar/systematic noise or the survey's trials.

## Weak and additional leads

TIC 431481555 is a previously recognized DSCT/GDOR/SXPHE pulsator. Its broad 5.635-day sector-1 feature has weak later-sector support and is nearly seven times the strongest broad variability period of 0.8057 d. Eight exploratory sinusoidal modes trained outside generous primary windows reduce its sector-1 raw local S/N from 9.73 to 3.94. This supports a pulsation/beat explanation without proving that the star cannot also host planets. Script, mode parameters and plots are in `tic431481555_multisector/pulsation_qa.*`.

TIC 139148478 has weak sector-28 support. Its light curve shares strong long variability with bright rotational neighbor HD 217027/TIC 139148479: native raw correlations are 0.989 in sector 1 and 0.855 in sector 28. The neighbor does not significantly reproduce the target's independently fitted approximately 0.70-day dips (neighbor S/N 0.12 and 1.43). The comparison therefore identifies shared variability but does not identify the dip source. No planet promotion is justified. See `tic139148478_multisector/neighbor_comparison.*`.

TIC 260820871's alternate approximately 0.673095-day signal warrants independent-sector testing: 28 positive raw events, S/N 13.25 and approximately 437 ppm depth. Its current radius is 2.408 solar radii and temperature 10942 K; initial odd/even difference is 3.2 conditional sigma. The two approximately 1.347-day folded eclipses alone cannot reject the shorter recurrence. TIC 290028406 has three raw events at approximately 6.772758 d, S/N 15.44, and current radius 1.088 solar radii; source/aperture and grazing-shape warnings remain unresolved. These are next-stage leads, not accepted planets.

Three other warning-free first-sector stars are weak and poorly constrained: TIC 410255802 and 79492853 lack current validated stellar radius/mass/temperature and require TIC source-disposition review; TIC 419589172 has a broad 9.6-hour feature and current contamination ratio 8.70. All have only three measured events.

The independent alternate-peak audit also records a strict quality-gap boundary in TIC 140803419: a real additional approximately 2% dip was omitted from the local secondary statistic when a gap barely exceeded three median cadences. Its zero secondary statistic is incomplete coverage, not proof of no eclipse. Related approximately 0.7012-day signals on TIC 139199396 and 152819299, nearly two degrees apart on the same detector, have strongly correlated high-pass flux and favor a shared instrumental/optical origin. See `extra_signals/independent_audit.*` for exact scope and source evidence.

## Reproduction and remaining work

Run all commands from the repository with the installed astronomy Python environment:

```bash
MPLCONFIGDIR=/tmp/tess-mpl /workspace/openrouter-check-venv/bin/python candidates/broad_s1/tic238202207_multisector/parity_phase_qa.py
MPLCONFIGDIR=/tmp/tess-mpl /workspace/openrouter-check-venv/bin/python candidates/broad_s1/tic269994642_p1.746821/parity_qa.py
/workspace/openrouter-check-venv/bin/python candidates/broad_s1/tic300814564_multisector/run_filtered.py
/workspace/openrouter-check-venv/bin/python candidates/broad_s1/render_multisector_notes.py
/workspace/openrouter-check-venv/bin/python candidates/broad_s1/summarize_review.py
```

Complete independent native-pixel source assignment and physical unequal-binary comparisons for TIC 238202207 and 269994642. Then investigate independent coverage for TIC 260820871 and 290028406, and obtain source/stellar clarification for the three weak unresolved stars before interpreting their planet sizes. Further broadening should target additional sectors, revised cadence-aware gap tolerance and noise-calibrated search completeness, rather than treating the present strict-queue rejection rate as the absence of planets in all searched stars.
