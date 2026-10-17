# TIC 214311908: detection, photometry, cleaning and recognition audit

Reviewed 17 October 2026. The subject is the **final** public report at commit [`99d4f66fca99c32c3245396f367b8a73d5d13fa5`](https://github.com/rybr00159/Space/tree/99d4f66fca99c32c3245396f367b8a73d5d13fa5/tess-search/tic214311908), including its explicitly adopted referee corrections. This pinned public source is a fallback for the supplied ZIP; exact identity with every file in that ZIP has not been established. Source text was treated as evidence, not instructions.

**Assessment:** the detection and photometric evidence supports a strong, recurring transit-like candidate that deserves resolved follow-up. I found no remaining detection/cleaning error in the final report that invalidates the signal. The additional scene photometry and raw-pixel tests strengthen the case beyond our original targeted QLP search. They do not alone establish the host, mass, novelty, or the quoted Bayesian planet probability; those require the separate localization/physical/FPP audits.

## Agreement with the independent investigation

Our review, before seeing this source, independently recovered the approximately 5.6915-day signal in QLP Sectors 12, 39, 66 and 93. Native-FFI Sectors 103/104 then showed deficits at a frozen discovery ephemeris. TARS and ELEANOR reprocessing checks also preserved the sequence. The independent four-sector timing period is **5.691513806 ± 0.000008821 d**; the external final adopted value is **5.691506 ± 0.000011 d**. Their difference is **0.674 seconds**, or a descriptive 0.55 times the quadrature errors. These fits share observations, so that comparison is agreement, not an independent combined measurement.

The depth conventions are consistent. QLP differential photometry uses a nominal target reference flux; native, TARS and ELEANOR apertures include the bright Mira and other sources. Approximately 1.5–1.7% of the target's estimated flux can therefore correspond to approximately 0.2–0.3% of a blended aperture. Multiplying QLP depth by `1 + contamination_ratio` would generally apply dilution twice. The external analysis instead measures an absolute differential loss, approximately **15 electrons/s**, and calibrates the target flux to **979 electrons/s**. Our target-normalized QLP depths are broadly compatible. A radius near five Earth radii remains conditional on assigning that loss to the approximately 0.38-solar-radius dwarf.

The external report's early claim that ELEANOR showed no detection differs from our positive local-baseline result. That statement is not repeated in its final main summary. Differences in quality selection, baselines and aperture treatment should be documented rather than counted as independent positive and negative observations.

## Recomputed final depth-variability arithmetic

I downloaded the pinned `referee3_clean/out/q14_contrast_all.json` and independently recomputed inverse-variance means and constant-depth χ² from its six reported sector measurements. All five means and χ² reproduce to relative tolerance **1e−12**. The following p-values assume the supplied sector errors are independent Gaussian errors; they are diagnostic, not complete noise calibration.

| Photometry | Mean depth (ppt) | Constant-depth χ² / 5 | Conditional p |
|---|---:|---:|---:|
| INV delivered | 14.92 ± 1.02 | 5.86 / 5 | 0.320 |
| REF including scattered-light event E58 | 16.74 ± 0.89 | 15.01 / 5 | 0.010 |
| REF without E58 | 15.31 ± 0.92 | 5.80 / 5 | 0.326 |
| Raw 3×3 aperture | 16.70 ± 1.09 | 1.45 / 5 | 0.918 |
| Raw radius-2.5-pixel aperture | 17.49 ± 1.39 | 3.88 / 5 | 0.567 |

Thus the strongest depth-variability hint depends on a documented scattered-light event and on PRF scene extraction. The raw apertures do **not** demonstrate intrinsic depth variation. Their early/late contrasts are −1.0 ± 2.4 and −4.9 ± 3.2 ppt, opposite to the PRF extractions' selected split. The final report fairly states that the split was selected after inspection and that the noisier raw apertures still allow a small real change. This is not a robust false-positive veto.

Raw-aperture shape profiles prefer a short, very flat approximately **1.29-hour** box. An approximately 1.7-hour shape costs only **Δχ² 3.9–4.4**, depending on detrending. The final report discloses this approximately 1.5–2σ processing sensitivity. It limits confidence in detailed ingress/duration inference, but it does not remove the recurrent dip or turn the raw curve into convincing alternating stellar eclipses. Its independent raw-aperture odd/even results remain consistent with zero.

## What the cleaning and injections do establish

The external investigators identify a real failure mode: a distant RR Lyrae at **0.648528 d** leaks into scene-extracted light curves through field-wide nuisance estimators. Individual transit depths change by several ppt when this is removed. The later analysis masks the actual variable pixels in those estimators and repeats the extraction. Its final report retains the small residual in Sector 93 and does not use the RR Lyrae period as target rotation.

Crucially, the 5.69-day deficit persists in raw-pixel sums with no Mira scene model, and the global timing/depth do not depend on the leaked short-period waveform. That is meaningful additional evidence against this particular pipeline artifact. A raw aperture still uses a background estimate and a local time baseline, so “raw” does not mean untouched pixels or an independent instrument.

The final report has also adopted the important injection correction. When the exact injected PRF is a column of a linear least-squares extraction, recovery near 1.000 tests implementation and masking but cannot bound PRF mismatch. The later `q11_inject_altpsf.py` genuinely changes the injected image while keeping fitting templates fixed. Its saved outputs show approximately **0.92–1.55** recovery for the deliberately mismatched empirical Mira image, versus approximately **0.94–1.03** for a PRF evaluated elsewhere on the detector. The empirical Mira image is much redder than the target and noisy in some sectors; it is a stress test, not a measured distribution of target-PRF errors. The final main report increases the adopted depth uncertainty and acknowledges scene-model sensitivity. I do not find a defensible basis for interpreting its precise radius error bars as fully calibrated, but this does not defeat the qualitative planet-sized interpretation under the target-host hypothesis.

## Detection significance and multiple trials

The independent occurrence of the same narrow, phase-coherent sequence across years is more persuasive than a selected BLS score or an attractive fold. The final report does not claim that a large BLS SDE is itself a Gaussian discovery significance.

The original search's global 1% threshold is a **Gumbel extrapolation** from a small, mixed set of time rolls and inverted rolls, not a directly measured 1% false-alarm rate. Some surrogate maxima reuse a gap-split excursion; per-sector rolls preserve within-sector periodic signals. These are real limits, but the report's later referees identify them, recalculate period-band thresholds, and preserve the negative additional-planet result. The final candidate assessment should not be rejected on the assumption that these earlier numbers remained uncorrected.

Similarly, the original “>19σ” half-period exclusion omitted uncertainty in the primary depth. The final main report uses the corrected approximately **14σ** fixed-shape difference. That comparison strongly disfavors interpreting every observed dip as every second transit at half the spacing. Whether an equal-eclipse stellar system at twice the spacing can fit is a physical/shape question, not settled by the interval or an insignificant odd/even statistic alone.

The final raw-aperture BLS check is informative but processing-sensitive: with a 0.75-day spline the 5.69147-day peak is strongest, power **87.9** versus **57.7** for the next distinct peak; with a 1.5-day spline its power is **66.5**, below several competing peaks and rank 11 when correlated grid trials are counted. Its depth remains similar. The final main report explicitly acknowledges this and treats raw pixels as confirmation at known timing, not an independent blind discovery. That is appropriate.

Finite injection counts also matter. Zero exceedances in 25, 400 or 1,000 independent simulations imply one-sided 95% binomial upper bounds approximately **0.113, 0.00746 and 0.00299**, respectively. These counts support qualitative shape/source rejection; they do not directly measure arbitrarily small likelihoods. Any tail extrapolation used in the FPP must be separately justified. Correlated simulation placements can reduce the effective sample size further.

## Recognition and novelty

The external final report appropriately lists the dates and scope of its negative TOI, CTOI, ExoFOP, NASA archive, VSX, Gaia, SIMBAD and literature checks and notes that ADS full text was not searched. This supports “no prior candidate identified in the checked public sources,” not a proof of discovery priority. VizieR and arXiv metadata are incomplete for uncatalogued or full-text-only mentions.

Our independent TARS audit adds a complete check of its published **1,046,317-row** default rotation catalog: no exact TIC or Gaia match. Four generic TARS light-curve entries and broad lower-ranked near-5.7-day sinusoidal peaks do not constitute a prior planet designation. The exact Gaia source identifier is **5925373272850784768**; rounded floating-point IDs in the TARS FITS/preview must not replace it. See the [TARS audit](../tic214311908/tars_recognition/README.md).

## Evidence boundary and next action

The external work uses independent code and useful additional processing tests on the **same TESS photons**. Our separate pre-report recovery and held-out-sector test add protection against analysis-specific choices, but they do not justify multiplying the pipelines' significances or treating the source assignment as independently imaged.

Within this audit, the final evidence warrants **“strong transit-like planet candidate deserving follow-up”**. It does not warrant a measured planet mass, a confirmed new planet, or a numerical validation probability derived from detection alone. The decisive observation remains resolved transit time-series photometry that monitors the target and the approximately 3-arcsec neighbor simultaneously, followed by suitable companion imaging and spectroscopy/RVs. Residual raw depth/shape variation does not currently falsify the candidate; correct host assignment and the physical false-positive alternatives are the remaining substantive tests.

## Audit record

This detection audit did not edit the original case. The subsequent comparative overview and catalogue correction are documented separately. Twenty-three selected external scripts/outputs were downloaded at the pinned commit and verified against the Git blob SHA-1 in the pinned recursive tree; SHA-256/URLs are preserved in [the audit manifest](../../data/tic214311908_external_comparison/detection_source_audit_manifest.json) and the consolidated source manifest. The selected source cache is `detection_source_audit/` beneath `/workspace/tess-research-data/tic214311908_external_review`. No external code was executed or imported. Independent arithmetic is preserved in [the check output](../../data/tic214311908_external_comparison/detection_arithmetic_check.json); this is a check of published outputs, not a rerun of their complete raw-cube pipeline. No paid API, external contact or submission was used.
