# TIC260820871: physical morphology and source ambiguity

Observed photometric recurrence P=0.673332598 d has alternating depths. A free integrated trapezoid on the same native photons gives even/odd depths 631.7/395.6 ppm, T14=2.131 h, ingress fraction=0.408. The physical orbital interpretation of alternating eclipses is approximately twice this period. This diagnostic does not measure companion mass.

The nominal current TIC scenario is R=2.40763 Rsun, M=2.85 Msun, Teff=10942 K. R/M20% and temperature1000K widths are explicit conservative assumptions, not catalog errors. The circular a/R at the photometric period is about1.9. Exact native30min/10min exposure durations, raw byte hashes and TDB time references are verified. The small depth and unusually prolonged ingress require testing physical geometry rather than identifying radius with sqrt(depth).

| Primary hypothesis | k | b | T14 h | Conditional R/Rjup | Data χ² | Density/physical prior χ² | Photosphere clear |
|---|---:|---:|---:|---:|---:|---:|---|
| identical_primary_at_photometric_P | 0.60000 | 1.5850 | 2.171 | 14.06 | 2818.6 | 0.36 | True |
| 2P_even_primary_e0_omega90 | 0.03306 | 0.9858 | 2.128 | 0.77 | 1455.8 | 6.01 | True |
| 2P_even_primary_e0.15_omega270 | 0.03347 | 0.9869 | 2.086 | 0.78 | 1455.7 | 4.18 | True |
| 2P_even_primary_e0.3_omega270 | 0.03503 | 0.9904 | 2.078 | 0.82 | 1455.5 | 2.48 | True |
| 2P_even_primary_e0.45_omega270 | 0.03915 | 0.9979 | 2.085 | 0.92 | 1455.2 | 1.12 | True |
| 2P_even_primary_e0.3_omega90 | 0.12000 | 1.0930 | 5.287 | 2.81 | 1454.4 | 8.60 | True |
| 2P_even_primary_free_density | 0.12000 | 1.0941 | 2.226 | 2.81 | 1454.4 | 0.01 | True |

Every-other primary fits use the physical2P orbit. Eccentricity scenarios test apastron/periastron velocities and enforce a(1-e)>Rstar+Rcomp; an arbitrary high-e long-duration fit that intersects the photosphere is inadmissible. The free-density scenario deliberately relaxes target stellar constraints. These are scenario comparisons with diagonal errors, not a mass validation, formal model selection, or rejection of all planetary possibilities. Even-primary χ² values use fewer photons than the all-event empirical comparator and must not be directly subtracted from it.

A separate unknown-source detached binary uses two equal radii, unequal surface brightness, and physical2P orbit, with a/R≥3 for detached Roche geometry. It fits dataχ²=2755.3, compared with empirical all-event χ²=2753.4; brightness ratio=0.626, binary fraction of aperture light=0.00207, b=0.788, a/R=8.22. Source density/temperature are free; target properties are not incorrectly imposed on a faint blended binary. Its illustrative diluted ellipsoidal amplitude is 5.5 ppm (alpha1.5,q1 approximate scaling). A faint binary can therefore have undetectable broad orbital humps. The catalog contamination ratio0.00246 is a prediction, not a source-localization measurement or a guaranteed exclusion of faint deep eclipses.

A thermal/reflection budget cannot substitute for a feasible, localized planetary primary geometry. Optical eclipse brightness of a hot or intrinsically luminous object never establishes mass. Pulsation/rotation, hot-star geometry, blended stellar eclipses, uncertain catalog radius and source attribution remain alternatives. **No credible new exoplanet is accepted from these tests.** Independent pixel localization, actual stellar constraints and source-specific spectroscopy/photometry are needed.

Reproduce: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python src/fit_260820871.py`. Frozen earlier two-lead source hash, original raw identifiers/hashes, fitted parameters, model arrays and diagnostic plots are preserved. No observations were deleted or paid models used.
