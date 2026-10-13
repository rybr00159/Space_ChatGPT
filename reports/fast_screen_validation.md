# Fast broad survey validation

Recorded 2026-10-13T18:12:47.773679+00:00; deterministic RNG seed 20261013.

Reproduce with `python tests/test_fast_screen.py --benchmark`; run integrity tests with `python -m unittest discover -s tests -p test_fast_screen.py`.

## Search and cadence handling

This searches actual light curves across 0.3–13 days with a logarithmic period grid. The nearest trial period drifts by at most half the shortest trial duration over the observed baseline, including both grid endpoints. Four nominal durations (0.025, 0.05, 0.10, 0.16 days) are filtered to those >=1.5 effective cadences. For sector-1 QLP 30-minute cadence this removes 0.025-day boxes; 0.05/0.10/0.16-day boxes remain.

Quality-zero data are split at gaps >=max(30 minutes, 2.5 median cadences), binned without crossing gaps, detrended with a one-day robust median, and conservatively assigned errors at least as large as observed residual MAD. Edge exclusions of 0.06 days lose some observing completeness. Known ephemerides are masked before binning. Two sampled local events, two positive local depths, <=0.15 duty cycle, and no event contributing >85% of summed positive event S/N squared are required by the guard flag.

Broad mode uses `n_refine=0, n_peaks=2, compute_full_stats=False`; it still evaluates the full documented period grid and retains an independent local single-event scan. Detailed followup uses the default top-20 local grid refinement and five detailed peaks.

## Recovery and runtime

Thirty-minute synthetic injections: 36/36 periods recovered within 0.5% in the top two peaks; 36/36 including simple half/double aliases. Injections span randomized 0.4–12.5-day periods, 1.2/2.4/3.84-hour durations, requested ideal white-noise S/N 8/12/20, weak smooth variability, a 1.2-day downlink gap, and 3.5% missing cadences.

Synthetic broad-screen runtime: median 0.0429 s, maximum 0.0519 s per star on this machine (one worker; excludes downloads).

- Requested S/N 8: 12/12 exact-period recoveries.
- Requested S/N 12: 12/12 exact-period recoveries.
- Requested S/N 20: 12/12 exact-period recoveries.

L 98-59 c known positive control: 3/3 available two-minute sectors recovered near 3.6906764 days. This verifies sensitivity to a known planet; it is not a new candidate.

QLP TIC 8195168, sector 1: 1282 raw, 1000 quality-filtered, 816 searched cadences, 2063 periods, 0.0405 s. Best score 9.42; this benchmark does not assign a planetary interpretation.

QLP TIC 8195216, sector 1: 1282 raw, 1000 quality-filtered, 816 searched cadences, 2063 periods, 0.0410 s. Best score 5.91; this benchmark does not assign a planetary interpretation.

## Statistical limitations

64 white-noise nulls produced median/max best BLS S/N 4.69/5.55. These are runtime/behavior controls, not a sufficient survey-wide false-alarm calibration. A many-star search requires empirical null controls preserving correlated/instrumental noise, trial accounting, and independent confirmation of recurrence and source. No reported BLS score is a false-alarm probability.

This deliberately broad stage can miss shallow, shorter-than-resolution, single/double transit, boundary, very long-period, or strong-variability signals. A recovered harmonic must be resolved by event inspection. Failures at requested S/N 8 do not exclude planets. The injection test is small and synthetic; real-data completeness is not established.

Every measurement and case is in `results/fast_screen_validation.json`.
