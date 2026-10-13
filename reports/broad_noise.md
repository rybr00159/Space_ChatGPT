# Real-observation screening calibration

Updated 2026-10-13T21:37:50.773705+00:00. Frozen broad-screen source SHA-256:
`26eede8aaeb9fa2a624b73a5c3966366f77c4cc73738153eba7e46d13eb10307`. Screening config is exactly `{'n_refine': 0, 'n_peaks': 2, 'compute_full_stats': False}`.
No active survey source, database, or selection was modified.

The first 128 successfully screened, catalog-unflagged
stars in the pre-flux randomized survey order were selected before inspecting
noise outcomes; no S/N, variability, magnitude, or radius selection was applied.
Excluded host flags: is_confirmed_host, is_toi_host, is_ctoi_host, is_eb_host, is_diamante_host, is_published_host. Snapshot selection exclusions and
catalog hashes are in `results/broad/noise_assessment.json`. Unknown variable,
eclipsing, and planetary signals may remain. Selection conditional on successful
screening and catalog exclusions limits extrapolation to all 160,000 stars.

Completed 128 original real-light-curve screens, 512
surrogate screens, and 128 raw-flux injections. TESS magnitude
range/quantiles: `{'minimum': 8.56, 'p10': 10.3676, 'median': 12.563500000000001, 'p90': 13.3137, 'maximum': 13.426}`; detrended residual MAD in ppm:
`{'minimum': 138.15587736641723, 'p10': 366.8014070755763, 'median': 1572.3048458078915, 'p90': 3535.6671654256456, 'maximum': 19037.820353696123}`; adjacent-point correlation:
`{'minimum': -0.05546811425359081, 'p10': -0.013813425138145439, 'median': 0.04250863504850627, 'p90': 0.503620408290983, 'maximum': 0.957047192646234}`.

Each target has two within-gap-segment block permutations and two random block
sign reversals. Residuals/errors are transformed relative to the original 1-day
median trend and reconstructed in raw flux units; the frozen pipeline then
detrends and searches them afresh. Blocks are approximately half a day (12
points at 30-minute cadence). Permutations move paired residuals/errors;
sign flips retain errors. Original timestamps, quality flags, gaps, and nominal
exposures are retained. Blocks preserve internal short correlations but create
boundaries and alter longer variability, so these are approximate surrogate
noise/astrophysical-mixture tests, not exchangeable null realizations.

Initial-guard plus diagonal S/N >= 8 lead rate: `{'count': 9, 'trials': 512, 'fraction': 0.017578125}`.
Stronger local recurrence guard rate: `{'count': 13, 'trials': 512, 'fraction': 0.025390625}`.
The stronger guard requires local diagonal S/N >= 8, >=3 sampled and positive
events, >=70% positive events, <=65% event-power dominance, positive depth,
and the original guards. It omits catalog and stellar-physics exclusions.
The BLS-S/N and local-S/N sets are not nested: the local guard does not impose
a separate BLS S/N threshold, consistent with the survey prioritizer.
Rates count at least one of the top two peaks in each surrogate, rather than
counting every trial period as independent. They are empirical lead fractions,
not candidate false-positive probabilities. Repeated trials per star are
dependent. Null maxima are `{'minimum': 3.5533616162630675, 'p10': 4.290620086226022, 'median': 4.836390835716215, 'p90': 7.9285506119395075, 'maximum': 200.77270869351844}`.
Star-level counts with any of four transformed trials passing are
`{'any_initial_guard_pass_snr_ge_8': {'stars_with_any_of_four': 5, 'stars': 128}, 'any_strong_local_guard_pass': {'stars_with_any_of_four': 6, 'stars': 128}}`. Sign reversal can turn a positive flare/outlier
normally removed by preprocessing into an artificial dip: the largest maxima
(TIC 214299961) fail the recurrence guards, and must not define a threshold.
The transformed distributions are not exchangeable with the original data;
their rates quantify screening vulnerability to structured/asymmetric noise.

The first 64 selected stars each receive matched
1000 and 3000 ppm box injections with 2.4-hour duration, independent frozen
random phase, and log-uniform period 0.5--12 days. The box is integrated across
the 30-minute exposure before multiplying the actual raw KSPSAP flux. No gaps
are filled; existing noise, contamination, variability, and exposure loss remain.
Error arrays are unchanged. Injection parameters are fixed in the sample file
before trial outcomes. Recovery requires period within 1% and epoch agreement
within half the sum of injected/recovered durations plus one exposure; separate
period-only, simple-alias, initial-guard, and stronger-local metrics are retained.

Recovery counts by depth (denominators are explicit):

| Depth ppm | Injection pairs | Period + epoch | With initial guards | With local recurrence guards |
|---:|---:|---:|---:|---:|
| 1000 | 64 | 26 | 20 | 11 |
| 3000 | 64 | 49 | 41 | 27 |


Recovery is also stratified in the JSON by fixed descriptive period boundaries
3/6 days, residual-MAD boundaries 1000/3000 ppm, TESS magnitude boundaries
11/12, and observed-event-count boundaries 3/6. Small strata are descriptive;
they are not a general occurrence-rate or completeness calibration. The
stratification was added after the first aggregate results, for descriptive
heterogeneity assessment; it was not used to select targets or injection pairs.
The
exposure duration is verified from the QLP TIMEDEL header as 1/48 day.

In this fixed sample, for residual MAD below 1000 ppm, 1000-ppm period/epoch
recovery is 16/18
(11/18 also meet local recurrence).
Above 3000-ppm residual MAD, recovery is 0/8
at 1000 ppm and 0/8 at 3000 ppm.
For periods longer than 6 days, exact period/epoch recovery is
0/7 at 1000 ppm and
2/7 at 3000 ppm;
0 and 0, respectively,
pass the local recurrence cuts. These small descriptive groups show that coarse
screening is strongly incomplete for weak, noisy, and long-period signals.
Only two top peaks are retained, so stronger existing variability can displace
injections. A negative survey screen cannot rule out planets in these regimes;
multi-sector data, multiple detrending scales, and searches retaining more
independent peaks are scientifically justified next steps.

Even an ideal 512-trial empirical null could resolve a plus-one tail probability
only to 1/513 ~= 0.00195; the 160,000-star Bonferroni 5% target is 3.125e-7.
This experiment cannot establish survey-global significance, and trial/star
dependence and null mismatch make that bound more optimistic than justified.
Candidate-specific raw-flux, independent-sector, pixel, contamination, and
catalog/literature checks remain necessary. No discovery or global FAP is claimed.

Run/resume with `/workspace/openrouter-check-venv/bin/python src/broad_noise.py`.
Independent atomic checkpoints: `/workspace/tess-research-data/broad/noise`. Raw files are reused and SHA-256
checked against the survey records; exact target URLs/hashes and input catalog
hashes are retained in the repository JSON. Nothing is downloaded or deleted.
Paid external requests: zero.
