# Illustrative empirical BLS search maxima

These permutations compare each observed light curve with maxima over the exact same production period/duration/epoch grid. They offer a caution about ranking statistics; they do not calibrate discovery probabilities or the full pilot search.

Residuals are obtained with the 1-day and 2-day median detrenders, 10-minute bins, known-planet masks, and 0.25-day removal around quality-filtered raw-data boundaries/gaps. Half-day contiguous blocks are randomly reordered onto the original timestamps. Original and shuffled data use identical constant errors equal to median diagonal errors, avoiding relocation of heteroscedastic errors. The original diagonal-error maximum is also retained in the JSON for comparison.

|TIC|Sector|Width (d)|Observed maximum S/N|Null n|Median null maximum|95th percentile|Null ≥ observed|Plus-one tail fraction|
|---|---|---|---|---|---|---|---|---|
|278727563|94|1.0|5.12|24|4.93|5.24|3|0.160|
|278727563|94|2.0|5.23|24|4.96|5.33|3|0.160|
|307210830|2|1.0|5.76|24|5.25|5.69|1|0.080|
|307210830|2|2.0|5.60|24|5.90|6.26|21|0.880|

Limitations:

- Only selected individual sectors of two targets are examined. These fractions are not global false-alarm probabilities across stars, sectors, contiguous-sector searches, detrending choices, or selected candidates. Choosing a quiet comparison star after inspecting data also limits generalization.
- With 24 permutations, the smallest plus-one fraction is 1/25 = 0.04; tail probabilities cannot support high significance. The empirical 95th percentile is unstable at this sample size.
- Permuting already flattened blocks fixes the detrending realization. It does not propagate uncertainties in stellar variability removal or adaptive masking. Correlations within a half-day block survive; longer correlations and their relationship to gaps are disrupted, and block joins introduce new discontinuities.
- The procedure preserves observed residual excursions, including possible real events; it is not a null guaranteed to contain only stochastic noise. It assumes approximate exchangeability of blocks and cannot simulate spacecraft artifacts tied to particular timestamps.
- A scalar median uncertainty slightly changes the observed ranking relative to production heteroscedastic weights. Inference uses observed and shuffled maxima under the same scalar convention.
- Widths use corresponding seed/permutation when sampling agrees. Per-width summaries include the full period grid; they do not automatically correct for selecting the maximum across both widths. Paired maxima across widths are summarized below where available.

TIC 278727563, sector 94: maximum across both widths is 5.23; 5/24 paired null maxima are at least this large (plus-one fraction 0.240). This includes the two tested widths only and remains subject to the limitations above.

TIC 307210830, sector 2: maximum across both widths is 5.76; 18/24 paired null maxima are at least this large (plus-one fraction 0.760). This includes the two tested widths only and remains subject to the limitations above.
