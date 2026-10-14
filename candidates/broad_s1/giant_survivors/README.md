# Stale-radius giant-star survivor checks

These four initial flag-free weak leads used obsolete QLP stellar-radius metadata. Current TIC v8.2 gives radii of 10–12 solar radii. The previously inferred small planet sizes are unsupported. Cached additional observations are preserved and were checked with independently optimized local periods and epochs.

|TIC|Trial P (days)|Old QLP R|Current TIC R|Sectors checked|a at 3 solar masses / current R|
|---|---:|---:|---:|---:|---:|
|260352037|0.88013|0.781581|10.1065|25|0.55|
|441386747|7.48010|0.733908|12.1823|2|1.91|
|406280158|10.12590|0.781807|10.5196|3|2.70|
|271697430|3.05937|0.742671|10.2269|25|1.25|

Keplerian axes are an illustrative 1–3 solar-mass range, not fitted masses. An orbit inside the photosphere cannot produce an external planetary transit on the giant. A signal originating on a blended dwarf remains possible and requires pixel source attribution. Where independently selected broad-search peaks recur, the optimization and giant correlated noise prevent a planetary interpretation from the score alone.

Per-sector depths, event sampling, raw phase coefficients, archive source URLs, and SHA-256 identifiers are in the QA JSON files. All original FITS remain cached outside Git. Reproduce with `python candidates/broad_s1/giant_survivors/recurrence_qa.py`.
