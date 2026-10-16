# Separate empirical-uncertainty recovery search

This is an explicitly labeled alternative-assumption analysis of observations
the frozen primary multi-sector survey could not search. The primary database,
original FITS products, producer error columns and quality flags are unchanged.
These are **additional period-searched observations**, overlapping the primary
survey's selected targets; do not count them as new unique targets or as a
second independent confirmation of any signal.

## Why recovery is necessary

Numerous selected Sector 95 QLP products contain thousands of quality-zero,
finite positive corrected-flux samples but no finite positive `DET_FLUX_ERR`.
The frozen screen correctly excludes unusable uncertainties and therefore
reports `insufficient_data`, even though the flux could be examined under an
explicit empirical-noise assumption.

For example, TIC 234490958 has 8,517 quality-zero finite positive flux samples;
all 8,517 associated producer uncertainties are nonfinite. The recovery estimate
uses 8,494 short-gap adjacent differences, giving a constant uncertainty of
0.0018190254 in corrected-flux units, or 1,818.9 ppm of median corrected flux.
The original FITS SHA256 and the diagnostic are in `estimator_validation.json`.

## Selection and uncertainty assumption

`src/multisector_missing_errors.py` opens the main SQLite database in read-only
mode, copies a consistent set of existing `insufficient_data` records, and
closes the read transaction immediately. It considers only observations with:

- At least 400 finite positive corrected-flux samples with finite times and
  `QUALITY == 0`.
- Fewer than 20 of those samples having finite positive producer uncertainties.

Other insufficient-data records are audited as `ineligible` and are not searched.
Products having valid positive producer errors remain in the primary method.
For each eligible target-sector observation independently, sort the good flux
samples by time. Form adjacent flux differences only where their separation is
positive and no greater than 2.5 native cadences. Require at least 100 such
differences. The cadence is the positive FITS `TIMEDEL`, with an explicitly
implemented median-positive-time-difference fallback if that header is absent.
Then estimate

```
sigma_flux = 1.4826 * median(abs(delta_flux - median(delta_flux))) / sqrt(2)
```

This single positive value supplies an analysis-only constant error array for
the entire target-sector observation. It is **not** written into the original
product or presented as a repaired producer error column. The frozen
`multisector_screen.screen_arrays` performs the same two detrending searches
(0.75 and 1.5 days), retaining its observed-scatter uncertainty floors. Quality
flags are passed through unchanged.

The adjacent-difference estimator assumes comparable, locally independent
scatter. Correlated noise can bias it low, rapid variability can bias it high,
and a constant error loses heteroscedastic information. Scatter floors reduce
some weighting problems but do not establish calibrated statistical significance.
The reported S/N values remain ranking statistics without a survey-level
false-alarm probability. Two detrenders of the same flux are correlated.

## Outputs, integrity and resume

- `summary.json`: completed recovery counts, assumptions and exact method hashes.
- `targets.csv`: usable-flux/error diagnostics, noise estimates and raw hashes.
- `signals.csv`: every retained trial peak from both frozen detrenders.
- `candidate_seeds.csv`: peaks passing the same permissive seed criteria, still
  requiring false-positive checks, source localization and recognition review.
- `all_scans.jsonl.gz`: complete scan records, including individual events,
  original-record SHA256 and source/raw-product hashes.
- `estimator_validation.json`: a synthetic 1,000 ppm white-noise control with a
  masked gap and extreme outlier recovered 999.57 ppm; a normal curve with
  positive producer uncertainties was excluded. This validates units and basic
  robustness, not the real-data noise assumptions or astrophysical nature.

The separate incremental SQLite checkpoint is
`/workspace/tess-research-data/multisector/error_recovery.sqlite`. It records each
observation immediately, exports every 20 records, and refuses mixed source
versions. The original raw-product SHA256 is checked before each analysis.
`all_scans.jsonl.gz` preserves the full numerical records in Git if the local
database is unavailable in a future environment. Producer raw-file URLs allow
re-retrieval, subject to the saved exact-byte checks.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  /workspace/openrouter-check-venv/bin/python \
  src/multisector_missing_errors.py --report-every 20
```

Running again adds newly completed insufficient-data rows from the primary
survey and skips records already in the recovery checkpoint. There is no
network access or paid-service call in this recovery. Do not append results
from a changed method to the same database.

Any independent follow-up must preserve this uncertainty label. Reading the raw
file with `read_curve` still returns the original invalid uncertainties; use
the saved empirical-noise estimate to reconstruct the alternative error array
explicitly. Discovery, novelty and source attribution have not been established
by a retained seed.
