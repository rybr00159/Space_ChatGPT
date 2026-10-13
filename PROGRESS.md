# Progress checkpoint

Updated 2026-10-13: **first research phase complete; no credible new exoplanet candidate.** This is a negative result in a limited pilot, not an exclusion of additional planets.

## Completed

- 19 stars: eight known two-planet hosts, ten independently selected field dwarfs, L 98-59 control; 66 public SPOC curves.
- 135 per-sector variants (three unmasked controls), 16 additional untrimmed joint variants from an incremental archive snapshot, 58 edge-excluded variants over the final sample; 775 ranked peaks in `results/signals.csv`.
- All three control planets recovered. Eight blind injection/control scans, 96 full-grid block permutations and 64 selected-star regional null runs completed. Five utility tests pass.
- Stronger leads investigated through raw SAP, local linear/quadratic models, recurrence, gap checks, variability/harmonic models, literature and contamination context. See `RESEARCH_REPORT.md` and `candidates/README.md` for results and precise limitations. No promoted candidate or physical planetary fit.
- Raw FITS plus intermediates retained in `/workspace/tess-research-data` (~154 MB); exact identifiers, processing metadata and hashes in `data/observations.json`. Compact candidate curves and plots tracked in Git.
- Install script and research start instructions saved in environment draft; installed script tested. Publication and fresh-task restoration remain unverified.
- Research-phase OpenRouter calls/spending zero; no other paid services used. Earlier connectivity test predates project and its charge was not recorded. Maximum authorized cumulative OpenRouter budget is $10; no paid calls planned.

## Resume from this checkpoint

1. Use existing checkout `/workspace/Space_ChatGPT`; do not create a worktree. Check Git status and read the report/latest research log before new work.
2. Interpreter `/workspace/openrouter-check-venv/bin/python`; activate it or run `bash scripts/install_environment.sh` if needed. Export `TESS_DATA_DIR=/workspace/tess-research-data` and `MPLCONFIGDIR=/workspace/tess-research-data/matplotlib` before imports.
3. `python src/retrieve_manifest.py` verifies/restores exact recorded products with SHA-256 and raw FITS checksums. Preserve mismatches for diagnosis; do not rewrite expected hashes or disable verification. All 66 products were verified at final checkpoint.
4. Existing scans and candidate tests are complete: reuse results rather than repeat them. `src/summarize.py` reconstructs searchable tables. `src/search.py` and `src/robust_search.py` skip completed variants. `src/joint_search.py` can fill additional untrimmed groups under the final expanded manifest; that is optional additional analysis, not missing final gap-edge coverage.
5. Historical null reproduction: `python src/null_calibration.py --frozen-search`, which loads the exact hash-verified `search_null_v1.py` from commit `7ecc52c`. It resumes all 96 completed permutations without repeats. Default production mode intentionally rejects that older source checkpoint.
6. Utility readiness: `python -m unittest discover -s tests -v`. Update log/report/tables and commit/push meaningful new work to main. Preserve raw observations and scientifically important intermediate evidence.

## Exact next scientifically justified experiment

Expand evidence on the quiet field subset before merely adding many similar stars:

1. Start with TIC **150359842, 278727563, 260414764, 260476305**. Their present global residual maxima are weak, with no promoted lead. Retrieve their already recorded **sectors 87–90** using `src/archive.py` query metadata (add explicit sector selection first; do not assume `--limit` selects those sectors). Retain the current sectors 93–96 as a separate era.
2. Add a gap-aware detrending alternative and sensitivity to **30–60 minute transits**, preserving higher-cadence observations before binning. Current production durations start at 1.08 hours and 10 minute bins attenuate very short events.
3. Before interpreting new ranks, inject across randomly selected period/epoch/radius/duration and real gap/known-mask coverage, with phase/detrending selection frozen on training data. Estimate completeness and empirically assess all tested star/period/detrender choices. Existing favorable injections do not measure population completeness.
4. Search each era separately and jointly within appropriate baseline/grid budgets. For any survivor, get exact TPF products and independent pixel reductions, inspect centroids/difference images and neighboring catalog stars, test odd/even and secondary features, fit a physical model with uncertainty, and refresh TOI/CTOI/confirmed/literature checks. No surviving signal currently requires TPF downloads.
5. If this yields no survivor, broaden to a new southern cone with selection before flux inspection or investigate newly released sectors on small known systems. Longer-period/single-transit methods remain a separate research track.

This phase provides no statistically valid yield comparison between the two purposeful tracks. Broader searches are future experiments, not evidence that a discovery is imminent.

## Preservation contract

`/workspace` is overlay, not a verified durable volume. Cloud publication can snapshot prepared files; no running processes are needed for this project and none should be assumed persistent. Review/save and publish the environment draft to preserve prepared runtime state. Exact Git manifests and retrieval scripts provide recovery if the cache is lost; restoration in a new task was not tested. Never put keys or raw credentials in files, logs or Git.
