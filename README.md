# Autonomous TESS exoplanet research

**Current result: no credible new exoplanet candidate in the completed 2026-10-13 pilot.** We searched real observations and investigated stronger signals; apparent detections failed raw-flux, local-baseline, recurrence, variability, or noise checks. This result does not exclude additional planets in the targets.

Read [the research report](RESEARCH_REPORT.md), [research log](RESEARCH_LOG.md), and [resumption checkpoint](PROGRESS.md). The original `Reddit Post.zip` is retained unchanged; its publicly reported signals are excluded from novelty claims.

## What was investigated

| Quantity | Completed |
|---|---:|
| Targets: eight known hosts, ten field dwarfs, one control | 19 |
| Public SPOC TESS light curves | 66 |
| Per-sector BLS variants, including three positive-control scans | 135 |
| Additional untrimmed contiguous-sector variants | 16 |
| Gap-edge-excluded variants across the pilot | 58 |
| Ranked screening peaks retained in the signal table | 775 |
| Blind injection/control search runs | 8 |
| Illustrative block-permutation full-grid searches | 96 |
| Additional selected-star regional null searches | 64 |
| Utility tests | 5 passing |
| Promoted new candidates | 0 |

Search additional planets after masking known ephemerides, and compare with deterministically selected bright cool dwarfs outside current planet catalogs. Two median detrenders, gap-edge checks, individual-event models, raw SAP comparisons, and recurrence tests attempt to falsify the strongest leads. More details and limitations are in the report; scores are not global false-alarm probabilities.

## Installation and execution

Python 3.12 with pinned packages in `requirements.txt`. The prepared cloud environment uses `/workspace/openrouter-check-venv` (the name comes from its initial connectivity check). From `/workspace/Space_ChatGPT`:

```bash
bash scripts/install_environment.sh
export MPLCONFIGDIR=/workspace/tess-research-data/matplotlib
export TESS_DATA_DIR=/workspace/tess-research-data
source /workspace/openrouter-check-venv/bin/activate
python -m unittest discover -s tests -v
python src/retrieve_manifest.py
```

The installation script was executed successfully, with dependency checks and astronomy imports. Optional Lightkurve PRF functionality is unavailable without oktopus; the analysis uses Astropy FITS/BLS and does not require it. No web service or paid model is needed.

Raw cache verification uses exact recorded URLs, SHA-256 and raw FITS HDU checksum validation. Download only missing files; preserve integrity mismatches for diagnosis. Product paths in the manifest reflect this cloud workspace. If working elsewhere, map local paths to the restored cache before invoking analysis; archive IDs and hashes remain unchanged.

The completed screening commands are repeatable and skip existing results:

```bash
python src/search.py --width 1 2
python src/joint_search.py
python src/robust_search.py
python src/summarize.py
```

`joint_search.py` would fill additional untrimmed contiguous groups if rerun with the final expanded manifest; the saved 16 variants reflect the earlier incremental download checkpoint. The gap-edge search already covers the final pilot. New scans are additional experiments; update the log and tables if executed. Do not repeat all scans unnecessarily.

Validation and dedicated vetting scripts are in `src/`. `python src/null_calibration.py --frozen-search` resumes the exact historical null experiment using the hash-verified search implementation from commit `7ecc52c`, without repeating completed permutations. Production defaults deliberately reject incompatible checkpoint versions.

## Results and evidence

- [Analyzed targets](results/targets.csv) and [searchable signal ledger](results/signals.csv): all ranked peaks are screening entries, not accepted planets.
- [Candidate/evidence index](candidates/README.md): stronger investigated signals, measurements and diagnostic plots.
- [Known-transit and injection validation](reports/validation.md) and [empirical noise comparisons](reports/null_calibration.md).
- [Catalog and literature scope](reports/catalog_access.md), [field selection](reports/field_selection.md), and [TIC 150299576 prior checks](reports/tic150299576_prior_check.md).
- [Exact observation manifest](data/observations.json): TIC, sector, MAST observation/product identifiers, processing metadata, retrieval URL and SHA-256.

## Storage and scientific limits

Raw FITS and reusable intermediates are retained outside Git at `/workspace/tess-research-data`, about 154 MB at the checkpoint. Compact candidate curves, plots, code, manifests and catalog snapshots are tracked. Nothing scientifically important was deleted. `/workspace` is overlay, not a separately verified durable volume. Environment publication can snapshot prepared files; new-task restoration was not verified, and processes are not persistent. Manifest-based retrieval protects reproducibility if the cache is lost.

The pilot incompletely samples shorter durations, periods beyond 25 days, single transits, mask/TTV uncertainties, stellar variability and population completeness. No surviving signal justifies a physical planetary fit or discovery claim; pixels/centroid localization and a full multi-target false-positive analysis remain unrun. Known recoveries and favorable injections do not establish survey completeness. Catalog/literature absence does not establish novelty.

No OpenRouter calls were made during research; current research spending is $0 against the $10 authorized cap. An earlier environment connectivity call predates this project and its charge was not recorded. Never log or commit API keys.
