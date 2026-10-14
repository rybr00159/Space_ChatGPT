# Autonomous TESS exoplanet and stellar-system research

**Stellar-system phase:** [TIC 273879360](candidates/stellar_systems/tic273879360/README.md) is a strong potentially unreported background eclipsing binary: approximately 21–27% primary and 7–9% secondary eclipses repeat on this source in four usable pixel sectors. A second coherent short-period eclipsing signal remains unresolved in the image of giant [TIC 260352037](reports/stellar_systems/hidden_sources.md). Existing source positions were already catalogued; component masses and first-discovery status are not established. Read [the stellar research report](STELLAR_RESEARCH_REPORT.md).

**Completed broad screening: 159,505 stars actually period-searched; 495 additional light curves were insufficient. All 389 shortlisted signals have received detailed checks. No credible new exoplanet candidate has been established.**

The earlier 19-star pilot found no credible new candidate. We searched real observations and investigated stronger signals; apparent detections failed raw-flux, local-baseline, recurrence, variability, or noise checks. This result does not exclude additional planets in the targets.

Read [the broad research report](BROAD_RESEARCH_REPORT.md), [historical pilot report](RESEARCH_REPORT.md), [research log](RESEARCH_LOG.md), and [resumption checkpoint](PROGRESS.md). The original `Reddit Post.zip` is retained unchanged; its publicly reported signals are excluded from novelty claims.

**Targeted claim audit:** [TIC 71316387 / CTOI 71316387.01](candidates/tic71316387/README.md) shows persistent alternating eclipses in five public QLP sectors. A roughly 1.30367-day eclipsing system is favored over ordinary transits at the catalogued 0.65184-day period. The specific quoted 5,293/2,487 ppm and 19.6σ calculation is not independently certified; conservative estimates and independent pixel evidence are preserved. This is an audit of an existing candidate, not a new planet discovery.

## Broad survey and detailed follow-up

Official sector-1 QLP inventory: 469,735 targets. A fixed seeded random sample of 160,000 was selected before flux inspection. All 160,000 selected targets completed: **159,505 actual period searches, 495 insufficient-data cases, zero processing errors**. Inventory, selected targets, screened stars, peaks and detailed follow-up are distinct counts, audited in `results/broad/audit.json`. Download [the compressed target table](results/broad/targets.csv.gz), [signal table](results/broad/signals.csv.gz), and [observation manifest](results/broad/observations.jsonl.gz). Every downloaded raw FITS is retained under `/workspace/tess-research-data/broad/qlp_s0001/raw`; SQLite checkpoints completed stars. The compact 30-minute observations permit this sample to fit the available disk but attenuate short transits. Newer timing-flagged QLP sectors are excluded.

`src/broad_survey.py` pipelines pooled HTTP downloads and four CPU workers. The validated logarithmic BLS grid searches 0.3–13 days, with cadence-resolved duration trials and local single-event screening. Fine grids, alternative detrending, raw/aperture comparisons and additional sectors are reserved for leads. Novelty checks add published DIAmante, Olmschenk/Yu and genuine TESS EB catalogs to NASA/TOI/CTOI. Official sector-107 TCE triage runs separately; archive-reported signals are not our original detections.

The exact survey can resume missing targets from SQLite/cache using the command below; it has no pending targets at this checkpoint:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python src/broad_survey.py --download-workers 32 --analysis-workers 4 --queue 192 --report-every 2000
```

Do not launch a second writer while that command is active. The source hash is enforced to prevent mixing analysis versions. See `reports/fast_screen_validation.md` and `reports/broad_strategy.md`. No discovery claim is supported by broad screening alone. The original 187-target detailed batch was refreshed with raw-flux orbital-variation tests; current TIC properties reveal several old QLP dwarf estimates are giants. The full-sample review checked 389 exact signals on 357 TICs using current stellar properties and wider nearby-source checks. The two strongest newly investigated recurring eclipses admit equally good doubled-period binary models, with alternating depths supporting that interpretation. Another shallow short-period signal is independently localized to a known faint eclipsing neighbor. See [all detailed signal evidence](candidates/broad_s1/all_signal_evidence.csv). Multi-sector and pixel checks favor eclipsing binaries for TICs 370247818 and 141980042. The recent-sector track independently examined 17 SPOC curves on five stars and recovered previously reported candidate TOI-2417.11; this is not a new discovery. See [candidate evidence](candidates/README.md), [injection/noise limits](reports/broad_noise.md), and [raw-data restoration](reports/broad_preservation.md).

## Completed initial pilot

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

Raw FITS, pixels, SQLite and reusable intermediates are retained outside Git at `/workspace/tess-research-data`. The broad sample alone retains 160,000 FITS files totaling approximately 12.9 GB; exact bytes are in the audit. The pilot cache was about 154 MB before expansion. Compact candidate curves, plots, code, manifests and catalog snapshots are tracked. Nothing scientifically important was deleted. `/workspace` is overlay, not a separately verified durable volume. Environment publication can snapshot prepared files; new-task restoration was not verified, and processes are not persistent. Manifest-based retrieval protects reproducibility if the cache is lost. Four compressed full-database chunks in `data/broad_checkpoint/` preserve all original search results; a 160,000-row exact restore was tested without repeating analysis. Run `python src/survey_checkpoint.py verify` and restore to an unused database path when needed.

The pilot incompletely samples shorter durations, periods beyond 25 days, single transits, mask/TTV uncertainties, stellar variability and population completeness. The pilot had no survivor justifying a physical planetary fit. Broad-phase pixel localization and physical eclipse fits were executed for strong leads, but do not establish planetary origin. A global survey false-positive analysis remains unresolved. Known recoveries and favorable injections do not establish survey completeness. Catalog/literature absence does not establish novelty.

No OpenRouter calls were made during research; current research spending is $0 against the $10 authorized cap. An earlier environment connectivity call predates this project and its charge was not recorded. Never log or commit API keys.
