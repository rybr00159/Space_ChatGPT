# Research checkpoint — broad survey continuation

Updated 2026-10-14. The fixed 160,000-target QLP sector-1 sample completed: **159,505 actual period searches, 495 insufficient-data cases, zero errors**. The interruption at 117,638 searched stars was resumed successfully. Read `results/broad/summary.json` and the audit for the stable final export. No credible new exoplanet candidate has yet been established.

## Active work and exact commands

- Survey completed; no live writer remains. Reconstruct/resume only if needed: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python src/broad_survey.py --download-workers 32 --analysis-workers 4 --queue 192 --report-every 2000`. Check active processes before launching; never start a second writer. Existing successful/insufficient targets are skipped. Database `/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite`; frozen `src/fast_screen.py` hash `26eede8aaeb9fa2a624b73a5c3966366f77c4cc73738153eba7e46d13eb10307` is enforced. Preserve incompatible checkpoints and start a separate experiment if changing the scanner.
- After the final stable export: `python src/summarize_broad.py` audits unique TICs, fixed-sample membership, status/manifest alignment and generates `results/broad/audit.json` plus an overview figure.
- Current stellar properties: `python src/refresh_shortlist_stars.py --include-vetted` queries all strict-eligible targets and previously vetted stars, saving current TIC responses and hashes. Then `python src/broad_prioritize.py --stellar-table data/broad_catalogs/current_stars_20261014/current_stellar.csv --catalog-radius-arcsec 120` creates a refreshed queue. Confirm current CLI and output paths in the source/report; preserve previous queue snapshots.
- Detailed follow-up: `python src/broad_vet.py --top 200` uses exact saved scan/peak IDs and version-aware checks. Current v3 source hash is `23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333`. Original 187-star batch completed under v3 with zero errors; 18 raw orbital-variation warnings were added. Old v2 evidence is preserved under `candidates/broad_s1/batch_history/`. Run the updated current-property full-sample queue, then investigate survivors independently.
- Remaining recognized-source pixel assessment: TIC 49617264 is already reported as TIC 49617264.11 / TOI-2417.11 in public 2024 ground follow-up. Native sectors 106/107 favor that faint star as the dip source. Preserve it as a recovered candidate and source-localization control; no novelty claim.

## Scientific findings already secured

TIC 370247818: eight sectors, 39 event timings and three independent pixel reductions support a real recurring source near the target. A raw secondary at phase 0.47546, grazing geometry and conditional light-budget checks favor an EB. TIC 141980042: 23 sectors/192 timings, strong secondary and crowded astrometry favor an EB/blend. TIC 308600649: strong raw half-period modulation and Gaia RV variation favor a stellar binary/blend. None is promoted as a planet.

Current TIC radii classify the four other focused leads (260352037, 441386747, 406280158, 271697430) as giants rather than the small stars in obsolete QLP headers. Fifty-seven retained light curves were independently checked. The 0.88-day orbit on TIC 260352037 would lie inside the claimed giant photosphere over the illustrative mass range; a blended dwarf remains a separate hypothesis. Additional secondary leads 79354493, 231825824 and 259901286 also have giant radii. Current stellar context must be applied before ranking planet sizes.

Sector-107 metadata: 846 public TCE rows on 623 stars, 71 eligible rows on 63 stars. Seventeen actual SPOC curves on five stars were independently investigated, distinct from broad-survey counts. Three variability interpretations fail, one star is a published complex periodic variable, and the remaining recurring candidate (49617264) was already followed in 2024. See `candidates/tce107/` and the exact prior-report proof in `reports/tic49617264_prior_check.md`.

The realistic fast-search injection experiment recovers period+epoch in 26/64 at 1,000 ppm and 49/64 at 3,000 ppm; local cuts retain 11/64 and 27/64. Structured surrogates produce false-looking leads. These trials do not establish global significance or population completeness. Longer periods, weak events, short durations and variability-dominated stars remain incompletely searched.

## Environment, preservation and resumption

Use `/workspace/openrouter-check-venv/bin/python` in the existing `/workspace/Space_ChatGPT` checkout. Set `TESS_DATA_DIR=/workspace/tess-research-data`, `MPLCONFIGDIR=/workspace/tess-research-data/matplotlib`, and both thread variables to 1. Eleven utility tests pass; `pip check` and astronomy/BATMAN imports pass. `bash scripts/install_environment.sh` installs pinned dependencies when required. No web service is needed.

Raw products, pixels, SQLite and intermediates are retained outside Git. `results/broad/observations.jsonl.gz` and follow-up manifests record exact archive products and SHA-256. `python src/restore_broad.py --limit 20` verifies cached files without HTTP; explicit `--download --max-downloads N` restores missing files with integrity checks and a storage reserve. Existing mismatches are preserved and block analysis. QLP intrinsic checksum cards are absent, recorded honestly; products containing checksum cards are verified. The original pilot uses `python src/retrieve_manifest.py` separately.

This overlay filesystem has no independently verified durable volume. Updated start instructions are saved in the environment configuration draft. Review/save and publish the environment to snapshot prepared raw files; draft saving alone does not publish, and live jobs never should be assumed persistent. Git preserves compressed scientific evidence and reproducible retrieval instructions if the cache is lost; fresh-task restoration has not been tested.

Research OpenRouter usage remains zero. The earlier connectivity-test charge was not recorded; reserve for it before any future spending. User-authorized cumulative cap: $10. Never display, log or commit secrets.

## Next scientific phase after the current queue

Inspect the final full-sample current-property ranking and detailed-vetting survivors. Retrieve independent sectors and pixel evidence for any signal retaining raw recurrence, plausible stellar geometry, and no decisive prior recognition. Validate with candidate-specific alternative baselines, event/timing models and contamination checks; absence from one catalog is insufficient. If no survivor remains, prioritize a preselected quiet-star multi-sector experiment with several detrending scales, more retained independent peaks, and weaker/longer-period injection recovery. Do not simply repeat the same single-sector coarse screen and claim improved completeness.

Update `BROAD_RESEARCH_REPORT.md`, README, log, candidate assessment index and exact counts after each phase. Commit meaningful progress and push `main` as already authorized. Original 19-star pilot details remain in `RESEARCH_REPORT.md` and `reports/pilot_progress_20261013.md`; those historical next steps are superseded by the broad continuation.
