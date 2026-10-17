# Research checkpoint — multi-sector exoplanet and stellar investigations

## Ground follow-up scope and access — 17 October 2026

Prepared [an observer-facing brief and cost/access comparison](reports/tic214311908_ground_followup_brief.md): separate the target and approximately 3.15-arcsecond neighbor photometrically through a complete event, with roughly five hours of time series/baselines. One good red-filter event can localize the dip; multiband measurements add blend tests. This does not alone measure mass or eliminate unresolved companions. Southern-site visibility, seeing/sampling, the nearby bright Mira, precise timing and observer-verified setup matter.

A willing observer/research collaboration or suitable allocated time may entail no direct personal observing charge; access is not guaranteed. Current LCO rules are programme-specific and are not a universal free standard-proposal route. Commercial rental is optional: iTelescope's checked standard-plan prices span roughly US$19–118/hour plus US$39.95/89.95 membership per28days, but no suitable instrument or actual quote has been selected. Five public access/price sources are retained in `data/ground_followup_options/20261017/`. No observer contacted, application sent, booking made or money spent.

## Publication and ExoFOP route checked — 17 October 2026

Journal shortlist is appended to the publication note. OJAp is the recommended first option to investigate: broad exoplanet scope, conventional peer review, no author charges and explicit AI-use disclosure/validation rules; arXiv-overlay workflow. AJ is a conventional observational/exoplanet fit but 2026 publication charges start at US$1,425. JAAVSO is a potential focused-paper inquiry option; peer-reviewed status is confirmed from AAVSO, but its current author pages returned403, so fees and candidate-only scope remain unverified. No comparative acceptance probabilities are claimed. Six public source snapshots and hashes are retained in `data/journal_options/20261017/`; no editor was contacted or submission made.

Read [the current-policy and publication-readiness note](reports/tic214311908_publication_path.md). Official ExoFOP announcement dated **19 August 2026** requires publication in an online peer-reviewed journal before new Community Candidate upload, plus approval through Published Candidate Upload Request. Existing account access alone is insufficient. The detailed guidelines explicitly allow an RNAAS note if the detection/vetting methodology was previously published in an online peer-reviewed journal and cited in the note; RNAAS itself is not peer reviewed. Do not assume a GitHub report, preprint or generic BLS citation qualifies.

Current evidence justifies developing a candidate paper and seeking expert review/follow-up; acceptance is not established. Before journal submission, make final deblended products/calibration/injections reproducible, bound systematic and false-positive claims, and complete paper-level novelty/reference checks. Resolved multiband transit photometry remains the most useful new observation; it is scientifically recommended rather than a blanket formal requirement for candidate publication. Official page snapshots/hashes are in `data/exofop_policy/20261017/`. No credentials, login, candidate creation, request, email or external submission were used.

## Latest assessment — comparison of TIC 214311908 investigations, 17 October 2026

The user authorized access to the other model's findings after the independent review and supplied the GitHub directory. Reviewed `rybr00159/Space` commit `99d4f66fca99c32c3245396f367b8a73d5d13fa5`. **Strong, unvalidated planet candidate is now a reasonable classification.** Read [the comparative report](candidates/tic214311908_external_comparison/README.md), [host/FPP audit](candidates/tic214311908_external_comparison/host_audit.md), [detection audit](candidates/tic214311908_external_comparison/detection_audit.md) and [new independent physical checks](candidates/tic214311908_external_comparison/shape_checks.json).

The six-sector calibrated PRF/WCS/colour analysis is substantially stronger than our initial two-sector Gaussian quicklook. Their source-specific background binary tests address the faint neighbor; our simple 50%-eclipse flux-budget observation was incomplete. New fits to our original 13 QLP windows/1961 samples independently disfavor two specified circular K-star binary models by conditional light-curve Δχ² 74.7/102.3. This is not a complete binary-population calculation. Large grazing-body/eccentric alternatives and unresolved hosts remain relevant; no mass, validation, confirmation or novelty is established.

The supplied nominal FPP 0.0013–0.0033 is conditional on its custom model. Published posterior-table arithmetic reproduces, but a complete independent PRF/evidence/FPP rerun was not performed. Its strict companion-host definition gives approximately 0.005–0.009 at adopted priors; joint prior draws have strict median 0.0127 and 5–95% range 0.00222–0.106. Do not quote an exact 99.9% planet probability as a calibrated observational conclusion. Original blind detection and holdout timing records remain unchanged.

Corrected a real error in our catalogue dossier: the saved multi-table El-Badry response already lists target Gaia 5925373272850784768 and companion 5925373272850782208, separation 756.8 AU and catalogue R 1.993e-5. We initially read the empty shifted-control table. The original response is retained and a new all-table/string-ID rereading check reproduces the positive entry.

Public source snapshots: 98 exact files /3,216,352 bytes under `data/tic214311908_external_comparison/source_snapshot`, with pinned URLs and Git blob/SHA-256 identities in `source_manifest.json`. The ZIP exceeded the transfer limit; the user-supplied public directory was used instead, with no ZIP-byte identity claim. Cached selected responses remain outside Git at `/workspace/tess-research-data/tic214311908_external_review`. Reproduce compact checks with:

```sh
cd /workspace/Space_ChatGPT
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/compare_tic214311908.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python candidates/tic214311908_external_comparison/shape_checks.py
```

The first verifies pinned sources, geometry, catalogue table selection and all 1500 prior-draw summaries. `--restore-sources` retrieves only missing expected public files without executing them. The second fits our existing saved QLP windows and is not an external-code execution. Preserve original outputs and source hashes before deliberate regeneration. Two delegated reviews reached their execution quota after supplying numerical results; root completed the audits and preservation locally. No paid calls or external messages/submissions. Next decisive science remains resolved in-event multi-band photometry plus companion imaging and spectroscopy/RVs, not another attractive fold.

## Latest targeted review — TIC 214311908, 17 October 2026

The user requested an independent investigation and withheld the other model's findings. None was accessed. **A promising recurrent 5.6915-day transit/eclipsing signal survives; its host and planetary origin remain unresolved.** Read [the full evidence](candidates/tic214311908/README.md), [measurements](candidates/tic214311908/measurements.csv), [independent analysis](candidates/tic214311908/independent_variability/README.md), [source assignment](candidates/tic214311908/source_localization/README.md) and [TARS recognition audit](candidates/tic214311908/tars_recognition/README.md).

Original QLP Sectors 12/39/66/93 contain 13 positive local event windows. A frozen discovery ephemeris then predicts eight positive sampled windows in previously uninspected native FFIs from Sectors 103/104; another S103 prediction falls in a gap. The 21 positive windows are not all individually significant. TARS and ELEANOR reductions support the same sequence and preserve one additional early window, but share the original photons. No QLP 103/104 timing or flux was used because of documented product issues.

Independent discovery-only timing gives **P = 5.691513806 ± 0.000008821 days; T0 = BTJD 1632.611612 ± 0.002312 TDB**. Its worker never accessed held-out data, but saved the timing snapshot after root retrieved those data. The predictive test used the earlier preregistered BLS timing, P 5.691497929 days and T0 1632.614684206 BTJD; see [the chronology audit](candidates/tic214311908/holdout/timeline_audit.json). Keep the full covariance and systematic floor in `candidates/tic214311908/independent_variability/ephemeris_frozen_pre_holdout.json`, whose original filename is preserved. Convert barycentric TDB for a real observer before scheduling UTC. A nominal-target circular fit gives approximately 1.65 hours and a conditional radius of 5.1 Earth radii. Tested equal-star, doubled-period models fit worse by Δχ²≈15.7; this excludes neither other stellar configurations nor a background host. No mass, FPP or calibrated global FAP is measured.

Host ambiguity is decisive. Faint TIC 1515709327 at approximately 3.15 arcseconds has about 2.9% of the nominal target T flux; an approximately 50% eclipse can mimic the QLP deficit. The nearly comoving TIC 214311910 at 11.7 arcseconds is another possible host. Propagated Gaia/native pixel fits cannot discriminate these sources from the target. QLP already adds the expected target TIC flux to difference photometry: do **not** generically multiply the depth by 1+contamination ratio. A separate weak 0.15804-day feature matches a known neighboring binary; its changing phase relative to the deep events rules out a simple alias explanation. Finite catalogue/literature checks found no prior planet designation or matching coherent transit ephemeris, but public TARS headers already contained weaker near-5.7-day sinusoidal peaks. Novelty is unproven.

This instance began at old snapshot `2864163`, with no old raw cache or astronomy packages. The clean checkout was fast-forwarded to saved main `ff9216f`; the existing setup script restored pinned dependencies and passed import/pip checks and real FITS/BLS work. Reusable install/start instructions were saved as a configuration draft, not published. Previous broad/multisector live caches remain absent: restore their committed archives rather than rerunning finished surveys. The new cache is `/workspace/tess-research-data/tic214311908/`. All 13 original scientific TESS FITS are additionally archived for exact offline restoration. Redundant ZIPs and the 136 MB parent TARS catalogue remain cached with checksums/URLs. No originals were deleted. OpenRouter research spending remains $0; no external messages, submissions or observing requests.

Verify before regenerating artifacts:

```sh
cd /workspace/Space_ChatGPT
export PATH=/workspace/openrouter-check-venv/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/workspace/tess-research-data/matplotlib
python scripts/tic214311908_checkpoint.py verify
python scripts/tic214311908_checkpoint.py restore --destination /tmp/tic214311908_new_offline
# If the original cloud cache is absent, restore directly to its missing paths:
python scripts/tic214311908_checkpoint.py restore --destination /workspace/tess-research-data/tic214311908
python scripts/tic214311908_retrieve.py --download-qlp
```

Restore refuses existing mismatching files. The scientific FITS and compact queries restore without network. After a lost-cache restore, use `verify --science-only` to check all scientific FITS and artifacts without demanding redundant ZIPs or the unneeded parent TARS catalogue. To reproduce current derived checks deliberately:

```sh
python scripts/tic214311908_screen.py
python scripts/tic214311908_screen.py --pool
python scripts/tic214311908_vet.py
python scripts/tic214311908_pixels.py
python scripts/tic214311908_reductions.py
python scripts/tic214311908_holdout.py
python scripts/tic214311908_holdout_diagnostics.py
python scripts/tic214311908_summary.py
```

Read each independent method before rerunning: preserve the original assessment, main preregistration and separate discovery-only timing snapshot. The latter's original filename says pre-holdout; the chronology audit qualifies that wording. The stellar comparison retains exactly the same 13 events/1961 samples as the original physical fit. Earlier software attempts and expected hashes are retained; do not overwrite them to accept a mismatch. The checkpoint establishes byte identity and restoration, not astrophysical validity.

**Next decisive science:** resolved event photometry measuring the target, TIC 1515709327 and TIC 214311910 simultaneously, supported by multicolour fluxes and accurate PRF modelling. Target-only RVs could miss a background binary. Existing TESS resolution cannot securely assign a 3-arcsecond source. No observing request or submission was sent. This one-star review does not modify the frozen survey counts below.

## Latest checkpoint — completed six-sector batches, unresolved TIC 229146531

**16 October 2026 UTC.** Actually period-searched **5,058 target-sector observations /4,391 distinct TICs** across sampled Sectors83/86/95/96/98/105. 5,059selected observations, one genuinely unusable S83curve. Producer-error searches3819; separate empirical-error searches1239. First batch4095observations /3612stars; second963/850;71shared TICs counted once. Exactly3998TICs were outside the frozen159505-star Sector1survey. Do not add repeated pools, controls, metadata inventories or follow-up observations as new survey stars.

All planned pools completed:471producer-only plus12mixed first-batch repeated stars,113second-batch repeated stars and71cross-batch all-observation searches. These sets overlap. No survey, pool or follow-up process is left running. Read [the current scientific report](MULTISECTOR_RESEARCH_REPORT.md), [searchable deduplicated table](results/multisector/all_observations.csv), [accounting/index](results/multisector/campaign_summary.json), [ranked assessments](results/multisector/ranked_assessments.csv) and [last follow-up report](reports/multisector/batch02_followup.md).

**No credible new planet or previously unknown star established.** [TIC229146531](candidates/multisector/tic229146531/README.md) is the highest-priority unresolved signal: P≈0.935024d in six sectors, stable native target-area deficit, but dark-companion and doubled-period grazing stellar models have almost identical fits. Already a Gaia321.389d astrometric binary; that orbit does not explain the short-period signal. Odd/even differences reverse between sectors. No measured inner companion mass, reliable radius or novelty established. Initial coarse-period historical non-recovery was invalid because phase uncertainty accumulated. Do not revive that rejection.

Known-neighbor contamination, published SB1orbits, variability aliases and native background/image-motion tests weaken other leads. New final pooled checks added14stars /30observations, then external sectors, native cubes and a reference-image/gradient positive control. TIC183373155is background/redistribution affected; its S29predicted event is in a gap, not a null. TIC71274048has no convincing native source-localized deficit. The low-S/N Gaussian offset for293268396and ill-conditioned high-order baseline extrapolations are invalid diagnostics. Frozen full helper crashes for zero events have explicit coverage/local-statistic fallbacks; do not silently replace missing depths with zeros.

### Environment and verification

Repository `/workspace/Space_ChatGPT`; Python `/workspace/openrouter-check-venv/bin/python`;5CPUs, approximately23GBfree at the last check. Raw/static/prepared data, live SQLite and preserved bulk figures remain under `/workspace/tess-research-data/multisector/`. The overlay is not independently verified durable storage. No OpenRouter requests were made: project research spend remains$0. No external messages/submissions/observing requests. Worker usage limits were reached; root continued locally through additional astrophysical/pixel checks and preservation.

```sh
cd /workspace/Space_ChatGPT
export PATH=/workspace/openrouter-check-venv/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/workspace/tess-research-data/matplotlib
python scripts/multisector_checkpoint.py verify
python scripts/multisector_inventory_checkpoint.py verify
python scripts/latest_sector_probe.py verify
python scripts/multisector_supplemental_checkpoint.py verify
python scripts/multisector_candidate_checkpoint.py verify
python scripts/restore_multisector_raw.py --tic 229146531 --max-downloads 0
python scripts/summarize_multisector.py
```

The original core source hashes remain: `src/broad_vet.py`23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333; `src/fast_screen.py`26eede8aaeb9fa2a624b73a5c3966366f77c4cc73738153eba7e46d13eb10307. New scientific drivers keep independent identities; method snapshots are retained. Never replace expected hashes merely to accept a mismatch. The main backup-WAL bug was caught, fixed and meaningfully regression-tested; the failed archive is retained in cache and identified by `data/multisector/checkpoint/wal_failure_record.json`.

### If a future environment loses its cache

The main archive restores4096producer records; supplemental archives restore551first-batch recovery records,963second-batch producer records,689second-batch recovery records, all667completed pooled-star JSON files, and exact small follow-up responses. Both metadata query archives preserve exact queries. TIC229146531's7original FITSincluding the native cube restore exactly without network. Use a fresh destination for a dry run; these commands refuse different existing files:

```sh
python scripts/multisector_checkpoint.py restore --destination /tmp/new_multisector_main.sqlite
python scripts/multisector_supplemental_checkpoint.py restore --destination /tmp/new_multisector_supplement --archive-only
python scripts/multisector_candidate_checkpoint.py restore --destination /tmp/new_multisector_candidate
python scripts/multisector_inventory_checkpoint.py restore-queries
python scripts/latest_sector_probe.py restore-queries
```

For actual cache reconstruction, use `/workspace/tess-research-data/multisector/survey.sqlite` as the missing main destination and `/workspace/tess-research-data/multisector` as the supplemental/candidate destination. Restore completed results first; do not repeat thousands of searches. `python scripts/restore_multisector_raw.py --tic TIC_ID --max-downloads 20` retrieves missing static QLP FITSby recorded URL/hash; omit `--tic` only for a deliberate bounded broader restore. Existing mismatches stop restoration; mismatching new responses are retained separately. Dynamic TESSCut regeneration may change container/header bytes, so important229originals are archived directly and other cutouts have recorded science hashes/FFIidentities. Prepared NPZand bulk figures can be regenerated with the corresponding frozen code; source method changes must be separately versioned. Nothing was deleted simply to save repository space.

### Exact next scientific work

1. Independently reduce/classify TIC229146531. Reproduce its existing physical/parity/overview tests with `scripts/tic229146531_physical.py --models --timing`, `scripts/tic229146531_parity.py`, `scripts/tic229146531_overview.py`; these intentionally regenerate derived files. Test joint variability/eclipse models and both0.935/1.870d hypotheses; do not interpret a grazing radius or conditional chi-square difference as a mass/probability. Seek independent reduced archive data and expand alias/literature recognition before any novelty claim. Resolved multicolour photometry and RVs would discriminate physical hypotheses; no observing request or catalogue submission has been made.
2. Improve correlated-noise calibration and300–500ppm recovery in a new versioned method. Existing288injections/32surrogates show poor300ppm sensitivity and no globalFAP. Preserve frozen controls/results. Avoid merely lowering significance cuts; quantify completeness and structured false alarms on more noise backgrounds, especially QLP missing-error products.
3. Define another pre-flux sector sample only after selecting the method/control criteria, using `scripts/latest_sector_probe.py` and `scripts/run_multisector_batch.py` with a fresh batch name. Timing exclusions99–104and corrected74–79products remain mandatory under current evidence. The last106/107inventory requests were404, not searched sectors. The current6sectors are partial samples; many targets and sectors remain unexamined.
4. Preserve raw datasets, update the ledger/log/docs and commit meaningful progress. Never call another completed scan a newly searched star. Reuse the current immutable selection and named batch machinery; do not write two processes into a live SQLite database.

## Latest targeted comparison — TIC 4206066

15 October 2026, America/Halifax. Read the user's supplied Rabtsevich draft and checked exact TIC membership in the old complete Sector 1 parent inventory, selected targets, exported targets/signals and read-only survey database. All absent. The prior Reddit review mentioned this star but no previous photometric analysis was executed. The old 159,505-star count is unchanged.

Retrieved and preserved the paper's Sector 6 TESS-SPOC, Sector 32 SPOC and Sector 98 QLP observations. Unmodified screening with the actual original broad configuration finds near-3.18-day peaks at 3.179552/3.184918/3.181348 days. Local S/N is 5.37/8.83/6.32; the actual prioritizer retains only Sector 32. Sector 98 also retains a near-11.13-day peak, without establishing three independently clean transits or validating the tentative planet. Current QLP DET_FLUX needs a file-interface adaptation from old KSPSAP_FLUX; search functions and the old survey are unchanged.

Read [the comparison](candidates/tic4206066/README.md), [complete results](candidates/tic4206066/recovery.json) and [diagnostics](candidates/tic4206066/recovery.png). Reproduce with `/workspace/openrouter-check-venv/bin/python scripts/tic4206066_recovery.py`, adding `--download` only to restore missing originals by immutable hashes. Raw FITS and PDF remain in `/workspace/tess-research-data/tic4206066/`; compact curves and outputs are tracked. This control does not reproduce pixel localization, physical transit fitting or FPP tests. The paper explicitly claims a candidate, not statistical validation. No paid API calls. Return to the broader discovery checkpoint below after this requested comparison; no additional discovery claim follows.

## Latest fresh search for uncatalogued stars and companions

Evening 14 October 2026, America/Halifax; archive manifests use UTC 15 October. The user requested a new search beyond the old two candidates and then asked to continue. **No genuinely uncatalogued star is established.** Four complementary methods were actually executed and promising signals actively falsified. The investigation remains open; no claim that all public astronomical data are exhausted. Read [NEW_STAR_SEARCH_REPORT.md](NEW_STAR_SEARCH_REPORT.md), [search counts](results/new_stars/search_summary.json) and [ranked leads](results/new_stars/ranked_leads.csv).

- Direct imaging:12-field pilot plus 64-field stratified expansion;1035 recurring peaks, 495 Gaia-missing but deeper-catalogued.27 residual leads reduced to one, then rejected as inside all g/r/z native saturation masks (value 18462). One separate faint, broadened pilot peak at RA 216.428278/Dec 23.612155 remains unverified after stationary/common-motion PS1, deep-stack and SDSS checks. Not a star candidate.3832 static injection trials quantify only the detection step; at 4 arcsec, magnitude 21 recovery is 2.85%, so the negative result does not exclude close faint companions.
- Gaia:4103 fixed primaries, 1086 neighbor rows, 11 unique screened associations.5 were already recognized. Of 6 remaining, pairs 1/3 have cleaner quality and independent 1993/1996 paired images. Gaia IDs 4977459112829640832/4977459112829654400 (4.034 arcsec) and 4923799509418222080/4923799509418221952 (4.763 arcsec). Existing source positions, possible new associations only. No binding, masses or priority proven. Others have compromised photometry/dynamics. Catalogues' selection gaps mean null matches do not establish novelty.12309 shifted-primary trials yield zero passes, with finite-count/clustering caveats.
- Infrared:1946 selected CatWISE rows in 10 successful fixed fields, two failed fields;250 Gaia-unmatched and 10 matched controls received old-epoch tests.27 independently consistent no-Gaia cases are all already TIC STAR entries.28 optical tests leave one apparent mover, then SDSS falsifies its proposed optical counterpart as a stationary galaxy (0.126 arcsec over 6.87 years versus 2.64 arcsec predicted). Different blended movers remain unexcluded; no new star established.
- Spectroscopy:3430 repeated GALAH stars/13176 quality visits.39/40 strongest cases already flagged binaries; last strong case has reduction problems. A 1799-star strict subset leaves two weaker sparse-RV leads at Gaia 4690783166942834816/TIC 51939231 and 5437472031746145536/TIC 24731098. Same-field controls, four observed-flux spectra, exact rest-frame handling, known-catalogue checks and three QLP curves were examined. Wavelength/mask dependence prevents precise RV/orbit claims; no component mass or new star. The first is already a 20.766-day ROT variable.

All original observations and 218 intermediate image plots were retained.186 MB bulk plots moved intact outside ordinary Git; four diagnostics remain tracked. `data/new_stars/imaging/plot_preservation.json` proves identities. `/workspace/tess-research-data/new_stars/` has the cache. Small public-query snapshots are additionally archived in `data/new_stars/checkpoint/` for exact offline restoration. The overlay is not independently verified durable storage. No paid APIs or external messages/submissions/observing requests. Parallel workers hit usage limits; root continued locally rather than abandoning analysis.

Verified checkpoint: 10 methods, 479 scientific artifacts and all 5,955 cached files (377.9 MB) match their recorded hashes. Offline query restoration independently reproduced all 5,466 archived files exactly. All 10 new Python modules compile; Gaia derivation/controls and all five infrared stages reproduced from retained inputs. This verifies reproducibility, not discovery claims.

Reproduce without rerunning completed broad TESS searches:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/workspace/tess-research-data/matplotlib
python scripts/new_star_checkpoint.py verify
# If a future cloud task loses raw query responses:
python scripts/new_star_checkpoint.py restore-queries
python scripts/new_star_checkpoint.py download --max-downloads 20
python src/new_companion_search.py
python src/infrared_star_scout.py fetch
python src/infrared_star_scout.py gaia
python src/infrared_star_scout.py counterparts
python src/infrared_star_scout.py novelty
python src/infrared_star_scout.py optical
python src/spectroscopic_companion_scout.py
python scripts/rv_lead_recognition.py
python src/rv_spectral_validation.py
python src/rv_photometric_veto.py
python scripts/summarize_new_stars.py
```

Image runs use `python scripts/run_uncatalogued_imaging.py --expanded --run-name expanded64_v4`, which preserves bulk plots outside Git. Frozen version snapshots/earlier outputs are in the image data directories. After deliberate regenerated artifact changes, `python scripts/new_star_checkpoint.py freeze` refreshes a verified checkpoint; never change expected raw hashes to accept a mismatch.

Exact next science: independently review pairs 1/3's literature/other catalogues, then obtain differential RVs or stronger astrometry to test physical association. For the faint image peak, split original DECam exposures by night/band, model the primary's motion and require persistent PSF-like detections; existing coadds/PS1 are insufficient. For RV leads, seek genuinely reduced/independent spectral epochs and resolve line-profile/calibration discrepancies before fitting or claiming a stellar companion. The available ESO/RAVE/APOGEE cone queries are retained and bounded, not proof that no spectrum exists anywhere. Two timed-out infrared fields can be attempted under new response names without overwriting old errors. Do not count reruns or metadata as new stars.

## Latest two-system deep validation

Updated 2026-10-14. **Both eclipse signals survive further falsification; no new star is confirmed. TIC273879360 was already publicly machine-flagged for eclipse-like variability.** The prior description of Kostov2025 Table2 as a generic target inventory was incorrect; it is a published NN eclipse-feature collection. This correction changes novelty assessment, not the measured source localization/eclipses.

- TIC273879360: independently re-extracted four native pixel sectors, 22 accepted local event windows and 6234 cadences. Shared two-luminous-sphere fit gives P7.586890837d and conditional e0.32145; baseline/aperture/limb variants retain approximately e0.32. Archived ASAS-SN independently shows nine primary events, eight before TESS. Combined-source catalog properties do not measure component masses or unique radii. See [deep validation](STELLAR_VALIDATION_REPORT.md).
- Giant TIC260352037 field: 16 genuinely later QLP sectors bring total to41 across7.45years. Old timing predicts both families; raw central depths1696±27/784±28ppm and conditional timing P0.880788521±0.000000330d. Four cubes localize the blend but empirical-image nuisance fits fail secure source attribution. Independent Gaia spectroscopy confirms giant radius≈9.45Rsun. Neither nearby-source identity nor gravitational binding is established.
- New context includes42 source-hashed snapshots, positive spectroscopy/wide-binary query controls, contaminated faint-neighbor colors, and explicit GaiaRV/MSC interpretation. Gaia assumed-binary MSC fits do not independently recognize binaries. Ground HJD values use field centers and approximately200s accuracy; do not merge them into a precision BJD(TDB) ephemeris without time conversion.

Resume this checkpoint with the prepared Python and cache variables below. New background analysis commands are `python src/validate_background_binary.py`, `python src/background_binary_controls.py`, and `python src/background_binary_alias_tests.py`. They regenerate the four-source raw fit, sensitivity/block-bootstrap controls, and missing-eclipse/early-period diagnostics. First restore missing native cubes from their recorded manifests using the older source commands below; never overwrite a mismatching raw file.

`python scripts/stellar_background_verify.py` verifies the four raw cubes, method/output identities and numerically reproduces all22event/6234cadence fit arrays. After intentional regeneration, use `--write-manifest` once to checkpoint newly verified artifact bytes. The original early-only holdout fit did not converge; its large penalty is not a physical rejection. The final broadened profile converges and includes the full-data period inside its conditional delta-chi-square3.84 crossings. Missing half-period eclipse controls are negative. Current and superseded failed/profile artifacts are preserved.

The current two-system searchable table is `results/stellar_systems/deep_validation_measurements.csv`; ranked assessments now reference the new fits/dossiers. Their old measurements remain preserved in historical reports and Git. Verified background checkpoint:4raw cubes,4methods,13artifact hashes plus exact fitted-array reproduction; giant checkpoint:106code/raw/output hashes. Source modules compile and `pip check` passes. The frozen broad-vetting SHA remains unchanged.

Giant reproduction/retrieval and independent manifest checks are in [its new report](reports/stellar_systems/tic260352037_validation.md). `python data/stellar_systems/tic260352037_validation/verify_manifest.py` checks106 code/raw/output hashes. `python scripts/stellar_deep_context.py` verifies retained snapshots and derives current properties; `python scripts/stellar_ground_photometry.py` reproduces the archived ground analysis. Add their documented retrieval flags only for missing observations. The pinned environment now includes pyarrow26.0.0 for preserved public ASAS-SN index data.

The next useful science is resolved photometry of the giant and its4–8arcsecond neighbors during both eclipse families, and resolved multiband photometry plus spectra for TIC273879360. Such new observations are outside currently available public time series; no observing requests/submissions were sent. The earlier64-target reserve can still support an independent discovery search, but was not part of this two-candidate follow-up.

## Latest stellar-system investigation

Updated 2026-10-14. Completed an EB-oriented reclassification of the preserved survey, selecting 1,253 stars and actually remeasuring a fixed 192-star pilot. All 192 raw hashes checked; 153 repeated raw leads, 40 half-phase secondary warnings and 48 parity warnings. Expanded recognition falsified the first promising scout batch: known EB/SB2 hosts and known neighbor contamination predominate. This phase does not add those catalog rows to the original 159,505 analyzed-star count.

**Strong source: TIC 273879360 / Gaia 4692910966756908288**, a background source 65.59″ from the original target TIC 273879358. Both eclipse phases localize there in S27/29/67/95. Core depths approximately 21–27% / 7–9%; earlier phenomenological period7.586891287d, secondary phase0.69179. No vetted prior binary classification found in the recorded finite checks, including both complete vetted TESS Ten Thousand tables and GALAH SB2 catalog. **Later correction:** its unvetted Kostov2025 row is a public NN eclipse-feature flag, not a generic inventory entry; first recognition is excluded. S68 is off science silicon, retained and excluded.

**Second candidate: unresolved short-period system near giant TIC 260352037.** 25 sectors, early-sector timing predicts dips in23later sectors; P0.880788466d, raw primary/secondary972±26/445±38ppm. Native pixels near giant image, exact source unassigned; close common-motion neighbor737466840 is a hypothesis. The other two giant/source tests failed heldout coherence. TIC197764088's apparent new EB is independently localized to known neighbor197764089, retained as a negative result. See [STELLAR_RESEARCH_REPORT.md](STELLAR_RESEARCH_REPORT.md).

Resume the stellar phase without repeating the192 completed curves:

The commands in this older section reproduce frozen historical analyses. Their legacy `unvetted_inventory` field names and old source-fit limitation text predate the Table2 correction; interpret them with the latest deep-context report. Current physical validation lives in separate folders so reproducing old numerical evidence does not replace the new assessment.

1. Verify `results/stellar_systems/scout_manifest.json` and source hash of `src/stellar_scout.py`. Exact completed reproduction is `python src/stellar_scout.py --limit 192 --deep-limit 48`; cached results skip unchanged targets. Do not change strata/defaults and overwrite the fixed pilot accidentally.
2. Main background source reproduction: `python src/stellar_pixels.py data/stellar_systems/tic273879358_large_config.json --sectors 27 29 --bootstrap 300`, then `python src/stellar_pixels.py data/stellar_systems/tic273879360_config.json --sectors 67 68 95 --bootstrap 200`, then `python src/stellar_background273.py`. The second command intentionally returns1 after recording S68 exclusion and completing usable sectors. Set the environment/cache variables below first. Keep all raw failed products; never accept prescan/overscan as a missing transit.
3. Current source prior checks: `python scripts/stellar_prior_audit.py`; `python scripts/stellar_scout_novelty.py --full-crossmatch --report`. Reuses verified snapshots. Fetch flags are documented in each script; full 36.4MB combined inventory is retained outsideGit with its manifest. Read `reports/stellar_systems/prior_audit_supplement.md` and actual vetted catalog rows, not only web-page existence or unvetted target inventories.
4. Strong source follow-up now needs resolved ground primary/secondary photometry on TIC273879360 and spectra; giant source needs finer spatial separation. No external messages/submissions are authorized. Public mass data are insufficient; do not invent masses or claim confirmed new stars.
5. The next64 reserve entries in `results/stellar_systems/scout_next_pilot.csv` already have raw scout measurements but incomplete live recognition. First twelve include370004684,31779604,394731519,277314291,33665798,150431384,52059035,39049265,238832404,140052356,29348301,425938488. Perform current VSX/ASAS-SN/Gaia/SIMBAD/spectroscopic recognition before new independent-sector/pixel work. This metadata task was not completed after parallel-agent usage limits; source claims remain unpromoted. Do not count this queue as64additional analyzed curves.

All new scientific data remain retained; reports, event measurements, source images, independent negative controls and retrieval manifests are tracked. No paid API calls. The earlier exoplanet and targeted-claim checkpoints below remain preserved.

Validation: `python scripts/stellar_restore_lightcurves.py` verified all60 retained QLP files for the giant tests and known-neighbor follow-up; add `--download` only to recover missing files. `python scripts/stellar_verify.py` verifies335 catalog snapshots, ten raw pixel files (including excluded/pilot evidence), six current pixel reports, all192 scout method identities and the background-source bootstrap attribution. Eleven existing utility tests passed. The latter verification excludes mass/novelty claims by design.

Updated 2026-10-14. **159,505 distinct TIC targets actually period-searched, 495 insufficient-data cases, zero errors. All 389 strict shortlisted scan/peak signals on 357 TICs received detailed v3 checks, including 32 alternate peaks. No credible new exoplanet candidate is established.** The historical 19-star pilot is not added to these counts and has no TIC overlap with the broad sample.

## Latest targeted request — TIC 71316387

Completed the user's alternating-eclipse claim audit: five retained QLP sectors 1/28/68/95/105, 134 primary depth fits and 99 accepted timings, plus an independent 136-event implementation and native sector 28 pixel localization. Both eclipse types lie near the target; persistent unequal depths favor a roughly 1.30367-day eclipsing system. The original precise 5293/2487 ppm/19.6σ calculation is not certified; the main common contrast is 2466 ± 286 ppm with conservative sector covariance. Current RAVEN entry is a PC, not validated; original RAVEN summary used S1/28 and includes odd/even diagnostics. No new planet or stellar-mass confirmation is claimed.

Read [the complete audit](candidates/tic71316387/README.md). Reproduce with `python scripts/tic71316387_retrieve.py --download`, `python src/validate_71316387.py`, and the independent catalog/pixel/event scripts listed there. No broad-survey frozen source or stored classification was changed. All raw observations/pixels/paper sources remain retained; no paid calls. Next useful steps for this target are spectra/resolved follow-up and obtaining the original quoted analysis to compare masks/noise definitions. No external communication is authorized. The broader research resumption steps below remain available.

## Secured scientific conclusions

- **TIC 269994642:** real approximately 2.7% eclipses independently recur in 2018/2025 pixels near the target. Photometric period 1.746372 days; a 3.492744-day grazing binary fits as well as a dark transit. Unequal depths, conservative parity about 3.33 conditional standard errors, and Gaia RUWE2.08 favor binary testing. Radius profile reaches a stellar-sized upper boundary. No mass or planetary classification. Sector27 coverage was prescan/collateral data and is excluded, with the failed experiment preserved.
- **TIC 238202207:** 17 sectors/140 event timings, severe Gaia unresolved-source diagnostics, alternating depths positive in16/17 sectors (about4.39 conservative conditional standard errors). A 4.643651-day unequal binary fits as well as a planetary model at apparent period2.321826days. Catalog radius disagreement materially changes inferred size. Binary favored, no accepted planet.
- **TIC 260820871:** shallow dips repeat across S1/27/28 and alternate strongly across93 fitted events. Independent S28 differences localize both eclipse phases to the faint known VSX eclipsing neighbor, rather than the hot target. Published neighbor period is approximately harmonic, not an asserted exact match. Physical diluted-binary geometry is viable. This is contamination, not a new planet.
- **TICs 370247818,141980042,300814564,273879358,270312735:** independent sectors reveal significant secondary eclipses; stellar sources favored. 370 has physical/pixel evidence;300 has23 usable sectors and224 event timings. Unusable original products remain retained.
- **TICs 139148478/431481555:** weak later support plus contamination/recognized pulsation respectively; no credible transit interpretation. **290028406:** later raw local support3.24, V-shaped fit and tentative half-phase feature; unresolved lower priority. Other clean-warning rows have missing stellar parameters or severe crowding/broad duration. Warning-free does not mean planet.
- Recent-sector track:17 actual SPOC curves on5 stars. Three variability/baseline aliases, a published CPV, and already-reported candidate TOI2417.11 (confirmed source localization, unconfirmed planetary nature). Public846TCE rows/623stars are metadata counts, not our analyzed light curves.

Read `BROAD_RESEARCH_REPORT.md`, `candidates/broad_s1/FOLLOWUP_REVIEW.md`, `candidates/broad_s1/all_signal_evidence.csv`, and physical/source dossiers. `results/broad/focused_assessments.csv` summarizes focused interpretations. All scores/parity tests are conditional diagnostics, not survey-wide significance. Absence from finite catalog/literature searches is not proof of novelty.

## Resume without repeating completed searches

1. Use existing checkout `/workspace/Space_ChatGPT`, interpreter `/workspace/openrouter-check-venv/bin/python`; no worktree. Set `TESS_DATA_DIR=/workspace/tess-research-data`, `MPLCONFIGDIR=/workspace/tess-research-data/matplotlib`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`. Check active processes and Git status first. No app service is needed.
2. Full SQLite state is `/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite`. If lost, first `python src/survey_checkpoint.py verify`, then `python src/survey_checkpoint.py restore --db /workspace/tess-research-data/broad/qlp_s0001/survey.sqlite` to an unused destination. Four committed chunks total~78.1MB, each<20MB. Actual160000-row restore/integrity/exact-content comparison passed; see `results/broad/checkpoint_roundtrip.json`. Restore does not download observations.
3. Verify/retrieve raw broad files with `python src/restore_broad.py --limit 20` initially. Downloads require `--download --max-downloads N`; preserve mismatches. Original pilot: `python src/retrieve_manifest.py`. Important follow-ups have separate manifests under `data/broad_followup_queries`, `data/broad_pixel_*.json`, and `data/tce107*`. Preserve all raw observations, intermediates and failed experiments.
4. The full survey has no pending targets. `src/broad_survey.py` source scanner hash is `26eede8aaeb9fa2a624b73a5c3966366f77c4cc73738153eba7e46d13eb10307`; existing results skip completed stars. Do not rerun full searches when the committed checkpoint suffices. `python src/summarize_broad.py` reconstructs the audit/overview from stable exports.
5. Detailed frozen v3 hash is `23c89080f4f10ef59b79d2896c86edb98dafeecc14081feb2c2a5c69fb37f333`. Completed main357 + alternate32 checkpoints are versioned; histories preserve earlier187-star batches. `python candidates/broad_s1/summarize_review.py` rebuilds exact389-row ledger and verifies hashes. Re-running complete review is unnecessary.
6. Source renderer contained a target-specific contamination sentence copied to other targets. Corrected prose/text hashes are preserved. After any future use of frozen `broad_vet.py --multisector`, run `python candidates/broad_s1/render_multisector_notes.py`; numeric arrays/methods remain unchanged. A gap-boundary corner case can leave a secondary window untested: check actual event coverage, never interpret zero/unavailable as absence.
7. Current-TIC table: `data/broad_catalogs/current_stars_20261014/current_stellar.csv`,3410exactrows,1541evolved/giantcontexts,270missingradii. No old-QLP fallback for missing current values. `shortlist_uncapped.csv` and its manifest are frozen input to the389-signal review. Preserve snapshots before changing selection or method.

## Exact next scientific steps

1. The final TIC260820871 pixel dossier is complete: both eclipse phases favor the known neighbor under four baseline/template choices and all 500 bootstrap draws. Use its preserved raw retrieval and reproduction commands only if verification/resumption is needed. Do not promote the bright target's inferred radius onto the known faint eclipsing neighbor.
2. Strongest ambiguous systems269994642 and238202207 need actual spectra/resolved follow-up to classify companions; available medians do not measure masses. Review finite recognition scope before any public novelty claim. No messaging to observers is authorized.
3. Remaining lower-priority source-warning reserves from the32-signal audit are TIC262608189 (originalP2.134621) and232074299 (originalP.966920). Retrieve independently available sectors before further ranking; assess current source identity, actual raw recurrence and both true/doubled-period hypotheses. Other warning/control groups were not all deeply followed.
4. For a new search phase, preselect quiet stars with current stellar/source validation and multiple genuinely usable sectors. Add several detrending scales and more independent peaks; measure real-noise injection recovery for weak/long periods and native short durations. Use saved isolated-event data as a separate longer-period track. Record selection before outcomes, fix the gap-boundary tolerance in a separately versioned method, validate doubled-period binaries, and enforce pixel detector/science-column checks. The coarse negative screen does not exclude further planets.

## Preservation, environment and budget

All160000rawQLP products are retained (exact12,902,400,000bytes), with follow-up pixels/light curves and SQLite outside ordinary Git. QLP intrinsic checksum cards are absent, explicitly reported; SHA256 identifies exact bytes, and products with checksum cards are verified. Nothing important was deleted. Eleven utility tests, scientific controls, model integration/binary-equivalence checks and database round-trip passed. Pinned installation is `bash scripts/install_environment.sh`.

This overlay has no independently verified durable volume. Git preserves full result state and retrieval instructions; environment publication can snapshot prepared raw files. Configuration saving alone does not publish, and live processes never persist reliably. Fresh-cloud-task restoration remains untested even though same-machine database restoration passed.

Research-phase OpenRouter usage remains zero; the earlier connectivity charge was not recorded. Enforce the cumulative$10cap with a conservative reserve if future calls become useful; no paid calls are needed for this checkpoint. Never display/log/commit secrets. Commit meaningful progress and push main as authorized. The historical pilot report and next steps in `reports/pilot_progress_20261013.md` are superseded by this broad checkpoint.
