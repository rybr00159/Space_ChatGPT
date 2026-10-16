# Six-sector exoplanet investigation

Updated **16 October 2026 UTC**. The defined batches and every planned blind repeated-star pool are complete. **No credible new planet or previously unknown star is established.** TIC229146531 is a real recurrent eclipse/transit signal in an already known astrometric binary, but its short-period companion class and publication novelty remain unresolved. The investigation remains open; this is a reproducible scientific checkpoint, not an exhaustive TESS search.

## Actual scope and results

**5,058 target-sector observations of 4,391 distinct TICs were period-searched**, from 5,059 selected observations. One Sector83 light curve remained unusable. This phase samples six sectors; it does not constitute a complete census of any one sector.

| Sector | Selected | Producer-error searches | Separate empirical-error searches | Unusable |
|---|---:|---:|---:|---:|
| 83 | 1024 | 1023 | 0 | 1 |
| 86 | 1024 | 1024 | 0 | 0 |
| 95 | 1024 | 566 | 458 | 0 |
| 98 | 1024 | 932 | 92 | 0 |
| 96 | 456 | 155 | 301 | 0 |
| 105 | 507 | 119 | 388 | 0 |
| **Total** | **5059** | **3819** | **1239** | **1** |

The first batch actually searches 4,095 observations /3,612 unique stars. The brighter/cooler second batch searches 963 observations /850 unique stars. Their 71 shared TICs are counted once in the 4,391-star total. Exact ID deduplication finds **3,998 stars outside the older 159,505-star Sector1 period-searched sample**, with 393 in common. The union of these two defined surveys is 163,503 TICs; targeted follow-up/control stars are not added to that union.

The [searchable observation table](results/multisector/all_observations.csv) has one row per selected target-sector, including failed original searches and successful uncertainty assumptions. [Campaign accounting](results/multisector/campaign_summary.json) indexes all eight searchable peak tables and their hashes. There are 80,928 single-sector peak rows and 10,752 pooled rows; these are **correlated trials, not planets or distinct independent detections**.

All **483** usable repeated stars in batch1 (471 producer-only, 12 mixed-error) and all **113** repeated stars in batch2 were blindly pooled, regardless of individual-sector scores. The **71 cross-batch TICs** also received all-successful-observation pools. These sets overlap; adding them as new stars is incorrect. All planned pools completed, without failed pooled stars.

## Strongest outcomes and active falsification

The [ranked assessment ledger](results/multisector/ranked_assessments.csv) records outcome classes and evidence links. Rank is follow-up priority, not planetary probability.

| Target | Scientific outcome |
|---|---|
| **TIC229146531** | Approximately **0.935024-day** dips recur in **S3/29/30/69/96/105**, with raw local S/N13–35. Native S96 pixels retain approximately 3000ppm target-area loss under sky, baseline, quality and image-gradient changes. However, dark-companion and equal-star grazing **1.870047-day** binary shapes fit almost identically (raw χ²3848.32/3847.63 for4183points). Companion radius/density/dilution are severely degenerate. Already a published Gaia **321.389-day astrometric binary**. Alternating-depth differences reverse sign between sectors and do not establish an inner binary either. **Unresolved transit/eclipsing signal; no planet or new-star claim.** [Full dossier](candidates/multisector/tic229146531/README.md). |
| TIC368248545 | The attractive recurring1.124-day peak accompanies0.16054-day modulation matching half a known neighbor's0.321179-day EW period. Native modulation maps lie0.019pixel from that neighbor versus1.194pixels from the nominal target. Strong known-source contamination/alias evidence. |
| TIC182517223 | Large primary and secondary eclipse sequence matches a published Gaia SB1 orbit, P2.207646d and K63.98km/s. Stellar mass function≈0.0599solar masses. Known stellar companion, not a new planet or star. |
| TIC353933809 | An isolated S86dip survives a same-detector ensemble control but fails stronger native pixel tests: background treatment changes its shape/sign, and image-scale+gradient fits leave only≈0.5–0.9ppt with≈0.4–0.7ppt uncertainty. Known long-period SB1; no secure inner transit. |
| TIC71274048 | Strong QLP1.359d trough reverses in another aperture; native pixel maps have maximum S/N<2.5 and reference-image flux loss agrees with zero. No secure astrophysical transit source. |
| TIC183373155 | Two deep QLP events suggest a nonunique15.71d alias. S96native differences are sky dominated; S105differences are bipolar. Reference-image/gradient models yield≈0±0.3–0.7ppt, preserving the positive TIC229146531control. No common stellar transit source established. |
| TIC307085108 | Broad1.54872d troughs are approximately4×a0.3865d waveform; the earlier1.16127d peak is≈3×that waveform. Out-of-event low-order harmonic modeling weakens residuals, supporting variability/alias interpretation, without excluding a smaller superposed transit. |

Other leading cases fail clean raw/local recurrence, are known binaries/variables, or are background dominated. TIC375090561's8.136d peak is a **3:2 alias** of a previously public5.424d signal, confirmed against the actual SPOC DV XML. TIC382190230 has previous binary recognition. TIC293268396's initially apparent image offset was **invalid**: low-S/N, boundary-limited Gaussian fits cannot locate a source. TIC238182368's sparse event maps are largely common background; an external event in a gap is not a nondetection.

Read the [first follow-up report](reports/multisector/followup.md), [northern pooled tests](reports/multisector/pooled_north.md), [southern pooled tests](reports/multisector/pooled_south.md) and [second-batch/final pooled report](reports/multisector/batch02_followup.md). The final queue alone added14stars /30observations of conditional raw checks and targeted external-sector/native-pixel investigations. The general follow-up table contains166diagnostic variants for66TICs /153distinct target-sector pairs; the separate final queue contains30observation tests. These overlap the survey and are not independent discoveries.

## Prospective selection and methods

The initial batch selected1,024observations per sector before viewing flux. From each complete official QLP inventory,8,192hash-ordered TICIDs received public metadata retrieval. Fixed salt `multisector-v1:` determines order. Eligibility: finite Tmag≤12.5, Teff3000–6500K, radius0.15–2.5solar. The first1,024eligible IDs per sector form the immutable selection. See [archive/version evidence](reports/multisector/archive_probe.md), [selection](data/multisector/selection.csv) and the original manifest.

After actual recovery/follow-up evidence, a separately fixed brighter/cooler batch used new S96/105inventories,16,384metadata IDs per sector, Tmag≤11.5, Teff3000–6000K and radius0.15–1.2solar. All456/507eligible targets were searched. Full parent inventories contain419,834/399,217entries; those parent counts are **not analyzed stars**. Official QLP106/107inventories returned404. Sector105inventory and actual FITS exist even though the documentation table stopped104. Timing headers/orbits were checked; this is not an independent timing calibration. QLP99–104warning sectors were excluded; correctedv02products are required where applicable for74–79. See [archive probe](reports/multisector/latest_probe.md).

The older Sector1 sources/results remain frozen. Separate wrappers use two gap-separated median widths0.75/1.5days,8peaks,32coarse refinement seeds and more durations. Recurring local diagonal S/N≥5.5 seeds receive follow-up. These looser cuts admit false alarms. Two detrenders share photons. Single-sector nominal period cap20d is additionally limited to halfbaseline; pooling extends to40d where the combined baseline permits. Blind pools include all repeated usable targets, plus residual searches after masking recognized transits where ephemerides permit. Fitted pooled sector consistency is conditional, not independent validation.

**1239 observations** have adequate finite Q0flux but unusable producer errors. The original error-based failures remain saved. A separate versioned branch estimates a constant uncertainty from1.4826MAD(adjacent differences)/√2, restricted to closely spaced cadences with sufficient pairs. Raw FITS/errors/flags remain untouched. Residual/event-scatter floors reduce some formal overconfidence, but this is not a calibrated correlated-noise model. The original survey is not retrospectively called successful on those files.

## Sensitivity, recognition and statistical limitations

[Actual-noise recovery controls](candidates/multisector_controls/README.md) include288distinct injections evaluated by both methods:576scans, not576independent injections. Retention in the SPOC background rises10/144→33/144, but only10→12with a common S/N8cut. Most gain is looser selection. **No300ppm injections were retained** in these tested backgrounds; the faint Sector1background gives0/144for either method. These two backgrounds/finite phases do not establish population completeness. Three genuine shallow planets recover their period/epoch in9/9control observations, with8/9retained. The supplied, previously reported TIC4206066signal is retained in S32/98; fixed known timing accumulates≈362±34ppm across3sectors. That is a conditional recovery control, not a new discovery.

[Thirty-two block-surrogate searches](candidates/multisector_controls/NULLS.md) from8preselected stars give7retained cases concentrated in2stars. Astrophysical fragments, block joins and nonstationarity prevent a formal survey null/FAP. No tiny p-value or planet probability is justified. A BLS maximum, large Δχ², attractive fold, catalogue absence or generic warning flag is insufficient for classification.

Fresh16October NASA confirmed/TOI and ExoFOPCTOI snapshots, RAVEN tables, published candidates, GaiaNSS/variability, VSX, SIMBAD, TESS-EB, literature and aliases support recognition. Leading cases have source-specific scope/results. Multiple finite nonmatches still do not prove novelty. Known hosts at different periods remain eligible. Integer harmonic matching misses3:2aliases; rational period agreement is a review flag rather than automatic exclusion because genuine multiplanet resonances exist.

All uncertainty estimates remain conditional on selected periods/templates/baselines and approximate noise. No transit-derived mass, precise component radius, gravitational binding or confirmed planet is measured. Native pixel localization resolves a field area, not every component of an unresolved system. Empirical gradient subtraction can remove blended genuine signals; its successful positive control is informative, not proof of a complete false-positive model.

## Preservation, verification and continuation

Raw observations, prepared curves, live DBs and rejected diagnostic bytes remain under `/workspace/tess-research-data/multisector/`. External persistence of the cloud overlay is unverified. Public URLs, exact raw hashes, TIC/sector/product identities, FITS timing/processing metadata, native FFI identifiers and query snapshots support retrieval. No datasets were discarded to make room.

The main compressed SQLite checkpoint contains all4096original records. A backup bug involving an uncheckpointed destination WAL was caught by exact restoration, fixed by closing/checkpointing the destination before hashing/compressing, and regression-tested with an active WAL writer. The failed archive/manifest are preserved separately; they are not called valid. The corrected archive restores all4096rows without sidecars and verifies raw hashes. Supplemental archives preserve three additional DBs, every completed pooled JSON and exact follow-up responses. Both metadata query archives have exact offline restore checks. [Diagnostic preservation](data/multisector/diagnostic_preservation.json) maps bulk PNG/NPZ bytes retained outside ordinary Git; important TIC229146531arrays/plots and representative falsification figures remain tracked.

Use the prepared Python environment and one BLAS/OMP thread per process. [PROGRESS.md](PROGRESS.md) provides exact verify/restore/resume commands. Do not rerun the finished surveys as though they added stars, modify frozen source hashes to accept mismatches, or treat gap coverage as a physical null. A known zero-event helper exception has explicit preserved local/coverage fallbacks. An ill-conditioned high-order Fourier baseline and weak Gaussian source offsets were rejected as invalid diagnostics.

The next scientific priorities are classification and independent detrending of TIC229146531, correlated-noise/shallow-signal sensitivity improvements, then a new pre-flux sector batch using the documented timing exclusions. Most attractive broad-search peaks have now received stronger falsification; obtaining a trustworthy300–500ppm search is more useful than scaling the same biased score alone. There is ample public data beyond this finite sample. Parallel worker usage limits were reached; local analysis, extra sector/pixel checks and preservation continued through a meaningful phase conclusion.

OpenRouter research expenditure remains **$0**. No paid external services, external messages, catalogue submissions or observing requests were used.
