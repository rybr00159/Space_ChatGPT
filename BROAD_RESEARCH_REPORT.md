# Broad TESS investigation — October 2026

This phase replaces the 19-star pilot with a fixed 160,000-star photometric sample and progressively more expensive validation. No confirmed planet or discovery is claimed. The completed broad export is audited in `results/broad/audit.json`; detailed current-property follow-up continues.

## Scientific strategy

The pilot compared known multiplanet hosts with a small field-dwarf sample. Known hosts offer established properties and a useful residual-search hypothesis, but have extensive prior analysis and known-event masks that reduce sensitivity. For broad coverage, compact public QLP sector-1 observations permit a six-figure sample within this machine's storage. A seeded random permutation selects 160,000 unique TICs from the official 469,735-target inventory before looking at flux. This is a random sample of the archive list, not a random sample of all Galactic stars. Sector 1 has already been extensively searched; its size and accessible data do not establish novelty.

QLP v01 30-minute products require approximately 12.9 GB for this sample. Representative newer 200-second TESS-SPOC products would require approximately 177 GB for a similar number of stars, exceeding this environment's disk. Higher cadence, more detrenders, and additional sectors are reserved for promising leads. Timing-flagged QLP sectors 99–104 are excluded. A parallel recent-sector track reanalyzes selected public sector-107 threshold crossing events (TCEs); these archive-reported detections are clearly distinguished from original detections by this project.

## Broad detection method

`src/broad_survey.py` combines pooled HTTP requests, a bounded queue, four CPU workers, and per-star SQLite commits. A frozen implementation hash prevents silently mixing changed search methods. Public product URLs, original byte hashes, sectors, processing/header metadata, and retrieval times are retained. All raw downloads remain outside ordinary Git.

`src/fast_screen.py` uses quality-zero finite positive flux/error rows, cadence-aware gap segmentation, a one-day robust median, and diagonal errors at least as large as measured robust scatter. Positive flare excursions are excluded; downward events remain available for inspection. The adaptive logarithmic BLS grid searches 0.3–13 days, subject to half the usable baseline, with trial durations resolved at the actual cadence. At 30-minute cadence the useful nominal trials are 1.2, 2.4, and 3.84 hours. Two separated peaks and local isolated-event summaries are retained. The initial guards require sampled positive events, limit duty fraction, and reject a single event dominating the positive squared score.

Known transit ephemerides are masked in a separate residual scan for eligible known hosts. Confirmed planets, TOIs, CTOIs, RV-only signals, published DIAmante/Olmschenk/Yu candidates, and genuine published TESS EB membership are checked separately. Catalog matching uses TIC/position and period/harmonic tolerances, with nearby objects marked as potential contamination. A missing catalog match is an unresolved prior-status result, not proof of novelty.

`src/broad_prioritize.py` applies stricter local recurrence cuts, preserves exact scan/peak evidence identifiers, and creates multiple warning/control groups. Current TIC stellar properties supersede old QLP headers for follow-up ranking. Several old approximately 0.7–0.8-solar-radius entries are currently classified as 10–12-solar-radius giants; ignoring that change produces misleading companion radii. Missing or negative sentinel properties never become physical estimates. Size and duration proxies are prioritization aids with dilution/eccentricity limitations.

`src/broad_vet.py` refines saved peaks with independent baselines, three detrending scales, four apertures including native SAP, event-level depths, odd/even tests, integrated approximate transit shapes, centroids, and secondary searches across all orbital phases. Its newer raw-SAP harmonic guard addresses stellar modulation removed from corrected flux. Every ranked result remains unvalidated until independent data and astrophysical alternatives have been examined.

## Validation and sensitivity limits

All eleven repository utility tests pass. Fast-search controls recover 36/36 idealized randomized cadence-integrated injections and all three known L 98-59 c controls. These checks establish specific pipeline behavior, not population completeness.

More realistic tests use 128 stars selected before outcomes, 512 structured block/sign surrogates, and 128 paired 1,000/3,000-ppm injections on 64 actual stars. Exact period-and-epoch recovery is 26/64 at 1,000 ppm and 49/64 at 3,000 ppm; adding local recurrence cuts reduces those counts to 11/64 and 27/64. Recovery is substantially worse for noisy stars and periods exceeding six days. See [the noise assessment](reports/broad_noise.md) for parameters, stratification, and raw hashes.

Nine initial-guard or thirteen stricter-local surrogate leads occur among 512 surrogate trials. The transformations are dependent and imperfect representations of real noise, so these are diagnostic failure rates, not global false-alarm probabilities. Even 512 ideal null trials cannot resolve a six-figure survey's extreme tail. Searches over many stars, periods, durations, detrenders, and chosen follow-ups incur selection effects. Formal BLS/local scores and conditional fit errors must not be interpreted as globally calibrated discovery significance.

The broad survey is incomplete for weak, short, long-period, single-event, highly variable, and lower-ranked signals. Keeping two peaks can bury planets beneath stellar variation. Thirty-minute integration attenuates short events. No negative screen excludes additional planets in a target.

## Strong-signal falsification

**TIC 370247818:** a genuine recurring eclipse was initially an attractive lead. Eight sectors and 39 event timings yield period 3.522194945921 ± 0.000005834 days. Three independent TessCut reductions localize the dip to the target neighborhood and disfavor the bright resolved neighbor 86.6 arcseconds away. An all-phase search then finds a secondary at phase 0.47546, with raw SAP depth 1,327 ± 232 ppm and alternate-aperture support. The half-phase-only test missed it. An exposure-integrated primary fit favors grazing geometry; the radius profile is broad and cannot precisely classify the companion. Conditional secondary light-budget comparisons favor a luminous stellar companion, and illustrative eccentric small-companion primary fits fail secondary-overlap geometry. A grazing EB is favored; no mass confirmation is available. Evidence includes [multi-sector measurements](candidates/broad_s1/tic370247818_multisector/README.md), [independent pixels](candidates/broad_s1/tic370247818_pixels/README.md), and [physical assessment](candidates/broad_s1/tic370247818_physical/README.md).

**TIC 141980042:** 23 sectors and 192 event timings support period 2.32207006 days, but a raw secondary near 894 ± 86 ppm and orbital variation favor an EB/blend. Current Gaia has a position-only solution with strong multipeak diagnostics; the old stellar radius is unreliable. [The dossier](candidates/broad_s1/tic141980042_multisector/README.md) preserves measurements rather than promoting this as a planet.

**TIC 308600649:** a first-sector shape survived the earlier automated checks, but raw SAP shows strong half-period modulation. Current Gaia RV variability and revised stellar/crowding properties supply further binary/blend warnings. This motivated the raw orbital-variation guard instead of treating a corrected-flux U-shaped dip as sufficient evidence.

**TICs 441386747, 406280158, 271697430, and 260352037:** current TIC radii of approximately 10–12 solar radii invalidate their initial small-companion interpretation. Independent-sector observations and current catalog evidence are retained. Giant-star variability and contamination must be addressed before planetary interpretation.

## Recent public TCE track

The official sector-107 table contains 846 TCE rows on 623 stars. Metadata triage retains 71 rows on 63 stars after known-period, sentinel, and within-star secondary cuts. Those counts describe public metadata, not 623 newly searched light curves. Actual native SPOC reanalysis initially covers nine curves on three stars: two signals are variability aliases, and TIC 178155030 is a published complex periodic variable (Bouma et al., AJ 167, 38; DOI 10.3847/1538-3881/ad0c4c). Its narrow recurring dips are not a new planetary discovery.

Eight additional curves on TICs 49617264 and 89187590 were retrieved and analyzed. TIC 89187590's nominal 8.726-day dips fail local-curvature/phase checks amid a roughly 1.1-day stellar waveform. TIC 49617264 has recurring dips across four sectors and plausible current dwarf-star properties. Wider prior investigation finds publicly uploaded 2024 ground follow-up identifying this exact star as an on-target candidate, TIC 49617264.11 / TOI-2417.11. The brighter star 43.18 arcseconds away was the original false-positive TOI-2417.01; assigning the dips to that neighbor merely because its catalog period matches would contradict the ground evidence. Independent target pixels support the faint target as the eclipse source. This is a recovered, previously reported candidate, with unresolved companion classification, and fails the project's novelty objective. Exact source notes and pixel diagnostics are retained in `data/broad_catalogs/` and `candidates/tce107/`.

## Preservation and reproducibility

The frozen sample, compressed target/signal tables, observation manifest, catalog source snapshots, candidate arrays/plots, and code are stored in Git. Raw observations, pixels, and SQLite remain under `/workspace/tess-research-data`. No important raw evidence was deleted. Raw QLP FITS checksum cards are absent; this is explicitly reported, while SHA-256 identifies the exact retained bytes. SPOC/TessCut products with checksum cards are checked against their original raw HDUs.

`src/restore_broad.py` verifies or restores recorded products with explicit download limits, immutable expected hashes, target/timing checks, and storage reserve. Verification of twenty cached products passed with no network requests. It does not recreate SQLite; the frozen survey can reconstruct searches from retained or restored raw files. `src/summarize_broad.py` audits unique TICs, completed statuses, sample membership, and manifest alignment before creating the survey overview. See [preservation instructions](reports/broad_preservation.md).

This workspace is an overlay filesystem. An independently durable volume and restoration in a fresh task have not been verified. Environment publication can snapshot prepared files; live processes must never be assumed to persist. Updated installation/start instructions are saved in the environment draft. Research-phase paid calls remain zero; the user-authorized OpenRouter cap is $10, with an earlier untracked connectivity-test cost requiring a conservative allowance.

## Final sample audit and next investigation

The stable final export passed independent consistency checks for unique TICs, fixed-sample membership, status counts, signal membership and observation-manifest alignment. The overview was visually inspected.

| Quantity | Completed broad phase |
|---|---:|
| Archive inventory, not analyzed count | 469,735 |
| Predetermined unique targets | 160,000 |
| Actual completed period searches | **159,505** |
| Downloaded but insufficient usable data | 495 |
| Processing/archive errors | 0 |
| Retained raw QLP products | 160,000 |
| Exact raw QLP bytes | 12,902,400,000 |
| Stored screening peaks | 319,496 |
| Stars with stored known-masked residual peaks | 243 |
| Original unique stars receiving detailed v3 checks | 187 |
| Focused giant-star follow-up light curves | 57 |
| Separate actual recent-TCE light curves / stars | 17 / 5 |
| Credible new candidates established to this point | **0** |

These broad counts do not add the historical pilot and do not treat an inventory row, peak, metadata TCE, or independent reanalysis as another unique searched star. Compressed [target](results/broad/targets.csv.gz), [signal](results/broad/signals.csv.gz), [raw observation](results/broad/observations.jsonl.gz), and [focused assessment](results/broad/focused_assessments.csv) ledgers preserve searchable evidence.

The full-sample current-property ranking and additional detailed queue are still being completed. Read `PROGRESS.md` for exact active commands and resumption state. The next scientific phase should inspect any remaining recurring dwarf-star signal with independent sectors and pixels. If none survives, improve sensitivity with a preselected quiet-star multi-sector sample, several detrending scales, and additional independent peaks, measuring recovery on actual noise. Repeating this coarse single-sector screen alone would not resolve its documented sensitivity limits.
