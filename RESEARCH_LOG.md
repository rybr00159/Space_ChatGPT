# Research log

## 2026-10-13 — Repository, resources, and source review

Fetched main into the previously empty local working branch. The repository contained only `Reddit Post.zip`; retained it unchanged. Runtime: five CPUs, about 32 GiB available RAM and 30 GiB free disk. `/workspace` is overlay; no separate durable volume identified. Installed Python astronomy packages in the existing virtual environment, with cache outside the repository.

The Reddit post concerns TIC 4206066 and publicly described signals near 3.18 and 11.13 days. These are prior reports and excluded from novelty claims. Its repeatability claims do not eliminate blends or instrumental effects. No paid research calls made.

Initial archive file/NASA TAP checks responded. A raw MAST invoke URL without a valid request returned HTTP 500; this did not test the query service. Proper CAOM filtered POST subsequently succeeded. An initial ExoFOP timeout and unsupported NASA SQL column were diagnosed; supported NASA TOI and ExoFOP CTOI requests succeeded later.

## 2026-10-13 — Pilot selection and real observations

Retrieved and hashed NASA snapshots: 4,811 confirmed transiting planets, 6,445 all-confirmed planets, 8,148 TOIs, and 5,138 CTOIs. Selected eight small known two-planet cool-star hosts: TOI-776, 1064, 1468, 1470, 2095, 4342, 1803, 904. Downloaded their latest three available SPOC sectors; TOI-1803 has two. L 98-59 sectors 2, 5, 8 are positive controls. Current TOI data already include TOI-2095.03 at 57.6869 days.

A product suffix bug (`_s_lc` versus actual `-s_lc`) initially selected zero light curves and was corrected before search. Astropy checksum warnings were diagnosed: SPOC files have CHECKSUM but no DATASUM; the verifier assumes zero data sum. Standard raw header-plus-data ones-complement checks independently validate the original HDUs. No files, expected hashes or verification settings were changed to pass. SciPy rejected FITS big-endian arrays; converting aligned scientific arrays to native floats fixed it.

Independent southern field selection used a 3° cone at RA 100°, Dec −60°, Tmag 9–11.5, Teff 3200–4800 K, radius 0.2–0.8 solar, finite contamination ≤0.05, catalog exclusions and deterministic brightness ordering. Initial six targets have 21–44 available SPOC sectors. Selection preceded flux inspection. Later expansion accepted all ten eligible stars with SPOC coverage (12 catalog-eligible, two without SPOC), four recent sectors each. TIC 150299576 has only earlier sectors, which were used rather than inventing recent coverage.

Final pilot: 19 stars and 66 raw light-curve products. Raw observations/intermediates remain outside ordinary Git; exact observation IDs, product URIs, processing metadata, hashes and retrieval details are tracked. No important scientific data deleted.

## 2026-10-13 — Controls, injection checks, and first screening

Free BLS recovers L 98-59 c at 3.6909, 3.6905, 3.6917 days in sectors 2/5/8, diagonal S/N 28.0/30.1/23.2. Half/double-period aliases demonstrate that rankings alone are insufficient. Direct ephemeris/local-baseline checks recover all three known transiting planets; sparse sampling reduces some measurements substantially.

Six nonzero pre-detrending blind injection scans (500/1000/2000 ppm × two median widths) recover a 4.731 day, 1.8 hour injected box as rank one, retaining about 97–98% differential depth. Two unchanged-data controls make eight total runs. A favorable epoch was selected by coverage, not noise depth; this is not completeness. An earlier shorthand calling these eight injections was corrected. Five tests cover checksums/corruption, timing, aligned indexing, masks, gaps and local fits.

Initial residual maxima are mostly diagonal S/N 5–7. Stronger TOI-1064 and TOI-4342 maxima select one-event/gap-boundary excursions. Searches continued after removing those edges rather than treating the first highest rank as a planet. Known epochs and BTJD/TDB headers were checked; gross time conversion error does not explain the leading artifacts.

## 2026-10-13 — Catalog/literature and contamination findings

Public host-label arXiv query yielded 37 papers. Literature finds prior candidates missing from confirmed/TOI catalog entries: TOI-4342 at about 47.5 days (RV), L 98-59 at 1.7361 days (RV). TOI-776/1470 rotation and TOI-1803 timing variations provide additional cautions. TTV sensitivity used widened linear masks; measured individual timing ephemerides were not reconstructed.

All ten field targets lack exact TIC/30″ matches in the checked confirmed/TOI/CTOI snapshots. This does not establish novelty. SIMBAD/Gaia identify activity/companions and solar-like variability flags for some targets. TESS EB web search hits were initially misread as membership by a reviewer; detail pages explicitly say `In catalog? False`. That mistaken interpretation was corrected and evidence preserved; no known-EB claim is made.

Focused TIC 150299576 catalog, live ExoFOP, nine-alias arXiv and eight SIMBAD-reference abstract checks find no prior planet report within that finite scope. A neighbor at 7.83″ with flux ratio about 0.0058 could dilute a 10–20% eclipse into the weak observed dips if its flux enters the aperture. No claim that the neighbor actually eclipses or that novelty is established.

## 2026-10-13 — Noise comparisons and lead falsification

Completed 96 half-day block-permutation full-grid scans: 24 per width for each of two stars. Paired-width illustrative tail fractions are 0.24 for TIC 278727563 sector 94 and 0.76 for L 98-59 sector 2. These are limited exchangeable-block comparisons, not global false-alarm probabilities. Search maxima near S/N 5 are ordinary under these tests.

TOI-1064/4342 original maxima fail gap, reduction and sampled-prediction checks. Revised TOI-1064 8.87/14.85 day leads remain reduction dependent; the 14.85 day last event has negative native SAP depths under multiple local baselines. TOI-4342's 6.67/13.34 day ambiguity lacks robust raw/recurrence support; baseline drift prevents treating a rolling-median SAP sign alone as definitive exclusion.

The final TOI-4342 19.45 day joint maximum has corrected-flux local formal S/N near 9 but native SAP linear/quadratic S/N −1.37/0.62. Its two supported cycles have only 2 versus 31 ten-minute samples and the same parity: first-event coverage is partial and odd/even agreement unavailable. No credible repeated transit interpretation is established.

TIC 348839788 has a broad ~0.27228 day waveform and BLS multiples. TIC 150359500 has a broad ~1.033 day, ~4% waveform agreeing in raw flux. Thirty-six frozen training ephemerides lack robust recurrence across both already-screened heldout sectors; these are not pristine blind validation. A ~0.68977 day residual fails wider-detrending/local-baseline tests, correlated beta ~3.5. Preservation injections retain roughly 80–86% of 1000 ppm boxes.

TIC 150299576's ~11.3 day two-event lead weakens in a four-sector coherent search. Thirty of 64 regional moving-block null maxima exceed its joint statistic (plus-one 0.477, not global FAP). Native SAP first-event local estimates are much smaller than corrected flux. TIC 150360085's 13.123 day score 6.86 becomes local-quadratic 2.04, with mixed event signs. TOI-904's 13.46 day score 7.93 becomes local-quadratic 1.18 and native SAP 0.51; its broad duration is a consistency caution, not proof against a planet.

TOI-1803's broad 4.59 day score 8.31 falls to local-quadratic 2.66 and SAP 0.88. A 945-combination training-only narrower-box safeguard peaks at event-scatter S/N 2.36. Its 1000 ppm preservation controls retain approximately 976–1035 ppm. Attractive broad joint scores do not become convincing local events.

## 2026-10-13 — Negative pilot checkpoint and reproducibility review

Completed 135 per-sector variants (132 residual + three unmasked controls), 16 additional untrimmed joint variants from the incremental archive snapshot, and 58 gap-edge-excluded variants covering the final 19-star sample. Retained 775 ranked screening peaks. No credible new candidate promoted. Lower ranks are not individual planetary exclusions; physical fits, pixel localization and global population false-positive/completeness analyses remain unrun.

Independent review confirmed all 66 raw hashes and 198 raw HDU checksum records, unique product URIs, TIMESYS=TDB and BJDREF=2457000. Remaining stale numeric-zero odd/even fields from earlier in-flight processes were marked null/unavailable when one parity lacked coverage, then tables rebuilt. This changes no search maximum. Field reproduction command now explicitly preserves the final limit of 12, yielding ten covered stars.

Preserved exact historical search source from Git commit `7ecc52c` as `search_null_v1.py`. Hash-verified `null_calibration.py --frozen-search` resumes zero pending permutations; all 96 measurements/paired summaries unchanged. Production default rejects incompatible checkpoints. No verification bypass.

Executed the saved installation script successfully; dependency checks, astronomy imports and five utility tests pass. Saved install_script/start_skill in environment draft. Saving is not publication; durable external storage and fresh-task restoration remain unverified. Cache including intermediates is about 154 MB, retained. Research-phase paid calls/spending zero.

`RESEARCH_REPORT.md`, `PROGRESS.md` and the candidate index document the limited negative result and exact next quiet-star/shorter-transit experiment. No discovery claim and no assertion that these stars lack further planets.

Final ledger review derives L 98-59 stellar properties directly from the same current confirmed catalog (Tmag 9.41242, radius 0.3155 solar, Teff 3415 K), replacing initial summary placeholders; photometric searches and measured depths are unchanged. Generated summary CSVs now use LF and include explicit odd/even availability. Original downloaded catalog bytes/hashes remain unchanged.

## 2026-10-13 — Broad survey replaces the small-pilot stopping point

The user requested much larger coverage, citing a 158,263-star comparison. Selected 160,000 unique targets by a fixed random permutation (seed 20261013) from the authoritative 469,735-target QLP sector-1 inventory, before flux inspection. Inventory counts are separate from processed photometry. Raw 30-minute files are approximately 80,640 bytes (12.9 GB for the sample), compared with newer SPOC sector-83 files of about 1.1 MB (approximately 177 GB for 159,997 stars). This supports the compact broad track within this environment. Every raw download is retained. Official timing warnings exclude QLP sectors 99–104; sector-1 timing headers are verified.

The fast method benchmarks actual QLP data at approximately 0.04 seconds per star, using cadence-aware segments, quality-zero rows, a one-day median, conservative diagonal errors, and an adaptive logarithmic period grid. Synthetic controls recovered 36/36 idealized randomized 30-minute injections and 3/3 known L 98-59 c controls; 64 white-noise nulls and eleven utility tests were also run. These are not real-population completeness estimates or global false-alarm probabilities. Prior checks include 1,160 DIAmante rows, 185 Olmschenk and 288 Yu published rows, and 4,578 genuine published TESS EB members, alongside NASA/TOI/CTOI catalogs. RV-only signals and nearby contaminants are distinguished.

The first 256 light curves were searched in nine seconds, then reanalyzed from cache after source-hash tracking was introduced. SQLite unique TIC keys avoid duplicate counts. The full run uses 32 pooled HTTP workers, four CPU workers, and a bounded 192-item queue. Initial throughput was approximately 66 stars per second; inadequate light curves receive a separate status. The recent official sector-107 TCE table (846 signals on 623 stars) supplies a parallel interpretation track: these are publicly reported detections, not original discoveries by this project.

## 2026-10-13 — Focused falsification during the wide survey

TIC 370247818 initially appeared strong. Thirty-nine event timings across eight QLP sectors (1, 5, 12, 13, 27, 28, 32, 39) yield period 3.522194945921 ± 0.000005834 days and epoch BTJD 1326.916296 ± 0.001278. An exposure-integrated approximate fit gives depth about 8,485 ppm and duration 2.14 hours. Raw flux, alternate apertures, and sector-depth consistency support a recurring astrophysical signal. Three independently retrieved TessCut cubes in different camera/CCD orientations localize the dip to the target neighborhood and disfavor the bright red giant 86.6 arcseconds away. Unresolved blends and PRF systematics remain possible. Raw cubes, original FFI identifiers, compact arrays, and 500-bootstrap localization diagnostics are preserved.

An all-phase search finds a secondary near phase 0.47546, approximately 2.1 hours before half phase. Its raw SAP depth is 1,327 ± 232 ppm, with alternate-aperture and multi-era support. A half-phase-only check missed this feature. The circular primary fit favors a grazing geometry: radius ratio 0.278 implies approximately 0.468 solar radii conditionally on catalog stellar properties. The broad profile gives no precise upper radius constraint. Illustrative primary-only eccentric small-companion solutions predict no secondary overlap and fail the combined eclipse evidence. The conditional secondary light budget further favors a luminous stellar companion. A grazing eclipsing binary is favored; no mass measurement or planetary classification is claimed.

TIC 141980042 recurs across 23 sectors and 192 event timings (period 2.32207006 days), but a pooled raw secondary of approximately 894 ± 86 ppm and orbital variation favor an EB or blend. Gaia has a position-only solution and 71% multipeak windows, so its old QLP radius is not secure. Neither lead is promoted as a planet.

Real-data validation used 128 catalog-unflagged stars, 512 block/sign surrogates, and 128 paired-depth injections on 64 stars. Guarded surrogate leads occurred in 9/512 initial or 13/512 stricter local tests; these dependent, approximate surrogates do not estimate global false-alarm probability. Exact period/epoch recovery was 26/64 at 1,000 ppm and 49/64 at 3,000 ppm, substantially poorer for noisy stars and longer periods. Negative screens do not exclude weak or long-period planets. Source hashes and raw integrity were verified; paid research calls remain zero.

## 2026-10-14 — Resume an interrupted survey and strengthen raw-flux checks

The session interruption stopped the survey after 118,000 completed statuses: 117,638 stars actually period-searched and 362 downloaded with insufficient usable data. No archive errors were recorded. Resumed the exact fixed selection from SQLite, leaving 42,000 targets to process and retaining the existing cache.

The initial ranked batch independently vetted 187 unique stars, with zero processing errors. Thirteen had no first-sector automated warning; this was not planetary validation. TIC 308600649 shows strong raw SAP ellipsoidal modulation suppressed by the corrected-flux detrending, motivating an explicit raw orbital-variation guard. Prior multi-sector checks already reject TIC 141980042. Requested additional QLP sectors for TIC 260352037, 441386747, 406280158, and 271697430, with respectively 24, 1, 2, and 24 additional observations retained. Their exact selected peaks, rather than automatically the highest BLS peak, must seed follow-up.

The sector-107 alternative track retained all 846 original TCE rows and selected 71 rows on 63 stars for potential follow-up. Actual reanalysis of nine SPOC curves on three stars rejects two variability aliases; TIC 178155030 is a published complex periodic variable in Bouma et al. (AJ 167, 38, DOI 10.3847/1538-3881/ad0c4c). Its catalog absence as a planet does not indicate novelty. Further clean TCEs are being examined independently.
