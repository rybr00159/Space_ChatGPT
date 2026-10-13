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
