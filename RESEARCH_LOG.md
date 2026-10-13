# Research log

## 2026-10-13 — Initial inspection

Repository main had only `Reddit Post.zip`; fetched it and switched the empty working branch to main without altering the archive. Five CPUs, ~32 GiB available RAM, ~30 GiB disk. No dedicated durable storage identified: `/workspace` is overlay. Public archive and NASA TAP connectivity initially work; raw MAST invoke URL without a request returned HTTP 500 (not a valid service test). Proper query test underway. ExoFOP request timed out; NASA TOI TAP is an alternative. Python astronomy installation underway.

The Reddit post concerns TIC 4206066 and publicly described ~3.18 d and ~11.13 d signals. These are prior reports and excluded from our novelty claims. No candidate discovery claimed. No research-phase OpenRouter requests made.

## 2026-10-13 — Real data pilot and integrity troubleshooting

Retrieved dated NASA confirmed-transiting, all-confirmed, and TOI snapshots. Eight chosen known two-planet cool-star systems: TOI-776, 1064, 1468, 1470, 2095, 4342, 1803, 904; most recent three available SPOC sectors per host (TOI-1803 has two). L98-59 sectors 2, 5, 8 are positive controls. Current TOI snapshot adds TOI-2095.03 at 57.6869 d; not novel. Raw MAST CAOM filtered query works; first suffix match bug (`_s_lc` vs actual `-s_lc`) corrected before search.

Astropy 8 FITS checksum warning diagnosed: without DATASUM it assumes a zero data sum. SPOC products include CHECKSUM but lack DATASUM. Implemented standard raw 32-bit ones-complement sum over header plus padded data; full HDUs verify against 0xffffffff. No header, expected checksum, signature, or TLS verification modified. Early scipy error due to FITS big-endian arrays corrected by converting scientific arrays to native float.

Known-transit BLS control: L98-59 c (~3.6906764 d) recovered at 3.6909, 3.6905, 3.6917 d in sectors 2/5/8 with diagonal-error S/N 28.0/30.1/23.2. Half/double-period aliases appear, illustrating why BLS ranks alone are insufficient. Residual search underway at rolling-median widths 1 and 2 days; separate direct all-known recovery and injection tests commissioned. New candidates: none established.

## 2026-10-13 — First screening checkpoint

Known controls recover all three L98-59 transiting planets at pinned catalog ephemerides (independent local-box validation underway). First residual peaks generally diagonal S/N ~5–7. Stronger TOI-1064 sector94 peak (~12.03 d,S/N20.9) is supported by only one observed event; not accepted as periodic. TOI-4342 sector102 ~9.38d peak has odd/even difference6.3sigma. Both delegated for falsification, gap/known-event/pixel diagnostic planning. The independent literature search also finds prior ~47.5d RV candidate onTOI-4342 outside current confirmed/TOI catalog entries.

Matched exploratory field: 6 deterministically selected TIC dwarfs from3° cone centered RA100°,Dec−60°, Tmag9–11.5, Teff3200–4800K, Rstar0.2–0.8Rsun, catalog contamination≤0.05, excluding confirmed/TOI hosts. TIC348839788,150359500,278727563,150360085,348839628,260476305. These have21–44 SPOC sectors; selected recent four for pilot. This is an exploratory comparison, not a statistically matched completeness study. Observation cache expanding.

## 2026-10-13 — Validation, controls, and falsification

Eight pre-detrending injected-transit blind scans atP4.731d/1.8h recover500/1000/2000ppm rank1 under both median widths; ~97–98% depth retained. Favorable epoch selected using mask coverage, so this is not a completeness study. Five integrity/time/indexing/local-baseline tests pass. 96 half-day block-permutation full-grid searches (24 per width per two stars) demonstrate ordinary trial maxima nearS/N5; paired-width illustrative tail fractions0.24 forTIC278727563s94 and0.76 forL98-59s2 are not global false-alarm probabilities.

TOI1064/4342 leading peaks are segment boundaries, reduction-dependent and fail fixed-epoch heldout checks. TOI4342 trimmed6.6712/13.3424d two-event alias also lacks SAP/heldout support; opposite rolling-median SAP signs alone do not prove astrophysical exclusion because baseline drift is present. Pixel re-reduction remains a possible future revisit, not a current credible candidate.

Field sample expanded by a deterministic rule to all10 eligible stars with SPOC coverage in original3°cone (12catalog eligible,2noSPOC). Four recent sectors each;66total rawLC products for19stars includingcontrol. Strong field waveforms: TIC348839788 fundamental~0.2722d with3/4/5×BLSaliases; TIC150359500 ~1.033d broad~4%modulation. SIMBAD activity/companions and Gaia solar-like variability checks add cautions. TESS EB web search hits were initially misread as catalog membership by a reviewer; detail pages show`In catalog?False`, so no known-EB claim is made. Corrected evidence retained.

Marked odd/even statistics unavailable when one parity has no supported event; formernumericzero from infinite-error propagation cannot be interpreted as odd/even agreement. Added explicit TDB/BTJD reference guard. Weak centralTIC150299576s13~11.3d two-event lead under further vetting.
