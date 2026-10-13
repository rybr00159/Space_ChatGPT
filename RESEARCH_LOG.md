# Research log

## 2026-10-13 — Initial inspection

Repository main had only `Reddit Post.zip`; fetched it and switched the empty working branch to main without altering the archive. Five CPUs, ~32 GiB available RAM, ~30 GiB disk. No dedicated durable storage identified: `/workspace` is overlay. Public archive and NASA TAP connectivity initially work; raw MAST invoke URL without a request returned HTTP 500 (not a valid service test). Proper query test underway. ExoFOP request timed out; NASA TOI TAP is an alternative. Python astronomy installation underway.

The Reddit post concerns TIC 4206066 and publicly described ~3.18 d and ~11.13 d signals. These are prior reports and excluded from our novelty claims. No candidate discovery claimed. No research-phase OpenRouter requests made.

## 2026-10-13 — Real data pilot and integrity troubleshooting

Retrieved dated NASA confirmed-transiting, all-confirmed, and TOI snapshots. Eight chosen known two-planet cool-star systems: TOI-776, 1064, 1468, 1470, 2095, 4342, 1803, 904; most recent three available SPOC sectors per host (TOI-1803 has two). L98-59 sectors 2, 5, 8 are positive controls. Current TOI snapshot adds TOI-2095.03 at 57.6869 d; not novel. Raw MAST CAOM filtered query works; first suffix match bug (`_s_lc` vs actual `-s_lc`) corrected before search.

Astropy 8 FITS checksum warning diagnosed: without DATASUM it assumes a zero data sum. SPOC products include CHECKSUM but lack DATASUM. Implemented standard raw 32-bit ones-complement sum over header plus padded data; full HDUs verify against 0xffffffff. No header, expected checksum, signature, or TLS verification modified. Early scipy error due to FITS big-endian arrays corrected by converting scientific arrays to native float.

Known-transit BLS control: L98-59 c (~3.6906764 d) recovered at 3.6909, 3.6905, 3.6917 d in sectors 2/5/8 with diagonal-error S/N 28.0/30.1/23.2. Half/double-period aliases appear, illustrating why BLS ranks alone are insufficient. Residual search underway at rolling-median widths 1 and 2 days; separate direct all-known recovery and injection tests commissioned. New candidates: none established.
