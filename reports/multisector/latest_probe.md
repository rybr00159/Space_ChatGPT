# Latest safe-sector availability and pre-flux target sample

Checked October 16, 2026 (America/Halifax; individual retrieval manifests use UTC). This reconnaissance produced target-selection inputs and four real-product interface checks. Its inventory and metadata counts are **not completed period searches**. Screening of the separately frozen `b02_096_105` batch is recorded elsewhere.

The live [official QLP page](https://archive.stsci.edu/hlsp/qlp) still labels Sectors 99–104 with asterisks, states their corrected timestamps have not yet been redelivered, and excludes them from this search. The page's availability table ends at Sector 104, although its year heading says 97–107. Actual archive access is more informative here: the complete **Sector 105** target list and two real light curves return HTTP 200. Sector 106 and 107 target-list URLs return HTTP 404, with the exact responses retained; absence of those inventories does not establish the absence of every possible data product. **Sector 96** is therefore the available second sector. No corrected timing release is assumed for 99–104, and no ephemeris recomputation independently certifies the later Sector 105 timestamps.

| Sector | Complete official inventory rows | Fixed metadata IDs | Metadata returned | Original broad eligibility | New bright/cool/small eligibility |
| --- | ---: | ---: | ---: | ---: | ---: |
| 96 | 419,834 | 16,384 | 16,384 | 3,765 | 456 |
| 105 | 399,217 | 16,384 | 16,384 | 4,022 | 507 |

Every selected metadata ID returned exactly once, with zero missing responses. Original eligibility is finite Tmag ≤12.5, Teff 3,000–6,500 K and radius 0.15–2.5 solar radii. The new prospective filter tightens those bounds to **Tmag ≤11.5, Teff 3,000–6,000 K and radius 0.15–1.2 solar radii**, improving the observational prior for detecting small planets around brighter, smaller hosts. These are catalog estimates, not confirmed stellar parameters. The tighter sets contain 963 sector observations on **850 distinct TICs**, including 113 targets in both sectors. They are selected before any flux inspection or signal ranking.

The frozen inventory implementation, `src/multisector_inventory.py`, is reused through imports without scientific changes. It selects the first 16,384 TIC IDs in ascending `SHA256('multisector-v1:' + decimal_TIC_ID)` order from each **complete** official target list. Public TIC requests are exact-ID batches of 512 and preserve their raw JSON bytes, request definitions, retrieval times and SHA256 hashes. The normalized files retain the original broader eligibility for audit; the root batch planner separately freezes the tighter selection.

## Actual raw observations

The first two deterministic hash-ordered TICs in each complete inventory were downloaded. This interface probe does not select the scientific sample using signal properties.

| Sector | TIC | Native rows | Quality-zero finite positive corrected flux | Usable positive producer errors | Orbits | Raw bytes |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| 96 | 234490958 | 10,819 | 9,727 | 0 | 199, 200 | 659,520 |
| 96 | 301868290 | 10,819 | 7,148 | 0 | 199, 200 | 659,520 |
| 105 | 234490958 | 11,723 | 8,688 | 0 | 221, 222 | 714,240 |
| 105 | 304986222 | 11,723 | 8,407 | 0 | 221, 222 | 714,240 |

All four parse as FITS, pass structural verification and match the requested TIC and sector. Their schema is `TIME CADENCENO SAP_FLUX DET_FLUX DET_FLUX_ERR QUALITY ORBITID SAP_X SAP_Y SAP_BKG SAP_BKG_ERR DET_FLUX_SML DET_FLUX_LAG SYS_RM_FLUX`. They explicitly declare `TIMESYS=TDB`, `TIMEREF=SOLARSYSTEM`, BJD reference 2,457,000, and `TIMEDEL=200/86400` days. Measured median cadence is 199.996–200.000 seconds. Approximate timestamp ranges are BTJD 3907.77–3933.45 in Sector 96 and 4205.41–4233.18 in Sector 105. Calibration is TICA v1.4.0, with archive product filename version v01 and TIC version 8.2.

**Every `DET_FLUX_ERR` value in all four probes is NaN.** A producer-error-weighted search therefore has no usable measurements in these particular files. The separately labeled, already implemented empirical adjacent-difference uncertainty branch can analyze retained corrected flux without altering producer flux, errors or QUALITY; it remains an assumption-dependent sensitivity analysis, not repaired producer data or calibrated significance. The four probes do not establish that all files in either sector share this issue.

No file contains native FITS CHECKSUM or DATASUM cards. Producer checksum verification is therefore unavailable; local SHA256 establishes the exact retained bytes. The explicit time reference does not independently prove the physical accuracy of barycentric corrections.

## Preservation and reproduction

Full successful inventories remain at `/workspace/tess-research-data/multisector/inventories/`, raw observations under `/workspace/tess-research-data/multisector/raw/s0096/` and `s0105/`, and all public TIC responses under the existing `tic_metadata/` cache. Nothing was deleted. The successful full inventory manifests record the public static URLs, exact sizes and SHA256 hashes; changed future responses must be preserved separately rather than silently substituted.

Machine-readable probe details and original failed HTTP responses are recorded in `data/multisector/latest_probe/` and its exact query archive. The compressed archive retains all 64 raw TIC response bodies with adjacent manifests, the official status page and failed later-sector inventory checks. It does not duplicate the much larger successful full inventory CSVs or raw FITS inside ordinary Git. The normalized target rows and full metadata request manifests are versioned under `data/multisector/inventory_s0096_*` and `inventory_s0105_*`.

```bash
# Cached probes are verified against their original bytes.
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/latest_sector_probe.py probe

# Reuse frozen metadata, or reconstruct it only if the named outputs are absent.
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/latest_sector_probe.py metadata --sectors 96 105 --size 16384

# Verify or restore exact archived public query snapshots.
/workspace/openrouter-check-venv/bin/python scripts/latest_sector_probe.py verify
/workspace/openrouter-check-venv/bin/python scripts/latest_sector_probe.py restore-queries

# Recover missing full successful inventories, enforcing their original SHA256.
/workspace/openrouter-check-venv/bin/python scripts/latest_sector_probe.py download-inventories
```

Existing named metadata outputs are reused; the helper refuses a different selection size or method. Rebuilding missing outputs should be done only when the batch writer is stopped. Do not change source hashes or overwrite a frozen selection to extend a run. `probe` repeats the structural checks using the same cached bytes and preserves the original descriptive checkpoint timestamp; it does not refresh the retained official status page. A new date's release-status check needs a new named cache if the original response is to remain frozen.

The final archive contains **136 exact response/manifest files** (10,130,616 compressed bytes). Verification passed all 136 retained cache members and 10 source/artifact hashes. An offline restoration to a fresh temporary directory restored all 136 files and passed the same exact-byte checks. Cached metadata reuse and full-inventory SHA checks also passed.
