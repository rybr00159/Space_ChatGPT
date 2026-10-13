# Broad-survey strategy and archive benchmarks

Research-design probe conducted 2026-10-13. The original 19-star pilot is too small to justify stopping an autonomous discovery search. A broad time-series survey is feasible, but catalog indexing and actual light-curve searches must be counted separately.

## Decision: compact QLP primary-mission light curves first

The best storage-efficient route is the public MIT Quick-Look Pipeline (QLP) archive for a primary-mission sector. Sector 1 contains **469,735 target identifiers** and representative light curves are only **80,640 bytes each**, at 30-minute cadence. A predetermined 158,263-target subset would require approximately **12.76 GB** for raw FITS files if this size is representative. The count is a resource estimate, not completed analysis. Selecting deterministically across the whole target list avoids the pilot's narrow cone and first-few-TIC spatial ordering.

TESS-SPOC full-frame-image sector 83 has **159,997 catalog rows**, but its representative 200-second light curve is 1,108,800 bytes. A complete sector would require approximately **177.4 GB**, exceeding this machine's approximately 29 GB free disk. TESS-SPOC sector 14 has a representative 30-minute light curve of 158,400 bytes, implying approximately 25 GB for 160,000 curves before intermediates. It is less comfortable than QLP for this machine. Do not request all pixel files: TESS-SPOC documents combined pixel, light-curve and DV volumes of hundreds of GB to TB per sector.

QLP sector 1 is substantially previously searched; searching it does not imply little prior study. The scientific value is a changed search followed by independent falsification: weaker/shorter events, alternative detrending, explicit multi-event support, broad sample coverage and recurrence in later sectors. Newer corrected sectors improve novelty but are significantly larger per target. Escalate selected targets to new sectors rather than attempting to cache every later-sector curve initially.

## Authoritative sources and verified archive paths

- [Official QLP HLSP documentation](https://archive.stsci.edu/hlsp/qlp), DOI `10.17909/t9-r086-e880`.
- [Official TESS-SPOC HLSP documentation](https://archive.stsci.edu/hlsp/tess-spoc), DOI `10.17909/t9-wpz1-8s54`.
- [Mission bulk-download entry point](https://archive.stsci.edu/tess/bulk_downloads.html).
- [Mission TCE-statistics download page](https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_tce.html).

Working target-list URLs retrieved during the probe:

| Target list | Rows indexed, not searched | Download size | URL |
|---|---:|---:|---|
| QLP sector 1 | 469,735 | 18,000,797 bytes | `https://archive.stsci.edu/hlsps/qlp/target_lists/s0001.csv` |
| QLP sector 99 | 2,149,786 | 98,607,483 bytes | `https://archive.stsci.edu/hlsps/qlp/target_lists/s0099.csv` |
| TESS-SPOC sector 83 | 159,997 | 8,118,429 bytes | `https://archive.stsci.edu/hlsps/tess-spoc/target_lists/s0083.csv` |

Lists contain TIC, RA and Dec, not all stellar parameters; rows cannot be described as scientifically vetted stars. Exact SHA-256 and retrieval URLs are in `/workspace/tess-research-data/broad_probe/catalog_manifest.json`. The lists are retained there, outside ordinary Git. Sector 99 was downloaded only to inspect archive scale and must not be the science default while its timing remains uncorrected.

Official shell manifests:

```text
https://archive.stsci.edu/hlsps/qlp/download_scripts/hlsp_qlp_tess_ffi_s0001_tess_v01_llc-fits.sh
https://archive.stsci.edu/hlsps/tess-spoc/download_scripts/hlsp_tess-spoc_tess_phot_s0083_tess_v1_dl-lc.sh
https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_99_lc.sh
```

The mission path requires `/missions/`; omitting that component returns 404. The QLP and TESS-SPOC manifests support HTTP `Range: bytes=0-4095` with status 206, useful for inspecting path formats without downloading scripts containing hundreds of thousands of commands. Do not execute downloaded shell scripts blindly; parse exact public product URIs and feed them to the cached downloader.

For TIC `8195168`, QLP sector 1 uses this verified URI:

```text
mast:HLSP/qlp/s0001/0000/0000/0819/5168/hlsp_qlp_tess_ffi_s0001-0000000008195168_tess_v01_llc.fits
```

Zero-pad the TIC to 16 digits and split it into four 4-digit directory components. The supported API URL is `https://mast.stsci.edu/api/v0.1/Download/file/?uri=` with the URI URL-encoded. TESS-SPOC uses the same directory split but inserts `target/`, and orders TIC before sector in the filename.

### Critical processing-version warnings

The official QLP page, checked during this probe, flags **sectors 99–104 with a barycentric timing error** and says corrected files are being redelivered. All six were still asterisked. Avoid those sectors pending authoritative correction status. QLP **sectors 74–79 require corrected v02 products**; v01 times can be wrong by approximately three minutes. Sector 94 onward uses TESS-Gaia Light Curve methodology, so comparisons across this boundary involve a reduction change.

The official TESS-SPOC page says sector 40–68 files lack CDPP values because of a processing bug. Missing CDPP does not imply infinite noise or a quiet star. TESS-SPOC FFI searches use a minimum two-hour pulse duration; independent searches of shorter events and older-sector multi-sector combinations are useful method differences. Its multi-sector FFI searches started with sectors 56–69; primary-mission and first-extended-mission FFI multi-sector searches are not planned by that team.

## Representative actual products and performance

Four retrieved raw observations are preserved under `/workspace/tess-research-data/broad_probe/products/`. No paid API was used. These are transfer/schema probes, not additional discovery results.

| Product | Bytes | Samples | Nominal cadence | QUALITY=0 and finite positive flux | Transfer latency |
|---|---:|---:|---:|---:|---:|
| QLP TIC 8195168, sector 1 | 80,640 | 1,282 | 1,800 s | 1,000 | 1.56 s |
| QLP TIC 8195216, sector 1 | 80,640 | 1,282 | 1,800 s | 1,000 | 1.56 s |
| TESS-SPOC TIC 196142994, sector 83 | 1,108,800 | 10,780 | 200 s | 10,462 | 2.24 s |
| TESS-SPOC TIC 950761231, sector 14 | 158,400 | 1,288 | 1,800 s | Not measured in this probe | 1.58 s |

The first three downloads ran concurrently. These are individual latency measurements, not a measured sustainable archive throughput. Benchmark a bounded concurrent download batch before estimating completion time; use connection reuse, limited retry/backoff and sensible concurrency. Preserve partial failures and completed downloads in an atomic resumable ledger.

The QLP FITS columns are `TIME`, `CADENCENO`, `SAP_FLUX`, `KSPSAP_FLUX`, `KSPSAP_FLUX_ERR`, `QUALITY`, `ORBITID`, `SAP_X`, `SAP_Y`, `SAP_BKG`, `SAP_BKG_ERR`, `KSPSAP_FLUX_SML`, `KSPSAP_FLUX_LAG`. Alternate apertures and background/centroids can help reject blends or systematics, but aperture agreement is not proof of an on-target planet.

QLP sector-1 sample times span BTJD 1325.329–1353.162, approximately 27.83 days. Headers explicitly state TDB and BJD reference 2457000. Both examples have 267 samples flagged 4096 and 15 flagged 4128; quality handling must be explicit. Stellar header Tmag/Teff are available, but radius/logg can be `-1`, requiring missing-value handling. Median-normalized robust flux RMS is approximately 1948 and 2087 ppm respectively; adjacent-sample scatter is approximately 1639 and 1877 ppm. These are example noise metrics, not survey population summaries.

The QLP examples have **no FITS CHECKSUM cards**. Retain transport TLS verification and record byte-level SHA-256, but do not claim intrinsic FITS checksum validation when absent. TESS-SPOC examples do have FITS CHECKSUM cards; reuse the repository's raw-byte verifier rather than altering expected checksums. Detailed hashes, schema/time metadata and per-quality counts are in the probe JSON files.

An Astropy fast-BLS timing probe on 1,000 good QLP samples used periods 0.5–13.5 days, durations 1/1.5/2/3/4 hours, oversample 5 and objective S/N. The simple median detrending used for timing is not a validated production treatment of data gaps.

| Frequency factor | Period grid | Seconds per star measured | Approximate wall time for 158,263 stars on five fully used CPUs |
|---|---:|---:|---:|
| 1 | 34,902 | 0.579 | 5.1 hours |
| 3 | 11,635 | 0.214 | 1.9 hours |
| 10 | 3,491 | 0.059 | 31 minutes |

These estimates exclude downloading, I/O, extraction, multiple detrenders and vetting, and assume the small benchmark represents the workload. Coarser grids reduce sensitivity: this example's maximum depth S/N fell from 9.72 to 9.16 at factor 10. The maximum was **not vetted as a planet**. A broad coarse pass is scientifically valid as exploratory screening only if its sensitivity limitations are explicit and tested; it cannot support a blanket statement that no planet exists in all searched stars.

## Executable survey hierarchy

1. Freeze a deterministic TIC sample from a corrected, manageable sector before looking at flux. Record the target-list checksum, sector/version, sampling seed/hash rule and all exclusions. Avoid first-N sorted TIC selection if claiming broad spatial sampling.
2. Begin actual flux screening immediately in checkpointed batches of at least 1,000 stars, then expand to tens of thousands and a 158k-scale subset if measured I/O and wall time permit. Indexing hundreds of thousands of catalog rows does not satisfy this step. Report downloaded products, unique stars with completed time-series searches, failed/quality-rejected stars and deeply vetted survivors separately.
3. Keep original FITS outside Git; retain aligned good-quality arrays or reconstruct them from exact product URIs. Hash results against data and method versions. Do not delete candidate data to make later batches fit. Stop expanding before the remaining disk required for candidate follow-up is consumed.
4. Use two scientifically reasonable gap-aware detrending alternatives or corrected QLP flux plus independent SAP local fits. Search finite durations appropriate to 30-minute integration and state that sub-cadence events are attenuated. Record more than one non-harmonic peak and actual event support. Single supported events must become a separate single-transit category rather than a periodic planetary claim.
5. Measure coarse-grid sensitivity using raw-flux injections at a frozen selection of periods, epochs, depths and durations. Include known controls and noise nulls distributed across star noise/activity bins. Infer no global false-alarm probability from uncorrected BLS S/N or per-star SDE.
6. Fine-search the strongest few percent with a denser period grid, multiple detrenders, local raw-SAP baselines and native cadence. Veto signals shared by many unrelated stars at identical spacecraft times, gap/sector/orbit boundaries, and stellar-variability harmonics. Keep rejection reasons in the signal ledger.
7. For genuinely recurring survivors, retrieve independent archival sectors and alternative reduction products, mask catalogued signals and search residuals, perform odd/even and secondary tests, physical transit fitting with uncertainties, contamination assessment and target-pixel/centroid localization. A deep eclipse or unmatched catalog period is not sufficient for planet promotion.
8. Refresh confirmed planets, TOI, CTOI, eclipsing-binary catalogs and literature for survivors using TIC, coordinate and aliases. Novelty scope/date and incomplete sources must be explicit. Archive-produced TCE tables can prioritize promising unmatched signals, but label this as re-vetting public search products, not a new analysis of every target the archive pipeline previously searched.

## Resumption and evidence preservation

The full lists and four schema products remain cached; no observations were discarded. Exact catalog manifest and local JSON schemas are retained in `broad_probe`. The parent survey should incorporate any reused observation into its durable URI/hash manifest. `/workspace` persistence remains unverified between tasks, so code, sample manifests, exact source identifiers and compact evidence must be committed. Large raw files remain outside ordinary Git. Re-fetching exact archived product URLs is the recovery path when a machine cache does not survive.

The survey must continue after the initial batch. A candidate-free thousand-star batch is a checkpoint for expanding coverage or changing a justified method, not a scientific reason to stop an otherwise feasible search.
