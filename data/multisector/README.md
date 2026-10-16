# Data for the six-sector research checkpoint

The survey samples Sectors83/86/95/96/98/105. `selection.csv` and `b02_096_105/selection.csv` are immutable pre-flux target lists; metadata and their query provenance are retained. Read `../../MULTISECTOR_RESEARCH_REPORT.md` for actual deduplicated search counts and scientific limitations.

Large raw/cache data remain outside ordinaryGit at `/workspace/tess-research-data/multisector/`. The last raw inventory records5152FITS+8ZIPfiles,4510156932bytes, each with an existing expected SHA256and public URL. External durability of this cloud filesystem is unverified. No datasets were discarded to reduce repository size.

Tracked preservation products:

| Path | Purpose |
|---|---|
| `checkpoint/` | Consistent main SQLite backup:4096completed records, including3545producer-error searches. Exact original download manifest. Failed WAL backup identified separately and preserved in cache. |
| `archive_probe_query_archive.json` | Manifest for132exact original metadata query responses and their compressed archive. |
| `latest_probe/` | New-sector availability, timing/version evidence and exact metadata-query archive. |
| `supplemental_checkpoint/` | Consistent551/963/689record DBs;667completed pooled-star JSONfiles; exact small follow-up responses; recorded-input and raw inventories. |
| `candidate_checkpoint/` | TIC229146531's six original QLP light curves plus native pixel cube:7exact FITS,7.85MBcompressed, tested offline restore. |
| `catalogs/` | Immutable public catalogue/alias/literature responses and manifests. Nonmatches have finite release/query scope. |
| `diagnostic_preservation.json` | Original/retained paths and hashes for bulk/curated diagnosticPNG/NPZfiles.376bulk files,123011858bytes, moved intact outsideGit; important arrays/figures remain tracked. |

Use the verified restore commands in `../../PROGRESS.md`. Supplemental restoration preserves completed analysis without rerunning thousands of searches; prepared curves and figures can be regenerated with frozen code as needed. Static QLPobservations can be retrieved in bounded batches with `scripts/restore_multisector_raw.py` using immutable hashes. Existing mismatches are never overwritten; mismatching fresh responses are retained separately.

TESSCut containers can change creation/header bytes when regenerated even from the same FFIs. The priority target's original cube is therefore archived directly. Other native cube manifests retain retrieval URLs, timing, original byte/science hashes and source FFIidentifiers; a regenerated container must not silently replace expected originals. Original ZIPs remain in local storage rather than being duplicated alongside their extracted cube inGit.
