# New-sector archive reconnaissance and pre-flux selection

Checked October 15, 2026 (America/Halifax; retrieval manifests use UTC October 16).
This report describes archive coverage and the target-selection input, **not a
completed transit search of the inventory stars**. Eight randomly chosen raw
products were inspected to establish the current data interface.

## Chosen sectors

| Sector | Official QLP inventory rows | Ecliptic hemisphere | Raw probe size | Raw probe rows | Quality=0 rows in the two probes |
| --- | ---: | --- | ---: | ---: | --- |
| 83 | 1,127,702 | North | 639,360 bytes | 10,505 each | 10,328; 10,328 |
| 86 | 977,790 | North | 682,560 bytes | 11,207 each | 5,763; 5,553 |
| 95 | 497,114 | South | 648,000 bytes | 10,626 each | 7,816; 7,433 |
| 98 | 628,730 | South | 1,460,160 bytes | 24,162 each | 12,783; 11,974 |

All four complete target lists returned HTTP 200 and are retained, without
modifying the old Sector 1 inventory. Their columns are **only `#TIC_ID`, `RA`,
and `DEC`**. They cannot support a brightness or stellar-radius selection by
themselves. Inventory counts overlap and must not be summed as unique stars.

The hemisphere characterization uses a fixed random sample of 2,000 coordinates
per sector, transformed with Astropy to barycentric true ecliptic coordinates;
all sampled Sector 83/86 targets were northern and all sampled Sector 95/98
targets southern. The selected fields cover different longitudes with useful
overlap near the continuous viewing zones. In the eligible metadata sets below,
83/86 share 236 TIC IDs and 95/98 share 350; north/south sets share zero. This
provides both fresh sky coverage and repeat measurements without labeling
repeat observations as new stars.

## Current official timing warning

The [official QLP page](https://archive.stsci.edu/hlsp/qlp) currently warns that
barycentric timestamps in **Sectors 99–104** were affected by ephemeris and
interpolation problems. The page still marks those six sectors with asterisks
and says their corrected products have not yet been redelivered. We therefore
exclude those sectors from the first blind-search tranche. A file declaring
`TIMESYS=TDB` does not by itself prove this known problem is corrected.

Sectors **74–79** had a different ephemeris reference-plane error, up to about
±3 minutes, now redelivered as **v02** products. Their v01 products should not be
mixed into a precision timing search. The four initial sectors (83, 86, 95, 98)
avoid both warning ranges. This is a current official-status check, not an
independent recalculation of all barycentric timestamps.

The exact HTML and readme responses, URLs, retrieval timestamps and SHA256
hashes are retained. Their compressed snapshot is committed with the metadata
query archive. The older official readme primarily documents the earlier data
releases, so the actual downloaded FITS files establish the new schema here.

## Actual FITS interface and storage implications

All eight probes returned HTTP 200, parsed as FITS and passed Astropy structural
verification. Their data extensions contain:

```
TIME CADENCENO SAP_FLUX DET_FLUX DET_FLUX_ERR QUALITY ORBITID
SAP_X SAP_Y SAP_BKG SAP_BKG_ERR DET_FLUX_SML DET_FLUX_LAG SYS_RM_FLUX
```

The usable corrected flux interface is `DET_FLUX` / `DET_FLUX_ERR`, rather than
the old `KSPSAP_FLUX` naming. `SAP_FLUX` and alternate-aperture flux columns are
retained for later checks. All have `TIMESYS=TDB`, `BJDREFI=2457000`, and a native
cadence near 200 seconds. `BJDREFF` is absent and contributes zero. `TIMEDEL`
equals 200/86400 days. There is no `PROCVER` header in these probes; the archive
filename records `v01`, not a complete internal pipeline version.
The primary header separately records the upstream calibration as
`CALIB='TICA v1.4.0'` and `TICVER=8.2` in these probes.

**No probe contains FITS `CHECKSUM` or `DATASUM` cards.** Astropy's corresponding
return value 2 means absent, not a verified checksum. Local SHA256 measurements
preserve the exact downloaded bytes. Producer checksums cannot be claimed.

Sector 98 probe timestamps span about 57.4 days despite many flagged/gapped
points. Sector numbers do not guarantee a uniform 27-day continuous series.
Quality masking and explicit segment boundaries are necessary, and live time
must be calculated from retained cadences rather than `TSTOP-TSTART` alone.
The Sector 98 probes contain orbit identifiers 205–208, preceding the current
official timing warning's start in Orbit 210; those identifiers are recorded
explicitly with cadence counts and ranges in the probe JSON.
For Sector 86, roughly half of the sampled cadences have nonzero quality flags.
These examples do not determine the distribution for all targets in a sector.

A 1,024-target-per-sector first tranche would retain approximately 3.5 GB of raw
products at these example file sizes. Downloading every listed target in the
four sectors would be a terabyte-scale task, exceeding this environment's
ordinary storage budget. Sector-by-sector checkpointing, explicit target
selection and candidate preservation remain necessary.

## Fixed pre-flux metadata sample

`src/multisector_inventory.py` takes the first 8,192 IDs in each complete target
inventory ordered by ascending
`SHA256('multisector-v1:' + decimal_TIC_ID)`. The same salt in each sector
encourages useful repeat coverage. No light curve, transit score, variability
measure or known candidate ephemeris enters this selection.

Exact-ID public TIC requests use `Mast.Catalogs.Filtered.Tic`, 512 IDs per request,
request columns `ID,ra,dec,Tmag,Teff,rad,contratio,lumclass`, page size 10,000.
The service returned full TIC rows despite the requested projection; all exact
responses are preserved. All 32,768 requested sector/target metadata rows were
returned without missing IDs. Repeated stars across sectors are not distinct
newly searched stars. The TIC `version` field in a sampled response is
`20190415`, whereas the HTTP retrieval date is October 2026; those dates have
different meanings.

The initial eligibility filter requires finite `Tmag<=12.5`,
`3000<=Teff<=6500 K`, and `0.15<=rad<=2.5` solar radii. Stellar parameters and
contamination are catalogue estimates and need survivor-specific verification.
Missing/outside-range entries remain in the full metadata CSV, with explicit
rejection reasons; this filter is not proof those objects cannot host planets.

| Sector | Pre-flux metadata rows | Initially eligible |
| --- | ---: | ---: |
| 83 | 8,192 | 1,026 |
| 86 | 8,192 | 1,043 |
| 95 | 8,192 | 1,761 |
| 98 | 8,192 | 1,535 |

The 5,365 eligible sector/target rows correspond to **4,779 distinct TIC IDs**.
This is target availability, not the number downloaded or period searched.
`data/multisector/inventory_sNNNN_eligible.csv` provides the pre-flux ordered
pool for the survey runner; `inventory_sNNNN_metadata.csv` retains every chosen
metadata row, including rejects. `sourcehash` identifies the exact TIC body;
`inventory_sourcehash` identifies the complete official target list.

## Reproduction and preserved evidence

```bash
OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python \
  src/multisector_inventory.py --sectors 83 86 95 98 \
  --size 8192 --batch-size 512 --workers 6
```

Existing response bodies are validated against their saved hashes and reused.
The normalized CSVs should reproduce exactly from the cached responses; a new
live catalogue query may legitimately change and must be a separately recorded
snapshot. Full official target lists and raw FITS remain outside ordinary Git
tracking under `/workspace/tess-research-data/multisector/`. These files have
not been deleted. Persistence of the workspace between future tasks is not
guaranteed.

`data/multisector/archive_probe.json` records every raw-product probe URL,
path, size, SHA256 and header/schema measurement. The adjacent geometry and
overlap JSONs record their sampling/counting scope. The compressed query
snapshot and its manifest preserve exact TIC responses plus the official
warning page/readme. Inventory URLs and source hashes support fresh retrieval;
if an archive product changes, preserve the changed response separately.

The exact query snapshot is 10,410,576 bytes in xz format. The checkpoint helper
checks all 132 archived response/manifest members and the 14 source/normalized
artifact hashes. To verify it or restore the original query bytes without
network access:

```bash
/workspace/openrouter-check-venv/bin/python \
  scripts/multisector_inventory_checkpoint.py verify
/workspace/openrouter-check-venv/bin/python \
  scripts/multisector_inventory_checkpoint.py restore-queries
/workspace/openrouter-check-venv/bin/python \
  scripts/multisector_inventory_checkpoint.py download-inventories
```

`restore-queries --cache-root /tmp/multisector_query_roundtrip` was executed as
an offline restoration test. It reproduced all 132 members byte-for-byte and
verified their recorded hashes. `download-inventories` checks the recorded
immutable source hashes, refuses to overwrite changed existing files, and saves
changed live responses separately instead of accepting them as original data.

The next scientific operation is an actual light-curve search on a fixed subset
of these eligible pools, followed by cross-sector recurrence checks and
injection-based sensitivity measurements. None of the inventory or metadata
counts above is a detection claim.
