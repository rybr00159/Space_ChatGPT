# Catalog recognition for the multi-sector search

Public core catalogs were retrieved on 16 October 2026 UTC (evening 15 October
in America/Halifax). The immutable snapshot is
`20261016T010726453469Z/`; `latest.json` identifies it explicitly. No earlier
catalog file or source method was changed. These are recognition checks, not
light-curve analyses and not proof of novelty.

| Source | Rows | Snapshot scope |
|---|---:|---|
| NASA confirmed planets (`pscomppars`) | 6,451 | Fresh; all discovery methods, including RV-only planets |
| NASA TOI | 8,151 | Fresh; all dispositions |
| ExoFOP CTOI | 5,138 | Fresh; all community reports and dispositions |
| DIAmante v2 | 1,160 | Reused immutable 2026-10-13 published catalog |
| Prša et al. TESS EB | 4,584 | Reused immutable published sectors 1–26 catalog; not all current EBs |
| Olmschenk 2021 and Yu 2019 | 473 | Reused immutable published candidate tables |
| RAVEN release 19661443: NSFP > 0.9 | 14,815 | Versioned public release; reused verified 2026-10-14 download |
| RAVEN vetted table | 2,170 | Same release; overlaps other RAVEN tables |
| RAVEN validated table | 143 | Same release; not independently validated by this project |
| RAVEN large-radius table | 207 | Same release; overlaps other RAVEN tables |
| User-supplied TIC 4206066 manuscript | 2 signals | Reported 3.182785-day candidate and tentative 11.13274-day signal |

The combined index contains 43,294 rows and 29,261 supplied TIC identifiers.
Rows, targets and confirmed planets are different quantities. Overlapping
RAVEN tables deliberately preserve the separate provenance of each report;
their sum is not a count of independent candidates.

Each fresh HTTP response has its original bytes, retrieval time, URL, HTTP
status, SHA-256, row count and exact column schema recorded. Parsing rejects
schema changes, HTML/errors, malformed CSV and unexpectedly short responses.
Incomplete snapshots remain preserved but cannot become the active snapshot.
Published support tables retain their original manifests and dates. The RAVEN
release record's immutable SHA-256 and each file's release MD5 are verified;
earlier failed 91-byte responses are never accepted as catalog tables.

The supplied manuscript is retained outside Git and its SHA-256 is recorded
in `candidates/tic4206066/recovery.json`. It is included explicitly because
neither reported signal appears in the checked bulk online catalogs. A
rediscovered signal at that TIC/period cannot be claimed novel.

## Reproduction

Use the configured Python environment and run from the repository:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python src/multisector_recognition.py
python src/multisector_recognition.py --controls
python src/multisector_recognition.py \
  --input results/multisector/signals.csv \
  --output results/multisector_recognition/signals.csv \
  --rank-unknown results/multisector_recognition/prioritized_unknown.csv
# New immutable public core snapshot, only when a refresh is wanted:
python src/multisector_recognition.py --retrieve
# A reproducible run can pin --snapshot data/multisector/catalogs/20261016T010726453469Z
```

Import API with `src/` on `PYTHONPATH`:

```python
from multisector_recognition import FreshCatalogIndex
index = FreshCatalogIndex()
flags = index.match_signal(tic, period_days, ra=ra, dec=dec)
searchable = index.transit_mask(tic, times_btjd, ra=ra, dec=dec)
```

Input columns include `tic`, `period_days`, and optional `ra`/`dec`. The
operational ranking additionally uses the survey's widths, event counts and
local metrics. Annotation preserves the exact input CSV as a compressed file,
records input/output and source-code hashes, and writes distinct derived
output. Positive controls are in
`results/multisector_recognition/positive_controls.json`.

## Interpretation and limitations

An exact TIC or <=2-arcsecond coordinate match is a provisional host match;
other reports within 30 arcseconds are possible contaminants. Coordinates are
not propagated for proper motion. Within 1% of a reported period is a coarse
recognition flag; integer harmonics/subharmonics 2–5 are retained separately.
These flags do not prove that events or sources are identical. Matching
epochs, raw events, source localization and literature remain necessary.

Known hosts with a different period remain eligible for additional-planet
searches. Nontransiting RV planets are included even when they have no transit
epoch. Missing epochs are never invented for transit masks. Some included
RAVEN detections are not vetted planet candidates; their membership establishes
a prior public report, not a physical classification.

`no_report_in_checked_snapshots_pending_literature` means exactly that. It
does not establish novelty. Coverage does not include a complete all-sector
SPOC/QLP TCE/DV inventory, every later EB, all citizen-science reports, all
aliases or all published literature. The optional rank is a follow-up queue;
agreement across two detrending widths uses the same underlying photons and
is neither independent validation nor a calibrated false-alarm probability.

`--live-tic TIC` additionally saves a public ExoFOP target JSON, a 30-arcsecond
SIMBAD cone, a 30-arcsecond VSX cone, and an exact-TIC-text arXiv API query.
These small lead-only checks retain errors explicitly. Empty results from a
failed query are never interpreted as absence. Live target pages may expose
prior files and notes missing from the bulk tables; literature searches are
finite and incomplete. All original responses remain in timestamped `live/`
subdirectories. Derived positions/separations and finite-check results are in
`results/multisector_recognition/live_context_summary.json`.

After a live target check, `--deep-tic TIC` retrieves Gaia DR3 properties,
60-arcsecond neighbors/classifications, actual target NSS and selected
variability tables, an ASAS-SN published 2021 cone, and a 60-arcsecond SIMBAD
alias/bibliography query. These retain source names and coordinates so that
nearby objects cannot inherit the target's classification. Existing successful
responses are reused only after URL and SHA-256 verification; changed or failed
responses are never silently accepted. The derived context and finite-check
limitations are in `results/multisector_recognition/deep_context_summary.json`.
