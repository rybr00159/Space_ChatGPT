# Broad-survey observation preservation

The raw broad-survey products are retained outside ordinary Git tracking in
`/workspace/tess-research-data/broad/qlp_s0001/raw/`. No independently durable
external volume has been verified. A future session must therefore check whether
these files survived before relying on the local cache. The committed manifest,
fixed target selection, scientific outputs, and source code provide a recovery
route when they did not survive.

`src/restore_broad.py` verifies or restores the **recorded observations** from
`results/broad/observations.jsonl.gz`. It does not count catalog entries as
observations, regenerate transit detections, or restore the SQLite checkpoint.
The observation export is an atomic snapshot; while the survey is running it can
lag the live database. Recovery uses whichever manifest snapshot was committed.
Any selected targets not yet present in that snapshot remain the responsibility
of the checkpointed survey, rather than this restoration utility.

## Manifest schema and integrity

Each gzip JSONL record contains `tic`, `ra`, `dec`, the original public MAST QLP
HTTPS `url`, the original absolute `path`, expected `sha256`, `raw_bytes`,
`fits_checksums`, FITS `metadata`, and `checked_utc`. The current broad products
are sector-1 QLP v01 light curves; their `metadata` includes `TICID`, `SECTOR`,
`TESSMAG`, `TEFF`, `RADIUS`, `PROCVER`, `DATE`, `CAMERA`, and `CCD`.

The utility uses the recorded SHA256 and byte count as immutable expectations.
It remaps the old absolute path into the chosen cache using the product filename,
sector, and TIC. This produces the same directory layout as `broad_survey.py` and
allows recovery onto another machine. It verifies FITS structure, raw
ones-complement HDU `CHECKSUM` and data `DATASUM` when present, TIC/sector identity,
TDB time system, and the BTJD reference of 2457000. QLP products in this survey
have no intrinsic FITS checksum cards; their per-file SHA256 still identifies
the exact retained bytes. An absent FITS card is reported as absent.

The script hashes and reads the same open manifest descriptor so a concurrent
atomic export cannot change its input midway through a run. Each run saves its
manifest hash, source hash, settings, per-observation results, and summary to a
unique file under `<cache>/restoration/`. These logs preserve verification
history without overwriting an earlier run.

## Commands

Use the project Python environment after installing the pinned dependencies.
Default operation performs verification only and makes no HTTP requests:

```bash
python src/restore_broad.py --workers 4
```

Inspect a small sample or one target:

```bash
python src/restore_broad.py --limit 20 --workers 2
python src/restore_broad.py --tic 370247818
```

Restore at most 1,000 missing observations per run to a chosen cache, retaining
at least 5 GiB free space before starting each transfer:

```bash
python src/restore_broad.py \
  --cache /workspace/tess-research-data \
  --download --max-downloads 1000 --workers 4 --reserve-gib 5
```

Both `--download` and an explicit positive `--max-downloads` are required to
enable transfers. The cap counts attempted transfers, including failures.
Workers are bounded between 1 and 16, with at most twice that many queued tasks.
Concurrent disk reservations include the recorded file sizes. Repeating the
same command resumes at the observation level: already verified products use
no download allowance, while missing products can fill the next batch. Partial
HTTP bodies are retained for inspection; they are not resumed by HTTP byte
range. No archive credential or paid service is needed.

An existing file that differs from the recorded size/hash is preserved and
reported as a failure. The utility never repairs it by overwriting it. A newly
downloaded body is accepted into the raw cache only after all verification has
passed. Changed archive content and incomplete bodies remain at a unique
`.restore-<run-id>.part` path. Installation uses an atomic non-overwriting hard
link, so a simultaneous run cannot replace an existing destination. Removal of
a verified temporary hard link leaves the same verified scientific bytes at
their final path.

Exit code 0 means every selected observation was verified. Exit code 2 means
the selection remains incomplete, for example because files are missing, the
download/storage limit was reached, the selection is empty, or a requested TIC
was not processed. Exit code 1 indicates malformed input, an integrity failure,
or a transfer/execution error. Read the JSON summary and evidence log before
continuing analysis. A requested TIC absent from the committed broad snapshot
cannot be restored with this manifest; use its relevant follow-up manifest when
one exists.

After restoring the complete broad manifest, run verification again without
`--download`. To reconstruct a lost survey database, rerun the preserved
`src/broad_survey.py` with the fixed selection and frozen analysis version from
the corresponding commit; the recovered raw files will be reused. Keep the
old committed results and their manifest alongside any reconstruction. Do not
substitute a changed archive product or changed analysis version and label it
the original experiment.

## Verification performed on 2026-10-14 UTC

The final script verified 20 retained real QLP light curves with two workers:
20 `verified_cached`, zero HTTP attempts, zero integrity failures. The checked
manifest snapshot SHA256 was
`73338b46f8831348132b4e9e0ac664735bab13bb88a8514eae2940c6da84676e`.
The per-run evidence is retained locally in
`/workspace/tess-research-data/restoration/20261014T050317-584305af.jsonl` and
its sibling `.summary.json`.

Additional temporary offline fixtures verified that verification mode cannot
construct a network session, a cap of one admits only one of three missing
products, an existing corrupt file is left unchanged, a mismatching downloaded
body cannot enter the scientific cache, both raw `CHECKSUM` and `DATASUM` are
verified on a synthetic FITS product, and an empty/missing-TIC selection returns
an incomplete status. These are targeted integrity checks rather than a claim
that the entire broad raw cache has been audited.

## Full scientific checkpoint, independent of raw-cache survival

The completed 160,000-row SQLite state is also preserved in `data/broad_checkpoint/`: four checksum-verified gzip JSONL chunks, approximately 78.1 MB total, with each chunk below 20 MB. This stores the original result JSON strings, statuses, retrieval/source metadata, timestamps and method identity exactly. It avoids rerunning completed searches if a future task loses the local database. This is scientific analysis state; the 12.9 GB raw observations remain outside ordinary Git. The separate single-file compressed backup also remains locally retained.

```bash
python src/survey_checkpoint.py verify
python src/survey_checkpoint.py restore --db /workspace/tess-research-data/broad/qlp_s0001/survey.sqlite
```

Restore requires an unused destination and never replaces an existing database. Raw-product verification/restoration is separate: retained result paths do not establish that the raw files survived. An actual round-trip into a separate SQLite file passed the SQLite integrity check, identical status counts and metadata, and exact SHA-256 equality of all 160,000 original sorted records. Evidence is in `results/broad/checkpoint_roundtrip.json`. Restoration was tested within this machine; new-cloud-task restoration and publication remain unverified.
