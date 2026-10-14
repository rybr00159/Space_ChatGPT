# Stellar companion scout: retained QLP observations

The new stellar search reuses the actual flux series already retained from the
fixed 160,000-star sector-1 sample. It does not count catalog rows as newly
analyzed observations and does not change the frozen planetary search.

## Scope and method

The original 319,496 stored peak rows were reconsidered with an EB-oriented
selection. A candidate needed an unmasked depth of 3,000–500,000 ppm, diagonal
BLS S/N ≥15, local S/N ≥10, at least three sampled events, at least 70% positive
local depths, duty cycle below 0.3, and no single-event-domination flag. Wide
eclipses rejected by the planetary morphology guard were allowed. Cached Prša
2022 EB hosts and confirmed planet hosts were excluded; TOI/CTOI candidates were
retained with a modest priority penalty because revised binary classifications
are allowed by the new mission.

These rules select **1,253 distinct stars**. Before inspecting new native-flux
results, the pilot was fixed to the highest-ranked 144 moderate-depth
(3,000–50,000 ppm) and 48 deep signals. **192 native SAP light curves were then
independently remeasured**, with exact raw-file SHA-256 verification. The selected
sample is deliberately enriched in obvious eclipses and cannot measure the
incidence of previously unknown binaries.

Each candidate's original period was locally refined ±1% using QLP flux and a
slow running-median baseline. The native SAP series was then measured with
independent local quadratic baselines, empirical per-event scatter floors,
primary and half-phase eclipse windows, alternate-event groups at 2P, and three
QLP aperture products. The retained `broad_vet.local_quadratic` and variability
functions are imported without changing their frozen source. All formal
contrasts are conditional on a selected period and shape; they are **not global
false-alarm probabilities**. A half-phase search misses eccentric secondaries.
Approximate box depths differ from central transit-model depths.

## Pilot results

- **153/192** have ≥3 positive repeated raw SAP events, raw S/N ≥10, and raw
  eclipse depth ≥3,000 ppm.
- **40** of those also have ≥3 measured half-phase events and conditional
  secondary S/N ≥5.
- **48/192** show an alternate-event contrast ≥5 times its conditional estimated
  error. Some are already known binaries, and cadence/shape tests remain needed
  before interpreting their precise contrasts.
- All 192 retain complete measurements and failed checks in per-target JSON;
  the sorted searchable summary is `results/stellar_systems/scout_raw_measurements.csv`.

Expanded full-table prior checks used Shi/Qian/Li (2022, **1,626** EA entries),
Tang et al. (2025, **514** eccentric EBs, CDS released 2026-09-14), and
Battley et al. (2022, **9,400 YOUNGSTER** candidate EA/EB signals). YOUNGSTER is
the survey name: it must not be mislabeled as CDIPS. Raw retrieved TSV bytes,
URLs, retrieval times and SHA-256 manifests are retained in `scout_catalogs/`.
The published YOUNGSTER entries are candidates rather than uniformly confirmed
binary systems. Exact TIC matches and all neighbors within 120 arcsec are stored
in `scout_expanded_catalog_matches.json`. A nearby entry is a contamination
warning, not a host association.

A targeted 120-arcsec VSX, ASAS-SN and SB9 check of eight high-ranked deep signals
found six to be known binary signals or known nearby binary contamination:

| TIC | Result of targeted prior check |
| --- | --- |
| 260503860 | ASASSN-V J062601.94-554111.5, EA, P=1.7113 d |
| 115188534 | ASAS J210932-3958.1, ED, P=3.485592 d; scout chose P/2 |
| 300652458 | ASAS J073721-6729.1, ED, P=3.2586 d; scout chose P/2 |
| 115152243 | ASASSN-V J210752.55-404004.3, EA, P=4.7097 d; scout chose P/2 |
| 294206969 | ASASSN-V J071313.06-591848.8, EA, P=2.758 d |
| 89018908 | Nearby known DW Gru / TIC 89018905, EA, P=3.14283 d |
| 52285789 | No hit in these queried catalogs; 25.8% eclipse, Pphot≈3.30356 d, unequal alternating depths suggest ≈6.60712 d |
| 403150111 | No hit in these queried catalogs; P≈3.91134 d, primary≈6.94%, half-phase eclipse≈0.912% |

Subsequent recognition **rejects TIC52285789 as a new companion: it is a published GALAH SB2. TIC403150111 is already a recognized visual pair**, with any eclipsing subsystem still unclassified. The table above records the limited first-pass checks, not the final novelty assessment.
They need current Gaia variable/NSS and live TESS EB checks, source localization,
independent-sector confirmation, aliases and literature searches. Several
moderate-depth targets receive those checks in parallel; parent reports record
their resulting disposition. This scout alone establishes no discovery.

## Estimator controls and reproducibility

Synthetic cadence-integrated equal/unequal eclipse trains were injected onto
the actual quality-filtered timestamps of TIC 294272014, 38825979 and 281581467,
with a slow baseline and 500-ppm independent Gaussian noise. Twenty noise draws
per model per target give 60 equal-depth and 60 unequal-depth controls. The largest
equal-depth parity contrast was 2.26 times its estimated error; the weakest
10,000-ppm injected contrast was detected at 7.09 times its estimated error.
These controls check the local estimator on real cadence patterns. They do not
emulate all instrumental noise or calibrate the full 160,000-star search.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python src/stellar_scout.py --limit 192 --deep-limit 48
python reports/stellar_systems/scout_controls.py
python reports/stellar_systems/scout_catalog_audit.py
```

The source hash and catalog provenance are in `scout_manifest.json`. Native raw
products remain outside ordinary Git tracking. Their download URLs and SHA-256
are in each scout JSON and the original broad observation manifests. No raw
file or previous scientific result was deleted. Every completed target is saved
atomically, so rerunning the same frozen scout resumes rather than repeats it.
Paid-service expenditure for this phase is **$0**.

The final expanded recognition results are [scout_novelty.md](scout_novelty.md). Both vetted TESS Ten Thousand tables (7,936 new plus2,065 known EBs) and the complete Traven2020 SB2 catalog are now included in the reranking. TIC197764088 was independently localized to knownneighbor197764089. The first64 fresh reserve rows already have native scout measurements but need their own current live recognition review; they are not new candidates by catalog absence.
