# Catalog and source checks for TIC 71316387

Public sources retrieved on 2026-10-14 UTC confirm that **CTOI 71316387.01 is a RAVEN high-probability planet candidate, not a validated planet**. The current ExoFOP disposition is `PC`, with the upload dated 2026-08-27 and attributed to Lafarga et al. (2026). The relevant paper is [arXiv:2603.22597v2](https://arxiv.org/abs/2603.22597v2), with [Zenodo tables, record 19661443](https://zenodo.org/records/19661443).

| Quantity | Published value |
|---|---:|
| Fitted candidate period | 0.6518369465110034 d |
| Fitted epoch | BJD 2459067.327192476 |
| Fitted duration | 1.462566491476932 h |
| ExoFOP depth parameter | 5494.622875131619 ± 1338.466290297098 ppm |
| Published radius | 8.544629706471323 Earth radii |
| RAVEN score | 0.9309549774171888 |
| Mean planet versus non-simulated false-positive classifier probability | 0.9986432 |
| Number of full events in RAVEN table | 67 |
| Sectors in RAVEN target summary | 1 and 28 |

The exact target row is in `vet_table_retry.csv` and in the `raven_target_rows` field of the most recent `catalog_check_*.json`. This target does not appear in `val_table.csv` or `val_rp8_table_retry.csv`. The latter is the separate sample with RAVEN probability >0.99 and inferred planet radius >8 Earth radii.

RAVEN defines its final score as the **minimum of separate planet-versus-false-positive posterior probabilities**, rather than the joint posterior against all hypotheses. For this target the limiting comparison is the hierarchical eclipsing binary scenario (0.9309549774). The planet-versus-target eclipsing binary comparison is 0.9866898498; the background binary comparison is 0.9997331092. These model-dependent scores do not settle a contradiction found in original photometry.

The [RAVEN methodology paper](https://arxiv.org/abs/2509.17645v1), section 2.5.2 and Table 5, explicitly includes odd/even eclipse depth and duration ratios from separate trapezoid fits as classifier features. Its synthetic-event recovery also considers half/double-period aliases. The survey paper's section 3.4 describes additional visual checks of harmonic aliases. An independently demonstrated alternating-depth signal would therefore warrant reviewing this particular classification; it would not justify claiming that RAVEN lacks an odd/even diagnostic.

The ExoFOP target summary PDF (`raven_summary.pdf`, public file ID 1113786) supplies the sector list and scores. Its fitted depth is a model parameter from a single-period transit fit; it should not be directly equated with a median in-eclipse flux depression, especially when two event families or integration-time effects are under consideration. Sector 1 was observed at 30-minute cadence and sector 28 at 10-minute cadence. The target's public Villanova availability page also lists sector 68, which was outside the survey's sectors 1–55 and was not used in the published target summary.

## Finite independent catalog context

Current Gaia DR3, SIMBAD, VSX and selected arXiv aliases were checked. Gaia source ID is `6571343840265900544`. Neither Gaia's eclipsing-binary table nor its variability-classifier table returned a target row. None of the 22 Gaia sources in the 120-arcsecond cone has a Gaia eclipsing-binary summary flag. The VSX cone returned no objects. SIMBAD lists the target as the generic star `TYC 8001-154-1`, with zero linked references. A metadata search of arXiv for TIC number, selected 2MASS/TYC aliases and the Gaia source ID returned no entries. These are limited searches: papers' machine-readable tables, follow-up databases and unindexed reports need not appear in arXiv metadata. Absence here cannot establish a previously unreported binary.

Gaia RUWE is 1.0265, while its image-parameter determination reports multiple peaks in 41% of successful observations. The latter is context for possible image multiplicity or blending, not confirmation of a close binary. Gaia's reported robust RV amplitude is 9.93 km/s, but its RV constancy chi-square p-value is 0.275; the amplitude alone is not significant evidence for orbital RV variability. No Gaia non-single-star solution is flagged.

The closest resolved ExoFOP neighbor is TIC 71316388, 18.19 arcseconds away, TESS magnitude 16.0571. TIC 71316389 is 46.66 arcseconds away at magnitude 14.667; TIC 71316386 is 55.41 arcseconds away at magnitude 14.737. These sources can contribute to TESS apertures. Original pixel data, rather than catalog positional probabilities alone, are needed to attribute the dips.

## Reproduction and retention

From the repository root:

```bash
/workspace/openrouter-check-venv/bin/python scripts/tic71316387_catalog.py --large-papers
```

Existing responses are reused after SHA-256 verification; `--refresh` creates additional dated snapshots. Successful Zenodo table responses are additionally verified against the record's MD5 checksums. Exact public responses, including initial 403/503 errors and subsequent successful retries, are retained. Sidecar manifests record URL, retrieval time, HTTP status, byte count and SHA-256. The most recent `catalog_check_*.json` has no retrieval failures. Original paper PDF and source archive, each about 27 MB, remain at `/workspace/tess-research-data/tic71316387/catalog/` outside Git; their manifests and extracted paper text are tracked. This cache's survival across a fresh execution environment has not been verified. The script can retrieve the files again from their recorded public URLs.

The independent photometric validation is documented separately in the parent candidate directory. These catalog checks do not themselves verify either the quoted alternating depths or their significance.
