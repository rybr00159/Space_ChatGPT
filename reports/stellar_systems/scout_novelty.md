# Stellar scout recognition audit

Checked 2026-10-14T20:37:06.692046+00:00 (UTC). Eleven selected repeatable eclipse-like scout signals were investigated against multiple source catalogs. **Ten have a positively recognized binary host or a nearby known eclipsing binary with a matching period/harmonic. The remaining target, TIC403150111, is already a recognized visual double, with an eclipsing subsystem not yet fully classified. These results do not establish a previously unknown star.** Nearby matches are contamination hypotheses, not pixel-based source assignments.

The important falsifications were TIC38825979 and TIC52285789: they had no match in the original eclipsing-binary checks, but both are published double-lined spectroscopic binaries in [Traven et al.2020](https://doi.org/10.1051/0004-6361/202037484). TIC197764088 has a matching known EB, TIC197764089, 63.473arcsec away in [Kostov et al.2025](https://doi.org/10.3847/1538-4365/ade2d8). TIC394136010 is already recognized from [Nardiello et al.2019 PATHOS-I](https://doi.org/10.1093/mnras/stz2878). These cases show why absence from one EB catalog cannot establish a new companion.

| Scout TIC | Scout period (d) | Prior recognition | Assessment |
| --- | ---: | --- | --- |
| 294272014 | 0.521224 | Gaia DR3 5486571994791177344 / VSX3356349 | Known eclipsing star 2.33arcsec away at approximately half the detected period; TESS source attribution remains unresolved. |
| 355544861 | 0.492717 | Gaia DR3 6463554115486620544 / VSX7940975 | Known nearby eclipsing star with the same period; probable contaminant, requiring pixels for attribution. |
| 38825979 | 2.520145 | GALAH160124002601234 / Gaia4675527035087606144 | Known double-lined spectroscopic binary from Traven2020. TESS eclipsing characterization may add information, but a second star is already recognized. |
| 141866704 | 1.423306 | TIC141866705 / ASASSN-V J061547.24-751510.2 | Known eclipsing neighbor with same or twice-period ambiguity; candidate target CTOI is a separate existing recognition. |
| 141411137 | 1.757640 | Gaia DR3 4650316195543230592 / VSX4313616 | Host is already a recognized eclipsing binary in Gaia and VSX. |
| 281581467 | 1.995030 | ASASSN-V J074356.39-602347.2 / VSX838559 | Host classified EA+BY in VSX; catalog binary period approximately twice the scout's detected period. |
| 388181335 | 1.588086 | TIC388181338 / ASASSN-V J070310.13-735512.1 | Known eclipsing neighbor with matching period; probable contaminant pending pixels. |
| 394136010 | 1.406321 | UCAC4 096-000233 / PATHOS-I2019 | SIMBAD positively recognizes host as EB* from Nardiello2019, although several other EB catalogs have no match. |
| 197764088 | 2.682059 | TIC197764089 / TYC8820-40-1 | Kostov2025 vetted EB nearby has matching period and41/16ppt depths. Host Gaia RV is stable, whereas the neighbor is strongly variable; likely known-star contamination. |
| 52285789 | 3.303560 | GALAH170911004701327 / Gaia4692405049673289472 | Known double-lined spectroscopic binary from Traven2020; companion is already recognized. |
| 403150111 | 3.911343 | UCAC4 156-214940 / Gaia6411065522024603392+6411065526320580736 | Existing SIMBAD double-star recognition and two Gaia sources about1.6arcsec apart. One component has strong RV variability. A previously unclassified eclipsing subsystem remains possible; a new third star is not established. |

## Catalog scope and interpretation

The checks include 120arcsec cones in current VizieR VSX, ASAS-SN II/366/catv2021 (a2021 catalog version archived2024), Shi et al.2022 EA binaries, Tang et al.2025 eccentric binaries, Battley et al.2022 **YOUNGSTER** EA/EB catalog, Gaia DR3 EB/variability classification and SIMBAD. The earlier internal filenames `cdips2022.tsv` are preserved unchanged; J/MNRAS/511/4285/table2 is YOUNGSTER, not CDIPS. Exact bytes, URLs, parameters, timestamps and SHA256 hashes are retained alongside every response. One initial Gaia HTTP503 for TIC281581467 was preserved and a successful immutable retry reused.

Additional checks for TIC38825979,197764088,52285789,403150111 include SB9, WDS, GALAH DR3, GALAH double-lined binaries, all four Gaia DR3 NSS tables and finite arXiv alias metadata searches. Zero results from arXiv metadata searches do not exclude mentions in full text or supplementary tables. NSS=0 does not imply a single star: TIC38825979 has Gaia RV amplitude174.9km/s while its binary is positively recognized in GALAH.

The complete GALAH Traven2020 table has12,760 successfully analyzed entries. The complete **TESS Ten Thousand** publication comprises **7,936 new EBs in Table3 and2,065 previously known EBs in Table4**, verified using Astropy's CDS/MRT parser. Both were checked by exact TIC and <=120arcsec separation; they must not be confused with Table2, which contains about900,000 unvetted/unvalidated neural-network detections. Position matches alone do not establish that an aperture's eclipses originate on the nearby source.

TIC403150111 has two Gaia sources about1.6arcsec apart. Both have large IPD multi-peak fractions (47%,33%), and one shows RV amplitude119.7km/s with a zero rounded chi-square p-value while the other has a5.37km/s amplitude and p=.982. This favors investigating an unresolved eclipsing subsystem, but Gaia summary RV scatter is not a solved orbital model, masses, or proof of a third star. Existing visual-pair recognition must remain explicit.

## Retained hidden-source leads

The separate full-table audit includes the six earlier stellar leads and proposed sources TIC273879360,260352037,737466840,260352041. They have no exact or120arcsec-position match in the complete Traven2020 successfully analyzed SB2 table or either vetted TESS Ten Thousand table. This is a bounded negative recognition check. It does not establish novelty, and any presence in unvetted Table2 must be recorded separately as a prior pipeline flag. Source attribution, stellar classification and further literature checks remain necessary.

## Reproduction and evidence

Run `/workspace/openrouter-check-venv/bin/python scripts/stellar_scout_novelty.py --download --lead-checks --full-crossmatch --report` from the repository. Existing successful snapshots are hash-verified and reused. [Searchable assessments](../../results/stellar_systems/scout_novelty.csv), [full target/cone rows and provenance](../../results/stellar_systems/scout_novelty.json), and [complete-table matches](../../results/stellar_systems/scout_fullcatalog_matches.json) retain the scientific scope. The two TESS machine-readable supplements and complete Traven table are preserved under `data/stellar_systems/scout_novelty/`, with immutable request manifests. No authenticated or paid API was used.
