# Catalog access and initial pilot context

Retrieval date: 2026-10-13 UTC. All tables below were downloaded, parsed as CSV, and hashed with SHA-256; HTTP connectivity alone was not treated as success. Individual manifests record the exact retrieval time, query, URL, schema, row count, and checksum. These are snapshots of public reports, not complete proof of novelty or validation.

| Source | Cached file | Parsed rows | Scientific scope |
|---|---|---:|---|
| NASA Exoplanet Archive `pscomppars`, `tran_flag=1` | `data/catalogs/confirmed_transiting.csv` | 4,811 | Confirmed reported transiting planet ephemerides and stellar properties |
| NASA Exoplanet Archive `pscomppars`, all planets | `data/catalogs/confirmed_all.csv` | 6,445 | Includes known nontransiting/RV planets that must not be called new if a transit is detected |
| NASA Exoplanet Archive `toi` | `data/catalogs/toi.csv` | 8,148 | TOI reports of all dispositions, not just confirmed planets |
| ExoFOP TESS community submissions | `data/catalogs/ctoi.csv` | 5,138 | CTOIs, including unpromoted reports absent from the TOI table |

Public endpoints used: [NASA TAP](https://exoplanetarchive.ipac.caltech.edu/TAP/sync) and [ExoFOP CTOI CSV](https://exofop.ipac.caltech.edu/tess/download_ctoi.php?output=csv). No credentials or paid APIs were used. A first NASA query contained an unsupported `pl_refname` column and received an explicit TAP SQL error; the schema was checked using a real table response, and that column was removed. Scientific reference columns (`pl_orbper_reflink`, `pl_tranmid_reflink`, `pl_trandur_reflink`) were retained.

## Reproduction and preservation

Run `python data/catalogs/retrieve_catalogs.py`, then `python data/catalogs/select_pilot.py` and `python data/catalogs/write_known_planets.py`. Existing cache files are reused only after manifest checksum verification. Raw table snapshots are retained locally. Each small pilot derivative also has a manifest identifying parent files. Refreshing a catalog should preserve the old snapshot and its manifest under a dated name first; public tables can change between sessions.

`data/known_planets.csv` has the requested columns `tic,name,period,epoch_bjd,duration_hours,source`. Epochs are catalog BJD values without an implicit offset; TESS FITS TIME usually requires the explicit BJDREF offset, commonly 2457000. Confirmed ephemerides are preferred to TOI ephemerides for matching periods within 0.1%; unmatched TOIs are also included. This matching tolerance only removes catalog duplicate reports, and is not a transit-identification or novelty criterion. The original ephemeris uncertainties and citations remain in the raw catalog.

## Selected pilot

Eight hosts were selected as manageable known small-planet systems: two confirmed transiting planets each, TESS magnitude below 12, stellar temperature below 5000 K, and radius below 0.9 solar. This is a purposeful feasibility sample rather than an unbiased occurrence-rate sample. Availability and cadence must still be checked at MAST. Multi-sector coverage will determine whether searching for additional weak or longer-period signals is practical. These hosts have published prior analysis; they are not claimed to be unstudied.

| Host | TIC | Known confirmed periods (days) | Additional current report |
|---|---:|---|---|
| TOI-776 | 306996324 | 8.246620, 15.665323 | CTOIs at these same known periods |
| TOI-1064 | 79748331 | 6.443868, 12.226574 | No extra period in queried TOI/CTOI snapshots |
| TOI-1468 | 243185500 | 1.8805201, 15.532477 | No extra period in queried TOI/CTOI snapshots |
| TOI-1470 | 284441182 | 2.527093, 18.088160 | The second confirmed signal lacks a separate TOI row in this snapshot; confirmed catalog takes priority |
| TOI-2095 | 235678745 | 17.664840, 28.172320 | **TOI-2095.03 at 57.6869019 days is already publicly reported** |
| TOI-4342 | 354944123 | 5.5382592, 10.688662 | No extra period in queried TOI/CTOI snapshots |
| TOI-1803 | 144401492 | 6.293287, 12.885779 | No extra period in queried TOI/CTOI snapshots |
| TOI-904 | 261257684 | 10.877200, 83.999700 | CTOI at the same known outer 84-day period |

The L 98-59 positive control (TIC 307210830) has known transiting b/c/d at 2.253114, 3.6906764, and 7.450729 days. Known RV planets e and f at 12.8278 and 23.064 days also occur in the all-confirmed snapshot. A signal at either of these RV periods would require interpretation as a potentially new transit of a known planet, not a previously unknown planet.

The negative statements above mean only “no additional matching TIC row in the catalog snapshots”; they do not establish that a system has no prior residual signal or literature report. Before promoting any candidate, search aliases and coordinates, compare periods/harmonics and ephemerides, inspect published supplemental data, and investigate false positives.

## Public literature recognition check

A public [arXiv API](https://export.arxiv.org/api/query) query of the eight host labels plus `L 98-59` returned **37 papers** on 2026-10-13. The exact query, retrieval time, and checksum are saved alongside `data/catalogs/pilot_arxiv.xml`; parsed titles, dates, URLs, and abstracts are in `pilot_arxiv.csv`. This is an abstract-level host-label search, not a systematic review: aliases, journal papers without an arXiv deposit, supplemental data, conference reports, and recently unindexed papers remain unchecked. NASA period references are also available in the original catalog rows.

| Host | Published context and implications | Accessible source |
|---|---|---|
| TOI-776 / LP 961-53 | Two known planets originally found in TESS sector 10; a separate 34-day RV signal was attributed to stellar rotation using long photometric monitoring and spectral indicators. Atmospheric papers now characterize the known planets. A roughly 34-day residual requires strong activity checks. | [Discovery abstract](https://arxiv.org/abs/2009.08338) |
| TOI-1064 | The 6.44/12.23-day pair was studied using CHEOPS, TESS, ground photometry, HARPS spectra, and speckle imaging. Catalog period source: Wilson et al. 2022. | [Wilson et al.](https://arxiv.org/abs/2201.03570) |
| TOI-1468 / LSPM J0106+1913 | Two planets at 1.88/15.53 days were confirmed with CARMENES and MAROON-X plus ground photometry. Known inner-planet emission has since been studied by JWST. | [Discovery abstract](https://arxiv.org/abs/2208.10351), [2025 emission study](https://arxiv.org/abs/2503.19772) |
| TOI-1470 | Published analysis itself recovered an initially unannounced 18.088-day second transiting planet; the 2.527-day planet was the initial TOI. The star rotates at 29±3 days. The lack of a second TOI row does **not** make its transit novel. | [González-Álvarez et al. 2023](https://arxiv.org/abs/2306.08516) |
| TOI-2095 / LSPM J1902+7525 | Two independent 2023 papers studied the 17.7/28.2-day pair with up to 24 TESS sectors and follow-up. Current **TOI-2095.03 at 57.6869019 days**, disposition PC, was created 2026-08-27 and is already reported, regardless of absence from those 2023 abstracts. | [CARMENES validation](https://arxiv.org/abs/2304.09220), [Independent TESS study](https://arxiv.org/abs/2304.09189), current TOI snapshot |
| TOI-4342 | Initially found through QLP improvements, including a multi-planet search. A 2026 RV study reports an **additional likely nontransiting 47.5-day Neptune-mass candidate**, minimum mass 17.8±3.0 Earth masses. This candidate is absent from the confirmed and TOI rows queried here, illustrating the failure of catalog-only novelty checks. | [Original QLP validation](https://arxiv.org/abs/2301.01370), [Parc et al. 2026](https://arxiv.org/abs/2601.22115) |
| TOI-1803 | Known 6.3/12.9-day pair near 2:1 resonance has significant measured TTVs: roughly 10 minutes for b and **40 minutes for c**. Fixed-epoch masks and folded fits must accommodate this before interpreting residual events. | [Zingales et al.](https://arxiv.org/abs/2412.05423) |
| TOI-904 | Both the 10.88-day and 84-day planets are published and statistically validated using TESS, spectroscopy, imaging, and TRICERATOPS. CTOI notes explicitly report three long-period events and an 84-day solution. | [Harris et al. 2023](https://arxiv.org/abs/2310.15118), current CTOI snapshot |
| L 98-59 control | In addition to the three transiting planets and confirmed 12.83/23.06-day RV planets, a 2025 paper reports a **1.7361 (+0.0007/−0.0008)-day RV candidate**, minimum mass 0.58±0.12 Earth masses and model-dependent 2.9–4.2σ significance. This extra candidate is not in the confirmed catalog. Another 2025 analysis reports transit timing variations below 3 minutes for the c/d pair. | [Five-planet architecture](https://arxiv.org/abs/2507.09343), [Additional RV candidate](https://arxiv.org/abs/2507.06413) |

Crossref searches using only the hyphenated host labels initially returned unrelated `toi` documents; those results were rejected as scientific evidence. A corrected title query with `TOI 4342 planet` returned the relevant 2026 and discovery publications, but broad Crossref querying also encountered HTTP 429 rate limiting. The NASA ADS abstract page returned HTTP 405. The successful arXiv feed provided independently accessible publication metadata and abstracts; none of the blocked or irrelevant responses were counted as a completed novelty check.

Prior nontransiting candidate periods and rotation periods are saved in `data/catalogs/pilot_prior_signals.csv`. These entries often lack a transit ephemeris, so they should be checked against residual periods rather than silently added to a transit mask. No previously unreported candidate has been established by this catalog/literature stage.

## Independent field targets: catalog and source checks

The original field pilot was expanded by the selection workflow from six to ten targets while these checks were performed. The following checks use the ten targets in the saved `data/field_sample.json`, independent of residual BLS results. Run `python data/catalogs/check_field_catalogs.py` and `python data/catalogs/check_field_external.py` to reuse the current caches and reproduce the checks.

- **Confirmed, TOI, and CTOI:** zero exact TIC or sky-position matches within 30 arcseconds for all ten targets. `field_nearest_catalog_entries.csv` retains the nearest record from each table, so this conclusion can be audited. This radius covers a little more than one TESS pixel but is not a complete contamination search.
- **Potential neighboring source:** TIC 348839628 is 94.2 arcseconds from CTOI 348842168.01, a publicly reported 13.92-day event with epoch 2459045.44 BJD, six-hour duration, Tmag 11.493, and reported 34-Earth-radius size. The CTOI notes cite “Taaki et al. in prep 2022” and its user disposition is PC. Its large reported size does not establish an EB but makes a planetary interpretation questionable. A residual near this period should be checked for contamination. The separation is roughly 4.5 TESS pixels, so an explicit aperture/pixel check is needed to determine whether any flux reaches the target aperture.
- **TESS EB site:** exact TIC searches returned three objects, but inspection of their detail pages explicitly shows **`In catalog? False`**, with blank period/epoch. These are light-curve-only records, not positive EB classifications. Seven targets returned zero objects. The initial interpretation of search hits as known EBs was corrected after detail-page inspection; no EB exclusion is warranted on these results alone. Both query and detail HTML are preserved, and the script requires explicit catalog inclusion instead of treating an object page as an EB.
- **SIMBAD:** primary object types alone were mainly proper-motion stars. Querying all object types reveals published activity classifications for some targets, as detailed below. These are cautions for detrending and false-positive analysis, not proof that transits cannot occur. SIMBAD magnitude bands are not interchangeable with the TESS band.
- **arXiv:** a combined search for the ten numeric TIC identifiers returned zero entries in searchable metadata. This does **not** search all full text, aliases, coordinates, journal-only publications, or inaccessible supporting tables. The exact query and null-result XML are preserved. No claim of novelty follows.

| TIC | SIMBAD main identity | Additional relevant SIMBAD information |
|---:|---|---|
| 348839788 | PM J06353-5737 | Double-star and eruptive/X-ray classifications. PM J06353-5737B is only 1.84″ away, Gaia G=14.72 vs target G=10.29. The nearby companion must be considered for dilution and signal localization. |
| 150359500 | EUVE J0625-60.0 | Double-star, eruptive, rotating, variable, UV/X-ray classifications. Neighbor UCAC4 150-007762 is 9.98″ away, G=15.90 vs target G=11.03. |
| 278727563 | TYC 8543-780-1 | Proper-motion/NIR star; no variable or EB classification in retrieved SIMBAD all-types row. |
| 150360085 | UCAC4 145-007490 | Eruptive-star classification in SIMBAD. |
| 348839628 | SCR J0635-5717 | Proper-motion/NIR star; nearby CTOI contamination caution above. |
| 260476305 | UCAC4 157-007462 | Proper-motion/MIR/NIR star. |
| 150299576 | UCAC4 143-007406 | Proper-motion/MIR/NIR star; catalogued faint galaxy 23.36″ away, G≈20.59. |
| 150359842 | UCAC4 147-007455 | MIR/NIR star. |
| 260415085 | UCAC4 157-007349 | Proper-motion/MIR/NIR star. |
| 260414764 | 2MASS J06220997-5922057 | Proper-motion/MIR/NIR star. |

These are finite catalog checks, not an assertion that any target is unstudied or free of contaminants. Full false-positive assessment should be reserved for independently robust event patterns; no credible planetary candidate has yet emerged from this catalog stage.

The Gaia DR3 TAP check returned rows for all ten supplied Gaia identifiers. Their positions are consistent with the TIC stars allowing for proper-motion epoch differences. TIC **348839788, 348839628, and 260415085** have `phot_variable_flag=VARIABLE` and general `SOLAR_LIKE` variability classifications. Their classifier scores are 0.307, 0.285, and 0.681 respectively; these are not calibrated planetary or astrophysical false-positive probabilities. None has a positive Gaia eclipsing-binary, rotation-modulation, or short-timescale table membership flag. The other seven have `NOT_AVAILABLE` and null variable-table flags, which means missing classification, not evidence of constancy. All ten have `non_single_star=0` and RUWE 0.88–1.16. These astrometric summaries do not rule out the SIMBAD companions or photometric blends.

The first Gaia request used incorrect boolean column names with a `_table` suffix and returned a parse error. The official TAP schema was retrieved, the query was corrected to the actual column names, and the resulting CSV was parsed. Request details, checksums, classifier rows, and the failure diagnosis remain in `data/catalogs/field_gaia_*`. This is a successfully completed finite Gaia catalog check, not a failed access test.
