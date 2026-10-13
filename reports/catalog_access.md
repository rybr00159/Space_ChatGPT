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
