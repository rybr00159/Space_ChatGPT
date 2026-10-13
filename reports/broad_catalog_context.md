# Broad-survey catalog context

Checked on 2026-10-13. This work prepares fast catalog matching for the expanded
survey; it does not count as light-curve analysis and does not establish novelty.
The source NASA Exoplanet Archive and ExoFOP CSV snapshots are unchanged and their
recorded SHA-256 hashes were verified before use.

| Catalog | Cached rows | Scope and source |
|---|---:|---|
| NASA confirmed planets, `pscomppars` | 6,445 | All discovery methods; retrieved 2026-10-13, including 1,198 RV-only planets. |
| NASA TOI table | 8,148 | All dispositions; retrieved 2026-10-13. False positives are retained as prior reports. |
| ExoFOP community TOIs | 5,138 | All community submissions; retrieved 2026-10-13. |
| DIAmante v2 | 1,160 | Published candidate catalog accompanying Montalto 2022/2023, from MAST; retrieved 2026-10-13. DOI: [10.17909/t9-p7k6-4b32](https://doi.org/10.17909/t9-p7k6-4b32). |
| Published TESS EB catalog | 4,584 | Prša et al. 2022 ephemerides for 4,578 stars, sectors 1–26. [CDS J/ApJS/258/16](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJS/258/16). Retrieved 2026-10-13; not the full current live EB catalog. |
| Olmschenk et al. 2021 | 185 | Human-vetted candidate table; [CDS J/AJ/161/273](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/161/273). Row count refers to the CDS table, not newly discovered planets. |
| Yu et al. 2019 | 288 | Published TCE/planet-candidate table; [CDS J/AJ/158/25](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/158/25). Includes low-ranked events; table membership is a prior report, not validation. |

Combined index: 25,948 entries and 19,266 distinct supplied TIC identifiers.
There are 307 confirmed-planet records without TIC identifiers; these must be
checked by sky position or aliases. Position matching covers their coordinates
when the analyzed target has usable RA/Dec.

## Reproducible interface

```python
from broad_crossmatch import CatalogIndex

catalogs = CatalogIndex()
flags = catalogs.match_signal(tic, period_days, ra=ra, dec=dec)
keep_for_residual_search = catalogs.transit_mask(tic, times_btjd, ra=ra, dec=dec)
```

Run from the repository with `src/` on the Python path. Or annotate an existing
CSV without modifying it:

```sh
python src/broad_crossmatch.py input.csv --output output_catalog_flags.csv \
  --metadata data/broad_catalogs/catalog_index.manifest.json
```

Input columns are `tic`, optional `ra`/`dec`, and `period_days` or `period`.
Output contains scalar catalog flags, period flags, and JSON match evidence.
Initialize once per process and reuse the index. This uses a TIC dictionary and
a unit-sphere KD-tree rather than one remote catalog request per star.

Restore supporting files and derive missing normalized tables with:

```sh
python src/broad_crossmatch.py --retrieve-support
```

Raw sources, retrieval URLs, checksums, and source metadata are under
`data/broad_catalogs/`. Existing data are preserved; checksum mismatch raises an
error. A remote file that changes after an exact snapshot was recorded is not
silently accepted as that snapshot. The original NASA/ExoFOP tables can be
retrieved with the existing `data/catalogs/retrieve_catalogs.py` script.

## Classification rules and limits

- Exact TIC or a position match within 2 arcsec counts as a provisional host
  association. Other matches within 30 arcsec flag possible contamination.
  The two categories are separate in the evidence. Catalog coordinates have
  not been propagated using proper motion; ambiguous close blends still need
  pixel/astrometry checks.
- Periods within 1% of a reported period are flagged. Integer harmonics and
  subharmonics from 2 to 5 are also flagged. These coarse screening flags do not
  establish that two signals are identical; ephemeris alignment and raw data
  remain necessary.
- All confirmed planets are included, including nontransiting RV planets
  without transit epochs. A new transit at a known RV period is not a new
  planet discovery merely because it was absent from a transiting-only list.
- Known-host status alone is not a reason to discard a star. A different
  recurring period can motivate a search for additional planets.
- Residual masks apply only to associated host entries with a period and
  transit epoch. Missing epochs are never invented. Default width is 1.7 times
  the listed duration, with an explicit 3-hour duration fallback if only the
  duration is absent. Suspect masks covering >=30% of a period are skipped.
  RV-only entries and nearby contaminant ephemerides are not masked. Stale
  ephemerides, TTVs, and long baselines can still make simple masking incomplete.
- DIAmante uses linear ephemerides when complete, otherwise its fitted transit
  epoch and BLS period, all in BTJD. Other transit epochs are converted from
  BJD to BTJD only according to the source schema. All durations are in days.
- Catalog nonmatch produces `no_cached_period_match` and
  `unchecked_literature_and_all_aliases`, never a discovery label.

The live Villanova website provides light-curve-only entries as well as EBs.
An object-link hit is **not** EB membership. Multiple requests to its advertised
bulk export/search with explicit membership filtering returned HTTP 500; the
published CDS membership table is used instead. Later-sector and later-reported
EBs remain incompletely covered. Promising leads still require per-target live
checks with explicit `In catalog=True`, not a search-page occurrence.

The attempted CDS path for Eisner et al. 2021 (`J/MNRAS/501/4669`) returned
`Table or Catalog not found`; that response is retained and marked unusable in
its manifest. No entries from it are used. The initially attempted NASA TCE
documentation URL returned HTTP 404 and the NASA TAP schema exposed no TESS
TCE table through the checked table-name query. These checks **do not** establish
that prior TCE information is unavailable. A full all-sector SPOC/QLP TCE/DV
crossmatch, recent publications, Planet Hunters reports, additional literature
tables, all aliases, and ExoFOP notes remain follow-up work for surviving leads.

## Validation performed

Smoke checks exercised confirmed transit recognition, period harmonics,
masking of an actual epoch, recognition of a known RV-only planet with no
invented mask, exact positional association, and a nearby-period contaminant
that must not inherit the target's known-period label. All passed.

The 775 saved pilot peaks were annotated to a temporary file as an integration
check; the source table was untouched. The result classified 19 peaks as direct
known-period matches and 50 as harmonics, with 706 having no cached period
match. Those 706 are **not new candidates**: the earlier light-curve vetting
already rejected the strongest unexplained pilot signals. No prior-report host
association was added for the ten field-pilot stars from the extra publication
or published EB tables. This is a catalog scope result, not an assertion that
the stars have never been studied.
