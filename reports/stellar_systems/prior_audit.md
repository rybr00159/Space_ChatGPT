# Recognition audit of six retained eclipsing-source leads

Checked 2026-10-14 UTC. **No prior eclipsing-binary identification was found for these six sources in the expanded catalogs and finite literature checks below. This is not proof of novelty, and no newly uncatalogued star or companion mass is established.** These are known catalogued stars with evidence for potentially unreported eclipsing stellar companions or blended systems.

## Most useful leads

| Priority | TIC | Period assessed, days | Stellar-companion evidence | Main limitation |
|---:|---:|---:|---|---|
| 1 | 370247818 | 3.522194946 | Eight sectors; raw secondary 1,327 ±232 ppm at phase0.47546; independent primary difference images and conditional geometry/light-budget checks favor a luminous companion | No mass, joint orbit fit or resolved close-source exclusion; stellar radii/grazing geometry uncertain |
| 2 | 273879358 aperture → background TIC 273879360 | 7.586891287 | Four independent native pixel sectors now place deep primary and secondary eclipses on the background source, 65.59 arcsec away | Foreground Gaia RV variability belongs to a different star; see the source-specific supplement and dossier |
| 3 | 300814564 | 1.860044776 | Twenty-three usable sectors and a raw secondary about2,180 ±62 ppm | Close Gaia source1.855 arcsec away prevents assigning the eclipsing companion to the brighter component without resolving them |
| 4 | 238202207 | 4.643651162 | Seventeen sectors; alternating depths of the same sign in16/17; conservative parity contrast4.39 conditional standard errors; Gaia RUWE4.13/IPD multiple peaks22% | Target/aperture localization and unresolved-source structure prevent component assignment |
| 5 | 269994642 | 3.492744247 | Recurring approximately2.7% eclipses in2018/2025 independent pixels; doubled-period binary fits; parity3.33 conditional standard errors; Gaia RUWE2.08 | Only two independent sectors for parity, and no companion mass |
| 6 | 270312735 | 8.370105786 | Six sectors; raw off-half-phase secondary about5,584 ±619 ppm | Strong catalog contamination and a bright neighboring star require source localization |

Periods are conditional measurements from the earlier retained analyses, not new independent period fits in this audit. For238202207 and269994642 the table adopts the favored doubled binary interpretation. Different choices of parity label swap the deep/shallow labels without changing the evidence. The statistics above are targeted same-photon diagnostics, not search-wide significance or independent probabilities.

## Expanded recognition checks

New public response bytes, exact URLs/query parameters, timestamps, HTTP status and SHA-256 hashes are retained in `data/stellar_systems/prior_audit/`. The reproducible table is [prior_audit_results.csv](prior_audit_results.csv), with [source inventory](prior_audit_sources.json).

- **Actual Gaia DR3 published solutions**, not merely the `non_single_star` flag: `nss_two_body_orbit`, `nss_acceleration_astro`, `nss_non_linear_spectro`, `nss_vim_fl`, `vari_eclipsing_binary`, `vari_classifier_result` and `vari_summary` were queried by exact source ID. All returned zero rows for all six sources and three important nearby sources. No solution does not exclude a binary. A metadata discovery query failed with HTTP500; its body is retained and was not treated as a membership check. The seven actual science-table queries succeeded with HTTP200 and valid CSV headers.
- **Live ASAS-SN Variable Stars Database**: six120-arcsec cones with nonperiodic and unclassified variables included, and discovery-only filtering disabled, each explicitly returned zero light curves. This adds a contemporary lookup to the **VizieR II/366/catv2021** mirror, which reports2021 catalog version and CDS retrieval from the live database on2024-09-10. The mirror also returned zero rows. These checks do not claim complete g-band/current-survey detection coverage.
- **Other TESS binary catalogs**, beyond Prša2022: Shi/Qian/Li2022 EA systems (`J/ApJS/259/50/table1`), Tang et al.2025 eccentric systems (`J/AJ/170/356/table1`, first CDS publication2026-09-14), and Battley et al.2022 YOUNGSTER EA/EB candidates (`J/MNRAS/511/4285/table2`). All six120-arcsec cones are empty. Independent TAP queries by exact TIC IDs also give zero rows, protecting against coordinate/indexing errors in the cone lookup.
- **SB9** (`B/sb9/main`), **Washington Double Star** (`B/wds/wds`) and **Tian et al.2020 Gaia DR2 wide-binary catalog** (`J/ApJS/246/4/catalog`): all six120-arcsec cones return zero rows. These are selected catalogs; the nulls do not establish singleness or exclude prior recognition elsewhere.
- **ATLAS2018 candidate variables** (`J/AJ/156/241/table4`): all cones empty. Only269994642 lies inside that survey's declination coverage. The other five lie south of−30°; their nulls have no novelty value. Other selection/incompleteness still applies within coverage.
- **Prior retained public sources**, independently byte-hash verified: current ExoFOP target/nearby pages, actual Prša membership table, VSX, Gaia neighborhoods, SIMBAD and finite arXiv alias results. The earlier source reports explicitly describe finite literature scope; absence of indexed aliases is not absence from full text or supplementary tables. TIC300814564's five SIMBAD references concern TESS-HERMES, APOGEE and GALAH data releases/spectral analysis, not an EB classification. APOGEE/GALAH binary value-added products are not exhaustively checked in this audit. ADS and all survey TCE/supplementary products remain unexhausted.

The narrower earlier analyses are [370247818](../tic370247818_prior_check.md), [269994642](../tic269994642_prior_check.md), and [238202207](../tic238202207_prior_check.md). Shared response bodies were preserved, and no frozen broad-search or prior physical-fit source was changed.

## Previously available multiplicity evidence

For273879358, Gaia's RV p-value and robust amplitude are **already published clues**, not fresh stellar detections from this project. They are measured variability diagnostics, not an orbital semiamplitude, companion mass, NSS orbit or confirmation that the same object produces the TESS signal. Establishing recurring primary/secondary eclipses on that source could still be a useful new classification/ephemeris result if broader literature checks survive.

For300814564, the Gaia pair has nearly equal parallaxes and similar proper motions. Exact DR3 positions give **1.855 arcsec** separation; the earlier informal3.45-arcsec number was erroneous and is not adopted. Both flux sources were already in Gaia/TIC. The potential result is an additional unresolved **inner eclipsing companion**, possibly within a wide system, not discovery of the resolved Gaia neighbor. WDS/Tian nonmatches are insufficient to claim that the wide association itself is newly recognized.

## Reproduction and interpretation

```bash
/workspace/openrouter-check-venv/bin/python scripts/stellar_prior_audit.py --download
```

Existing snapshots are hash-validated and reused; source bodies are never overwritten. Without `--download`, the script reconstructs the six-row machine-readable results and verifies retained sources. A failed or malformed catalog response is an error, never an empty membership result. This audit runs no paid API calls.

The recommended language is **“potentially unreported eclipsing stellar-companion candidate; no prior EB identification found in the listed checks.”** An accepted new star/system claim needs reliable component localization, classification and broader recognition review. The next useful analysis is spatial and physical vetting of370247818 and273879358; mass-sensitive spectroscopy is the strongest eventual discriminator.

Latest source attribution: the eclipses originally associated with TIC273879358 originate on or near background TIC273879360. Foreground RV variability is not RV confirmation of that background system. Read [the updated source dossier](../../candidates/stellar_systems/tic273879360/README.md). Complete vetted Kostov2025 Table3/4 and Traven2020 SB2 checks also yielded no matching earlier source among these retained fields; scope remains finite.
