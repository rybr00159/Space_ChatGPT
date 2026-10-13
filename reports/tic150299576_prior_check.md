# TIC 150299576: focused prior recognition check

Checked on 2026-10-13 UTC in response to an **unpromoted weak two-event lead near 11.27–11.29 days**. This document covers public prior reports and nearby catalog sources. It does not establish a planet, a coherent orbital period, or novelty. The light-curve investigation and held-out-sector tests determine whether the lead warrants further work.

The target is UCAC4 143-007406 / 2MASS J06202884-6130443 / WISEA J062028.88-613045.4, Gaia DR2/DR3 5481468435475904256. SIMBAD also lists WISE J062028.87-613045.5, GALAH 151223001301180, AP J06202884-6130443, and TIC 150299576. The identifier list is saved in `data/catalogs/tic150299576_simbad_identifiers.csv`.

Finite checks completed:

| Source | Result | Scope and limitations |
|---|---|---|
| Current saved NASA all-confirmed, TOI, and ExoFOP CTOI snapshots | Zero exact TIC matches or matches within 30″ | Catalog absence alone does not prove novelty; not all unpublished or literature-only signals are catalogued. |
| Live [ExoFOP target](https://exofop.ipac.caltech.edu/tess/target.php?id=150299576) and JSON export | Confirmed planets blank; TOIs 0; CTOIs 0; planet parameters 0; uploaded follow-up files and time-series observations 0 | A populated stellar target page does not itself imply prior transit analysis or exclude reports outside ExoFOP. HTML and full JSON with URL/time/checksum manifests retained. |
| arXiv API, all nine SIMBAD aliases | Zero matching metadata entries | This searches metadata, not every full-text table, all publication services, or unreported events. Exact request and null-result XML retained. |
| SIMBAD linked references | Eight references retrieved, titles and abstracts reviewed | References concern TESS-HERMES spectroscopy, TESS habitable-zone target lists, Gaia nearby stars, GALAH/APOGEE/SDSS catalogs, and observing-capability simulations. No retrieved abstract reports an 11.3-day planetary signal. Their complete target tables were not systematically inspected for every candidate report. |
| Live TESS EB site, exact TIC search | Zero objects returned | Does not exclude an uncatalogued eclipsing binary or contaminating source. |
| Gaia DR3 | RUWE=1.132; `non_single_star=0`; `phot_variable_flag=NOT_AVAILABLE` | Missing variable classification does not establish photometric constancy or exclude blends. |

No prior planetary report near 11.3 days was found **within this checked scope**. The correct novelty status remains **not established**, and signal validity remains independent of these catalog results.

ExoFOP lists TIC v8.2 stellar values Teff=3953±157 K, radius=0.616875±0.018256 solar, mass=0.601643±0.0202403 solar, distance=70.4737±0.106 pc, and Tmag=11.2229. These are catalog estimates, not a new stellar characterization.

## Contamination requiring falsification

The live nearby-TIC export lists nine sources within one arcminute. The closest is **TIC 150299575, 7.83″ away, Tmag=16.8134**. It is fainter than the target by 5.5905 TESS magnitudes, corresponding to an approximate TESS-band flux ratio of 0.0058. If comparable fractions of target and contaminant flux entered the aperture, an eclipse of roughly **10–20% on this source** could create a diluted 600–1140 ppm event. This is a possible alternative explanation, not evidence that the source actually eclipses. Aperture throughput, crowding, and the actual pixel location of a coherent signal must be tested before any planetary interpretation.

A second brighter neighbor TIC 150299571 lies 51.68″ away at Tmag=15.398; it requires an aperture/point-spread-function assessment if it enters the stamp or aperture. The catalog contamination estimate 0.011575 does not eliminate these alternatives. The raw nearby-source table is `data/catalogs/tic150299576_nearby_tic.csv`; SIMBAD also places a faint G≈20.59 galaxy 23.36″ away.

Reproduction: `python data/catalogs/check_tic150299576.py` retrieves/reuses the ExoFOP exports and alias query. `python data/catalogs/check_field_catalogs.py` reproduces saved-catalog matches. Linked-reference query URLs are recorded in the `tic150299576_simbad_references.csv.manifest.json` file. Scientific snapshots are retained rather than discarded when queries change.
