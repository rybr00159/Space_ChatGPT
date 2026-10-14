# Finite source and prior checks for three strong warning leads

Checked on 2026-10-14 UTC. This investigation covers TIC 300814564,
273879358, and 270312735, selected for strong raw-flux dips together with
centroid/aperture or secondary-eclipse warnings. **None should currently be
promoted as a credible new planet candidate.** Source attribution and stellar
binary explanations remain unresolved. These checks have not changed the
frozen broad-screen database, ranking, model, or current-stellar table.

The three retained signals are unrecognized in the finite checked planet and
published-candidate snapshots. That is neither a discovery nor evidence of
planetary origin. Strong counterevidence described below takes priority over
catalog absence.

## Current stellar and Gaia source context

|TIC|Examined local period (d)|Current TIC radius (solar)|Gaia GSP-Phot radius (solar)|Gaia target RUWE|IPD multi-peak (%)|Current TIC contamination ratio|
|---|---:|---:|---:|---:|---:|---:|
|300814564|1.861107|1.60118 ±0.07070|1.5354 [1.5250, 1.5460]|1.02165|4|0.276102|
|273879358|7.581870|0.982455 ±0.04354|1.0234 [1.0173, 1.0280]|1.02821|0|0.016991|
|270312735|8.373627|1.32993 ±0.05879|1.3810 [1.2913, 1.4023]|0.92374|0|0.552935|

Current TIC values come from the public ExoFOP target response, whose TIC
parameter rows are labeled v8.2; they are current served catalog estimates,
not newly measured stellar parameters. Gaia ranges are supplied catalog
lower/upper bounds. Target Gaia identifiers are, respectively,
5273775087036017024, 4692911684016155008, and 4617279547612110208.
All three Gaia target records have `non_single_star=0`, empty queried
`vari_classifier_result` and `vari_eclipsing_binary` results, and
`phot_variable_flag=NOT_AVAILABLE`. Missing published solutions/classifications
are not evidence that a star is single or photometrically constant.

|TIC|Gaia robust RV amplitude (km/s)|RV transits|Constant-RV chi-square p|RV renormalized GOF|
|---|---:|---:|---:|---:|
|300814564|14.04164|20|0.37447|0.81344|
|273879358|17.51296|22|2.600762e-6|4.63856|
|270312735|12.49375|23|0.10174|0.88557|

**TIC 273879358 has independent significant Gaia RV variability in the
reported metric**, supporting investigation of a stellar companion together
with its recurring secondary dip. No Gaia orbital period was found in these
queries. A robust amplitude is not a fitted orbital semiamplitude; amplitudes
alone for the other two targets do not establish binarity. Gaia target RV
uncertainties are 1.08127, 1.34858, and 0.88954 km/s, respectively.

## Per-target source and photometric assessment

### TIC 300814564

The preserved sector-1 analysis finds a 1.861107-day dip, raw SAP local S/N
44.56, and 11 sampled primary events. Its approximately 13,629-ppm processed
local depth would imply roughly 20 Earth radii if hosted by the current TIC
star, before dilution or grazing corrections; this is a simple radius proxy,
not a companion-radius measurement. The fitted trapezoid is broad/grazing
(`ingress_fraction=0.386`) and does not establish a planetary shape. Half-phase
secondary local S/N is 8.21 in the candidate README; the separately selected
additional phase is near 0.50167 and requires trials-aware independent checking.
The QLP centroid baseline tests show x/y warning scores 7.18/14.32, explicitly
without pixel source attribution.

There is a **real Gaia-resolved neighbor at 1.855 arcsec**:
TIC 765053491 = Gaia DR3 5273775091332140800, G=13.76441 versus the target
G=11.91909, with G flux ratio 0.18276. ExoFOP lists T=13.281 versus target
T=11.5237. The neighbor has RUWE=1.37816, a six-parameter Gaia astrometric
solution (`astrometric_params_solved=95`), no supplied BP−RP or radial velocity,
and parallax 2.07150 ±0.01634 mas, consistent with target
2.06743 ±0.00986 mas. Proper motions are similar but differ measurably:
neighbor (−4.70972, +6.27115) versus target (−4.63682, +6.41417) mas/year.
This is compatible with a physically associated source but does not prove
binding, orbital geometry, or which source eclipses. A diluted stellar eclipse
in this close source is a viable alternative. TESS pixels alone may not resolve
1.86 arcsec; PSF/difference-image constraints need honest resolution limits.

The retained Gaia 120-arcsec cone has 69 rows, including the target. A further
G=13.57407 star lies at 44.044 arcsec. VSX Version 2026-10-05 returns two
distinct sources: Gaia 5273774988252920704 (ROT, 95.970 arcsec) and
5273775297491241344 (S, 108.916 arcsec), without catalog periods in the
returned rows. They are not recognition of the target's transit period.

SIMBAD associates the target with TYC 9180-825-1 (`otype=*`) and five
references, all broad spectroscopy/data releases: TESS-HERMES DR1 (2018),
APOGEE DR16 (2020), GALAH+ DR3 (2021), APOGEE Net (2022), and SDSS DR17
(2022). Their titles and abstracts are preserved. This finite reference set
contains no prior planet claim; reference-list scope is incomplete. Both the
target and close-neighbor ExoFOP responses have empty TOI/CTOI/planet/files
arrays, and their live Villanova exact-TIC pages say the TIC is not in the
catalog. Required next tests: independent sectors, physically plausible
two-eclipse modeling, and source localization/constraining the close pair.

### TIC 273879358

The existing independent multisector analysis reports recurring primaries
across sectors 1, 27, and 29, with a timing-fit period
7.586894227 ±0.000027025 d. The finite catalog search used the original local
7.581870-d period, with 1% tolerance; the refined period differs by only 0.066%
and does not change the cached-match conclusion. The existing targeted
additional-phase search finds a **2005 ±218 ppm secondary at phase 0.69368**,
raw SAP S/N=9.10. The usual half-phase test alone would miss this evidence.
Combined with Gaia's 17.51-km/s robust RV amplitude and p=2.6e-6 variability,
an eccentric stellar binary is a leading alternative. A phase search introduces
trials, and this conclusion still needs independently recurring secondary events
and physical modeling. These facts do not yet assign both eclipses to one star.

The Gaia 120-arcsec cone has 28 rows. Nearby TIC 273879357/Gaia
4692911684016266624 lies at 26.413 arcsec, G=18.45959; TIC 273879362/Gaia
4692911688311414144 lies at 28.532 arcsec, G=16.89394 with RUWE=1.84969.
Their G flux ratios relative to the target are 0.00380 and 0.01607. Full
eclipse/dilution feasibility depends on TESS colors and aperture coupling;
G ratios alone must not assign or rule out a contaminant. More distant
G=13.03209/13.90376 sources at 77.071/82.632 arcsec are also retained.

SIMBAD's 120-arcsec cone contains only a different AGN at about 62 arcsec,
Gaia 4692911722671104512; it is not the target's literature association.
No SIMBAD target row or reference list was obtained within the 2-arcsec
association threshold. Current VSX returns no rows. The live Villanova page
exists as a general data/source page, but explicitly states **“No associated
EB entries found.”** Page existence is not EB membership. The required next
investigation is binary/source validation, not promotion based on absence
from the planet list.

### TIC 270312735

The preserved multisector results find primary recurrence across six examined
sectors, with timing-fit period 8.370105786 ±0.000026062 d. The original
8.373627-d period used in the 1%-tolerance cached check differs by only 0.042%.
Existing targeted phase scanning finds a **5526 ±602 ppm secondary at phase
0.51262**, raw SAP S/N=9.02, comparable in depth to the roughly 6,000-ppm
primary. This is strong stellar-eclipse counterevidence requiring independent
phase-consistent confirmation and source attribution. The zero/negative
half-phase test at precisely 0.5 did not exclude this nearby-phase signal.

A brighter real Gaia source TIC 270312732 = Gaia DR3 4617279616331585408,
associated with SIMBAD TYC 9491-2697-1, lies at **38.049 arcsec**. It has
G=11.83056 versus target G=12.12600 (G flux ratio 1.31273), T=11.2998 versus
11.7036, and **RUWE=6.15916/IPD multi-peak fraction=58%**. Its astrometric
excess noise is 0.82303 mas (significance 984.59), so reliable single-source
astrometric interpretations are inappropriate. Its Gaia parallax is
3.01465 ±0.07370 mas versus target 2.34379 ±0.00982 mas; these Gaia estimates
do not support casually asserting a bound association. Neighbor robust RV
amplitude is 5.36849 km/s, but constant-RV p=0.93945 is not RV-variability
evidence. Its ExoFOP TOI/CTOI/planet/files arrays are empty and its live
Villanova page has no associated EB entries. High RUWE and brightness make
it a source-attribution priority, not proof that it produces the signal.

The Gaia cone has 32 rows. Other examples are G=19.48257 at 9.759 arcsec,
G=16.48402 at 41.178 arcsec, and G=13.63126 at 113.579 arcsec. SIMBAD has no
target association within 2 arcsec; its sole cone result is the brighter
neighbor, whose `nbref=0`. Current VSX returns no rows. The target's live
Villanova page states “No associated EB entries found,” not known membership.

## Finite recognition and reproducibility scope

All three cached checks use `CatalogIndex(position_radius_arcsec=120)` and
25,948 retained rows across NASA Exoplanet Archive confirmed planets
(all detection methods, including RV), NASA TOIs, ExoFOP CTOIs, DIAmante v2,
Prša et al. published TESS EB membership, and vetted candidates from
Olmschenk 2021 and Yu 2019. None yields a target or 120-arcsec positional
match; consequently no exact or 2–5x integer harmonic/subharmonic match is
present at 1% period tolerance. Each derived check preserves source catalog
checksums, dates, limitations, and query metadata. These are dated snapshots,
not an exhaustive current all-sector TCE or supplementary-literature index.

The current ExoFOP target responses have no confirmed-planet string and
empty TOI, CTOI, planet-parameter, imaging, spectroscopy, time-series, and
uploaded-file arrays. Finite arXiv metadata queries using all served ExoFOP
aliases and available SIMBAD aliases return zero results for each target.
This does not search all journal full text, figures, unpublished follow-up,
supplementary tables, or alternate unlinked designations. **Novelty remains
incompletely assessed; no catalog-absence claim of discovery is warranted.**

All 51 retained source/derived body hashes were verified against their own
manifests. The complete source inventory with URL, ADQL/alias query, UTC
retrieval time, byte size, and SHA-256 is
`reports/strong_warning_leads_20261014_sources.json`; raw responses are
`data/broad_catalogs/tic{ID}_*`. No secret or paid-model request was used.
TLS verification remained enabled. No existing scientific snapshot was
overwritten.

The reusable retrieval and neighbor-derivation script is
`scripts/strong_warning_catalog_check.py`. Run it with the repository's Python
environment from the checkout. Repeated retrievals retain timestamped new
snapshots rather than replacing older responses. It reads the cached catalog
index and writes only these three targets' prefixed response/derived files.

The two existing multisector README files repeat the sentence identifying
“CD-79 151 at about 86.6 arcsec.” This report does not endorse that copied
source statement: the three sky positions are distinct and the source cones
above are the correct target-specific evidence. The candidate-document owner
was notified to replace that statement without altering the frozen models.
