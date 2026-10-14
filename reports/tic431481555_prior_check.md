# TIC 431481555: Gaia variability, stellar context, and finite prior checks

Checked 2026-10-14 UTC. **The target is already publicly classified as a
pulsating variable, while the approximately 5.635-day screened feature remains
unexplained.** Its broad/scattered photometric morphology requires an explicit
pulsation, beating, and detrending falsification test. The catalog classification
does not prove that these particular dips are pulsations or exclude an
additional planet. There is currently no credible new-planet claim.

## Positive variable-star recognition

Gaia DR3 source **4639812801472335104** has
`phot_variable_flag=VARIABLE` and published best class
**`DSCT|GDOR|SXPHE`**, score **0.5475003**. This composite label encompasses
pulsating-star alternatives rather than identifying a unique subtype. The
classifier score is a catalog classification metric, not a false-alarm
probability for the screened transit-like feature.

The VSX VizieR mirror, explicitly **Version 2026-10-05**, independently returns
an on-target entry named Gaia DR3 4639812801472335104, **OID2879383**, type
**DSCT|GDOR|SXPHE**, at RA49.42662, Dec−74.22826. It is 0.09407arcsec from
the supplied screening position and has the exact Gaia identifier association.
It lists G-band maximum/minimum magnitudes **12.900 and 12.920**, spectral
field A, and **no period**. The 0.020-mag range corresponds to roughly1.8%
full flux variation, but these sparse catalog extrema are not a measured
TESS variability amplitude or a model of the dips.

Gaia rotation-modulation, main-sequence-oscillator, and short-timescale
solution table queries each return no target rows. Thus no frequency, rotation
period, or pulsation solution from those tables is available here; their
absence does not negate the positive classifier/VSX recognition.

## Stellar and astrometric context

| Parameter | TIC v8.2 via current ExoFOP response | Gaia DR3 |
|---|---:|---:|
| Effective temperature K | 7147 ± 138.148 | GSP-Phot7310.3364, reported bounds7295.1904–7322.1167 |
| Radius solar | 1.42057 ± 0.0608542 | GSP-Phot1.6229, bounds1.5984–1.6499; FLAME1.5705849 |
| Mass solar | 1.59 ± 0.266237 | FLAME1.5292932 |
| logg | 4.33454 ± 0.083111 | GSP-Phot4.1368 |
| TESS magnitude | 12.6546 | — |
| Gaia G | ExoFOP12.9174 | 12.908891 |
| Distance/parallax | 905.517 ± 16.842pc | 0.985126 ± 0.018786mas; no distance inversion used |

The TIC row is explicitly dated2019-04-15; a current retrieval is not a new
spectroscopic measurement. The old QLP header gives Teff6750K, radius1.59396,
and TESS magnitude12.603. Radius differences matter for duration/size tests;
the Gaia radii are somewhat larger than TIC8.2 and broadly near the old header.

Gaia RUWE is **1.512188**, astrometric excess noise0.198046mas with
catalog significance**40.922966**. These merit caution about unresolved
astrometry/companions or model inadequacy but do not diagnose a binary.
The target has five-parameter astrometry (`astrometric_params_solved=31`),
IPD multiple-peak fraction0%, duplicated_sourcefalse, non_single_star0,
BP/RP excess factor1.1767365, and355 G observations with flux/error3226.2583.
The non_single_star0 value means no published Gaia NSS processing flag here,
not exclusion of any companion.

Gaia median radial velocity is **48.377563 ± 4.81554km/s** from18 transits.
No robust amplitude, chi-square probability, or renormalized RV goodness of
fit is supplied. The imprecise median is neither a mass measurement nor a
binary veto. No row is present in Gaia's eclipsing-binary solution table, and
the vari_summary field in_vari_eclipsing_binary is false; these do not prove
the signal is planetary.

## Photometric context and concrete falsification hypotheses

The existing sector-1 local vetting report gives period**5.635007716976d**,
epochBTJD1328.133361576761, and three qualified local events, with KSPSAP
local S/N9.95 and raw-SAP local S/N9.91. Its fitted exposure-integrated
trapezoid has duration **0.358279 ± 0.050932d**, or **8.60 ± 1.22hours**,
and depth0.0036636±0.0005902. The BLS box duration is4.80hours. The fit uses
a fixed period and approximate diagonal/covariance errors and has no limb
darkening/source validation; its formal errors are not robust posterior bounds.

Using TIC mass1.59, radius1.42057, period5.635008d and the local depth0.0033457,
the illustrative maximum central circular-transit duration including
`sqrt(depth)` is about **4.17hours**. With Gaia radius1.6229 and mass1.5293
it is roughly4.8hours. The 8.6-hour fitted width is therefore a physical
warning, not a definitive planet veto: eccentric geometry, uncertain stellar
parameters, and the difference between a flexible trapezoid width and a
physical transit duration must be considered. A planet claim cannot treat a
broad variability minimum as a well-measured contact duration. The TIC
depth-size proxy is about8.96 Earth radii, conditional on an undiluted event
on this target and appropriate nongrazing geometry.

A useful independent test is suggested by the **other frozen coarse BLS peak**:
period0.80569075599d, box duration0.16d, rejected for excessive duty cycle
(approximately19.9%). **Seven times this period is5.63983529d**, within about
0.0856% of the refined long-period feature. This relation may reflect
variable-star harmonics, sampling/selection, or beating; it is not itself
proof of a shared physical period. Period matching in the prepared exoplanet
catalog uses only2–5 integer factors, so this internal photometric7:1 relation
is recorded explicitly rather than silently altering the frozen catalog flags.

Prioritized tests are an independent Lomb–Scargle/frequency analysis of raw
and alternative detrended flux, conservative pulsation modeling/prewhitening,
event timing/morphology by cycle, and independent-sector recurrence. A
prewhitening choice must avoid subtracting possible transit events. Injection
recovery around the inferred pulsation structure and correlated-noise/search
trials tests would establish how often such broad features survive the coarse
screen in pulsating stars.

## Nearby Gaia/TIC sources

The current Gaia120arcsec cone returns14 sources including the target. There
is no listed neighbor within30arcsec. The nearest catalog neighbor is Gaia
4639812423514125184 = TIC650185495, at43.9478arcsec, G19.8884, full catalog
G-flux ratio0.001615 relative to the target. ExoFOP gives its TESS magnitude
18.9738 versus target12.6546, full TESS-flux ratio about0.00297. Without
aperture coupling and different passbands, that number must not be substituted
for measured contamination.

A brighter neighbor **Gaia4639812389155476864**, at**113.211arcsec**, has
G14.72584 and full G-flux ratio**0.18760**. Its RUWE0.95675, non_single_star0,
and missing variability flag provide no eclipse veto. It is absent from the
shorter ExoFOP nearby list, illustrating why the independent120arcsec cone
is needed. Bright wings/neighbor variability and background subtraction must
still be assessed from pixels; the full distant-neighbor flux does not all
enter the target aperture.

TIC's contamination ratio is0.000221, a predicted quantity for an assumed
aperture, not a measured source-localization test. The derived Gaia neighbor
CSV includes every cone source, separations, exact identifiers, and explicitly
G-band full catalog flux ratios.

## Finite planet/EB/literature checks

The cached CatalogIndex with120arcsec radius has no exact-TIC, positional
host/nearby, period, or2–5 harmonic/subharmonic catalog matches for either
frozen coarse peak (0.80569076 and5.62460812d). The index includes NASA
confirmed planets of all discovery methods, TOIs, CTOIs, DIAmantev2,
Prsa2022 TESS EBs, Olmschenk2021, and Yu2019. Exact source snapshot
dates/hashes and per-period results are retained; computation date is not a
claim that every underlying frozen table is current.

The live ExoFOP target response has empty TOI, CTOI, confirmed planet,
planet-parameter, spectroscopy, imaging, time-series, and uploaded-file
lists. Live Villanova explicitly says “This TESS ID is not in the catalog.”
SIMBAD's120arcsec cone returns no objects, so no target-associated SIMBAD
aliases or bibliography were available. Finite arXiv metadata queries using
TIC, 2MASS, GaiaDR3, UCAC4 and APASS tokens return zero results. These
negative finite checks do not override the positive Gaia/VSX variable-star
recognition and do not establish that the5.635-day feature is a new planet.

Outstanding novelty checks include full TCE/DV listings, recent journal
supplements, other variable-star surveys, Planet Hunters, all alias variants,
and non-arXiv full text. The current classification is **unexplained feature
on a previously recognized pulsating variable, requiring variability and
independent-data falsification**.

## Evidence preservation and reproduction

Public responses, exact URLs/ADQL, actual UTC retrieval times, HTTP status,
sizes, and SHA-256 sidecars are under
`data/broad_catalogs/tic431481555_*`. New suffixes are
`cached_matches_current.json`, `exofop_current.json`,
`exofop_nearby_current.csv`, `gaia_target_current.csv`,
`gaia_stellar_current.csv`, `gaia_class_current.csv`,
`gaia_eb_current.csv`, `gaia_neighbors_current.csv`,
`gaia_neighbors_derived_current.csv`, `gaia_rotation_current.csv`,
`gaia_ms_oscillator_current.csv`, `gaia_short_timescale_current.csv`,
`tesseb_live.html`, `vsx_current.tsv`, `simbad_current.csv`,
and `arxiv_current.xml`. All16 saved source/derived body hashes were verified
against their manifests. Existing response snapshots, frozen rankings,
broad source CSVs, and scientific pipeline code were not overwritten.

For the Gaia/basic stellar/cached-index lookups, invoke
`handle(431481555)` from `scripts.secondary_leads_catalog_check.py` in the
repository Python environment; it reads the sector-1 survey SQLite database
in mode=ro and timestamps any existing response name. Additional queries can
be replayed using exact public URLs/ADQL in each manifest with TLS verification
enabled and fresh output names. Photometric quantities refer to
`candidates/broad_s1/tic431481555_p5.624608/vetting.json`; this report does
not rewrite its frozen measurements or flags. No paid APIs, credentials, API
keys, or private observer records were accessed.

