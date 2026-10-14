# TIC 238202207: recognition scope, stellar validity, and source warnings

Checked 2026-10-14 UTC. **A recurring transit-like feature remains
unexplained, with significant Gaia source-quality and stellar-model
warnings.** No prior planet/EB/variable identification was found in the finite
checks below, but this is not proof of novelty or planetary origin. The
short, approximately1.4%-deep event must be evaluated against blended and
grazing stellar-eclipse alternatives.

## Principal findings

- Current Gaia DR3 source **4629516356114834944** has **RUWE4.1278853**,
  **IPD multiple-peak fraction22%**, astrometric excess noise0.574783mas
  and significance187.23259. These provide substantial unresolved-source/
  astrometric-fit warnings. They do not alone diagnose an EB or identify the
  transiting source.
- Stellar radii disagree materially: **TIC1.62747 solar radii**, **Gaia
  GSP-Phot1.1929**, and **FLAME1.63394**. The unresolved-source warnings
  make precise single-star radius interpretation particularly uncertain.
- Gaia median RV is−40.869778±5.347638km/s over21transits, with **no robust
  amplitude or RV variability statistics**. It cannot supply a planet mass
  or exclude a stellar companion.
- A120arcsec Gaia cone contains two nearly equally bright distant stars:
  **TIC238202215 at96.66arcsec, T13.5422**, and **TIC238225954
  at115.37arcsec, T13.4381**. The target is T13.4232. The latter has
  Gaia RUWE10.8034. Aperture/difference-image tests must determine whether
  either contaminates the actual light curve.
- Current ExoFOP has empty target TOI, CTOI, confirmed-planet, prior
  spectroscopy/imaging/time-series/uploaded-file lists. Cached catalogs,
  live Villanova, mirrored VSX, finite SIMBAD and arXiv checks produce no
  identification within their stated scope. This supports continued
  investigation, not a discovery claim.

## Current stellar and Gaia measurements

The current ExoFOP response cites TIC v8.2, dated2019-04-15. It provides
temperature5588±122K, radius1.62747 solar radii, mass0.99 solar masses,
logg4.01067, and distance940.079±54.1765pc; no radius or mass error is
supplied. Its contamination ratio is the string **“none”**, not a measured
zero. The old sector-1 QLP header radius1.05189, temperature5565K, and
TESS magnitude13.326 should not be treated as current host measurements.

Gaia GSP-Phot gives temperature5244.236K (reported bounds5118.792–5341.797),
logg4.2753, and radius1.1929 solar radii (bounds1.1350–1.2707).
FLAME radius is1.6339417 and mass0.9527874 solar units. These are different
catalog estimates with different modeling assumptions, not independent
spectroscopic validation; radius/model disagreement and possible blending
need explicit propagation into any physical fit.

Target Gaia G is13.965893, BP−RP1.0212107, parallax1.097247±0.051539mas,
and its astrometric solution is five-parameter (`astrometric_params_solved=31`).
Besides the warnings above, Gaia reports duplicated_sourcefalse,
non_single_star0, IPD harmonic amplitude0.0102448, BP/RP excess factor1.24368,
and360 G observations with mean flux/error1154.98. The non_single_star0
flag means no published NSS processing flag here, not exclusion of a binary.
IPD multi-peak fraction22% is a source-image diagnostic rather than a
time-series eclipsing-binary classification.

The RV expected S/N is8.90195 and template temperature5000K. The low-S/N,
imprecise median lacks orbital/phase information and cannot be interpreted
as planetary or stellar reflex amplitude. Target phot_variable_flag is
NOT_AVAILABLE; no Gaia variable-classifier or eclipsing-binary solution
row was found. Those nulls do not establish constant flux or rule out a binary.

## Photometric interpretation requiring physical validation

The existing sector-1 vetting source reports refined period2.320193020765d,
epochBTJD1327.608120440063, local raw-SAP S/N17.66 and seven qualified
events. Its fitted exposure-integrated trapezoid has duration
**0.0755775±0.0092944d =1.814±0.223hours**, ingress fraction0.30097,
and depth0.0154184±0.0013132. The local KSPSAP depth is0.01401693.
These measurements demonstrate a substantial feature in that photometric
diagnostic, but the flexible trapezoid is not a limb-darkened transit or
binary fit, and its diagonal/covariance errors are approximate.

The illustrative proxy `Rstar * sqrt(depth)` is about**21.0 Earth radii
(1.88 Jupiter radii)** with the TIC radius, or**15.4 Earth radii
(1.37 Jupiter radii)** with GSP-Phot radius. This assumes an undiluted
event on the catalog target and appropriate nongrazing geometry; it is not
a companion-radius posterior. A larger grazing stellar companion can produce
a shallow short event, and unresolved blending can change both the apparent
depth and inferred host radius.

Using TIC mass0.99, radius1.62747 and the same simple depth proxy, a circular
central-transit duration is about**4.44hours**. Using GSP-Phot radius1.1929
with FLAME mass0.9528 gives roughly**3.28hours**, a mixed-catalog
illustration rather than a self-consistent stellar solution. Both exceed the
1.81-hour fitted width. The nominal short duration plus TIC radius implies
near-limb/grazing geometry in a simple circular estimate, where the depth-size
proxy becomes unreliable; eccentricity and stellar uncertainties remain
alternatives. The source must be localized before a planet/EB fit is assigned
to a particular star.

The separate current multi-sector investigation includes the discovery sector
plus16 additional QLP sectors and reports recurrence near2.3218days.
That is scientifically useful; this catalog report does not independently
reprocess those observations or adopt their planetary assessment. Recurrence
alone cannot distinguish a planet from a coherent binary or blended source.

## Gaia/TIC neighborhood

The current120arcsec Gaia cone returns29 sources including the target.
The MAST `Mast.Catalogs.Tic.Cone` query also returned29 TIC rows and linked
the important Gaia sources to TIC counterparts.

| Neighbor | Separation arcsec | TESS mag | Full TESS-flux ratio to target | Gaia/source context |
|---|---:|---:|---:|---|
| TIC650145356 / Gaia4629516351818897408 | 10.73 | 19.1318 | 0.00520 | G19.3535, closest catalog neighbor; full G-flux ratio0.00700 |
| TIC238225943 / Gaia4629516283099411328 | 33.32 | 17.7598 | 0.01843 | Very red Gaia BP−RP2.3204; possible eclipse contaminant |
| TIC238202212 / Gaia4629516321755096448 | 35.20 | 15.9500 | 0.09756 | G16.4810; potentially sufficient flux for a diluted stellar eclipse |
| TIC238202215 / Gaia4629516149956404480 | 96.66 | 13.5422 | 0.8963 | G14.0543; RUWE0.9170, current TIC radius1.1279 |
| TIC238225954 / Gaia4629517902303061632 | 115.37 | 13.4381 | 0.9864 | G14.0163; RUWE10.8034, excess noise1.67037mas/significance1356.48 |

These are full catalog flux ratios, **not measured aperture fractions**.
The nearest faint star's full catalog TESS flux is below the1.4% target
depth, but this simplified flux comparison does not replace spatial
photometry or account for passbands, aperture coupling, source variability
and extraction/background effects. The33- and35-arcsec sources are plausible
blend alternatives if their light contributes sufficiently. The nearly equal
distant stars are strong reasons to inspect broad PSF wings and common
neighbor light curves rather than assume the short ExoFOP list is complete.

TIC238202215 has Gaia radius GSP-Phot1.1424/FLAME1.15957, parallax
1.403674±0.013874mas, medianRV5.92712±3.40579km/s from23transits,
without RV variability metrics. Its current TIC radius is1.1279, mass1.0
and temperature5656K.

For TIC238225954, Gaia GSP-Phot radius1.5761 differs greatly from FLAME
0.855813; TIC supplies no temperature/radius/mass. Its parallax
2.07158±0.15811mas and extremely poor astrometric fit warrant caution about
any single-star estimate. Gaia medianRV55.01430±3.53023km/s over22transits
likewise has no variability statistics. These facts do not tie its variability
or any possible companion to the target ephemeris.

Live ExoFOP responses for TIC238202212,238202215 and238225954 have empty
TOI/CTOI/confirmed and uploaded-file lists. Live Villanova explicitly lists
none as a catalog member. Their absence does not exclude an unpublished EB.
The target's own unresolved Gaia source structure is distinct from these
resolved distant neighbors and remains an independent concern.

## Finite prior-recognition and novelty scope

The cached `CatalogIndex(position_radius_arcsec=120)` returned **no exact
host or nearby catalog object** for both frozen coarse periods2.32229783
and4.64135203d. It checked exact TIC, positional host/neighbor associations,
and2–5 integer harmonic/subharmonic matches within1% using NASA confirmed
planets of all discovery methods (including RV), TOIs, CTOIs, DIAmantev2,
Prsa2022 EBs, Olmschenk2021 and Yu2019. The approximately4.64-day feature
is close to twice the2.32-day feature and is not evidence of a separate
planet merely because it appears as another peak. Per-period outputs and
all underlying snapshot dates/hashes remain in the cached-match JSON.

Current target ExoFOP planet-related and prior-observation/file lists are
empty. Live Villanova explicitly says “This TESS ID is not in the catalog.”
The SIMBAD120arcsec cone is completely empty, so no target-associated
SIMBAD aliases or references were available. VSX's mirrored
Version2026-10-05 contains no variable-star rows in the cone. Finite arXiv
metadata queries for TIC238202207, 2MASSJ03532906-7347087,
GaiaDR3 4629516356114834944, UCAC4082-005492 and APASS33556210 return
zero entries. They search indexed metadata, not complete paper text.

A finite MAST CAOM query with `obs_collection=TESS` and exact bare
`target_name=238202207` returned `COMPLETE` with zero observation rows.
The response and full request are retained in
`tic238202207_mast_target_observations.json` and its UTC/hash manifest.
This tests that exact index key only: it cannot exclude public FFI-based TCEs,
reports or differently indexed products, and is not an exhaustive TCE search.

**No prior identification was found within this finite scope.** Remaining
requirements include full SPOC/QLP TCE/DV lists, recent journal supplements,
Planet Hunters and other candidate/EB surveys, full-text/all-alias literature
checks, and independent source attribution. Catalog absence cannot establish
that this is previously unrecognized or planetary.

## Preservation, validity notes, and reproduction

New public bodies, exact URLs or POST request/ADQL, actual UTC retrieval,
HTTP status, sizes and SHA-256 manifests are stored under
`data/broad_catalogs/tic238202207_*`. Principal suffixes are
`cached_matches_current.json`, `exofop_current.json`,
`exofop_nearby_current.csv`, `gaia_target_current.csv`,
`gaia_stellar_current.csv`, `gaia_class_current.csv`,
`gaia_eb_current.csv`, `gaia_neighbors_current.csv`,
`gaia_neighbors_derived_current.csv`, `gaia_bright_neighbors_current.csv`,
`tic_cone_current.json`, `tesseb_live.html`, `vsx_current.tsv`,
`simbad_current.csv`, `arxiv_current.xml`, and
`neighbor{ID}_{exofop_current.json,tesseb_current.html}`.
The additional `mast_target_observations.json` response and request are also
retained. All23 source/derived body hashes passed verification against their
manifests.

The first ExoFOP response for neighboring TIC238202215 was HTTP200 with
an **empty body**, so it is unusable for any recognition inference. It was
preserved unchanged. A fresh retry produced valid JSON at
`tic238202207_neighbor238202215_exofop_current_20261014T102913Z.json`
with its own manifest; only that valid retry supports the neighbor statements
above. HTTP success alone is not scientific schema validation.

Invoke `handle(238202207)` from
`scripts/secondary_leads_catalog_check.py` for repeatable current target
Gaia/ExoFOP and cached-index lookups. Run
`PYTHONPATH=. /workspace/openrouter-check-venv/bin/python scripts/tic238202207_source_context.py`
for derived neighborhood, MASTTIC cone, and bright-neighbor Gaia lookups.
Other checks can be replayed using the exact public URLs/requests in their
manifests, with TLS verification enabled and new dated output names. The
survey SQLite was opened in mode=ro. Existing scientific snapshots, frozen
rank methods/tables, source catalogs, and photometric measurements were
preserved. No paid services, credentials or API keys were used.

