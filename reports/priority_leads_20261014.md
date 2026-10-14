# Current stellar context falsifies four apparent small-star leads

Checked 2026-10-14 for TIC 441386747, 406280158, 271697430, and 260352037.
**All four targets are classified as giants by the current public TIC, with
catalog radii 10–12 solar radii. Their old QLP header radii below one solar
radius must not support small-planet interpretations.** Giant variability and
detrending artifacts are appropriate hypotheses to test. A signal from a
different blended star remains possible and requires spatial localization.
These checks do not independently prove the origin of every photometric dip.

## Stellar checks

| TIC | QLP R solar | Current TIC R solar | Gaia GSP-Phot R solar | SAP local depth ppm | Radius proxy Earth radii | Nominal absolute G |
|---|---:|---:|---:|---:|---:|---:|
| 441386747 | 0.733908 | 12.1823 | not supplied | 1363 | 49.1 | 0.64 |
| 406280158 | 0.781807 | 10.5196 | not supplied | 1132 | 38.6 | 0.45 |
| 271697430 | 0.742671 | 10.2269 | 11.6443 | 818 | 31.9 | 0.83 |
| 260352037 | 0.781581 | 10.1065 | not supplied | 1123 | 36.9 | 0.68 |

Current ExoFOP/TIC v8.2 and exact-ID MAST TIC rows agree on the large radii.
The MAST `lumclass` field is `GIANT` for all four. Radius uncertainties and
masses are absent in these TIC rows and must not be invented. Gaia parallaxes
are 1.10323 ± 0.01646, 1.48083 ± 0.01352, 1.37068 ± 0.01223, and
0.59728 ± 0.00959 mas, respectively. The displayed nominal absolute G uses
`G + 5*log10(parallax_mas) - 10`, without extinction or parallax-zero-point
correction; it is a luminosity-context proxy, not a full stellar inference.
These bright absolute magnitudes and K-like colors support the evolved-star
context rather than the old dwarf interpretation.

Radius proxies use current TIC radius times `sqrt(SAP local depth)` and the
solar/Earth radius conversion, without dilution or grazing corrections.
They are not fitted companion radii. All are substantially above the broad
queue's 22-Earth-radius warning threshold. A grazing or diluted eclipse can
change interpretation; it cannot validate the old isolated-small-star model.

Gaia's independent GSP-Phot radius 11.6443, log g 2.698, and FLAME radius
12.70247/mass 2.70437 solar units additionally support a giant for
TIC 271697430. At the approximately 3.059-day period, Kepler's-law semimajor
axis is near or within these inferred giant radii, depending on which
radius estimate is used. The catalog-host orbital model is inconsistent or
close to the photosphere; parameter errors and potential off-target sources
remain caveats. Missing TIC masses prevent applying this calculation to all
four through a uniform current-TIC mass model.

## Literature and finite recognition

SIMBAD associations are CD-34 14539, CPD-62 6263, TYC 9191-104-1, and
TYC 8542-438-1. The first three are associated with Hon et al. 2021,
ApJ 919, 131, [doi:10.3847/1538-4357/ac14b1](https://doi.org/10.3847/1538-4357/ac14b1),
the QLP catalog of 158,505 oscillating red giants. TIC 271697430 also has
TESS-HERMES/APOGEE and red-giant spectroscopic-age references. TIC 441386747
has RAVE references. These are independent prior stellar-analysis context,
not prior exoplanet claims. No SIMBAD reference rows were returned for
TYC 8542-438-1 in the checked query.

All four current ExoFOP TOI, CTOI, confirmed-planet, and planet-parameter
records are empty. The live Villanova exact-TIC details state that each TIC
is not in the catalog. The seven prepared planet/prior-signal/EB catalogs
have no target or positional entry within 120 arcsec at the examined periods.
Finite arXiv metadata queries using the supplied ExoFOP aliases return zero
results; this is not full-text or supplementary-table novelty checking.
No discovery inference follows from these absences.

VSX Version 2026-10-05 has no cone rows for the last three. The first target's
cone contains the neighboring variable Gaia DR3 6791599491309064192 =
TYC 7464-742-1, at approximately 62.45 arcsec, G 12.30665, VSX period
0.51882 days and type `VAR`. That is a distinct star and a distinct reported
period, not a recognition of the target's 7.48-day box peak.

## Gaia quality, RV, and contamination caveats

All four targets have Gaia RUWE near one (1.0448, 1.0393, 0.9939, 0.9695),
IPD multi-peak fractions zero, and no queried Gaia EB membership. These
facts do not make giant-star troughs planetary. Gaia robust RV amplitudes
are approximately 1.24, 1.97, 1.28, and 5.39 km/s. TIC 406280158 has
constant-RV chi-square probability 1.30e−5 and renormalized goodness of fit
2.53, indicating RV variability in that metric without an identified orbital
period. For the other three, the catalog RV variability probabilities do not
provide analogous evidence; an amplitude value alone is not a binary diagnosis.

The two short-separation ExoFOP entries TIC 441392404 (5.88 arcsec) and
TIC 260352038 (3.96 arcsec, apparently brighter than its target) are marked
**`disposition=ARTIFACT`** in current exact-TIC MAST queries and have no Gaia
identifiers. They must not be treated as established contaminating stars from
the nearby-table magnitude alone. Their exact current rows are retained in
`priority_leads_20261014_close_tic_entries.json`.

Actual Gaia-resolved neighbors still warrant scrutiny. Examples include
Gaia 6791599182071420416, G 15.42149, 11.92 arcsec from TIC 441386747;
TIC 271697433, TESS magnitude 12.7903, approximately 37 arcsec from
TIC 271697430; and a near-equal-G-brightness source
5496091252930328576, G 11.46783, 99.32 arcsec from TIC 260352037.
The saved 120-arcsec Gaia cones contain all checked nearby sources.
Catalog flux ratios do not measure PSF/aperture coupling and do not establish
which star produces a dip. Current ExoFOP contamination ratio `none` is
unknown, not zero.

## Reproduction and effect on the broad survey

Public responses and URL/hash/time manifests are preserved under
`data/broad_catalogs/tic{ID}_*` and
`data/broad_catalogs/priority_leads_20261014_*`. Derived per-target SIMBAD
rows reference the complete queried batch and its checksum. Exact current
MAST service/parameters are retained, including the multi-ID proof query.
Original QLP files, survey SQLite, and catalog snapshots were unchanged.

This evidence motivated the executed all-strict-eligible stellar refresh:

```sh
python src/refresh_shortlist_stars.py --include-vetted --batch-size 128
python src/broad_prioritize.py \
  --stellar-table data/broad_catalogs/current_stars_20261014/current_stellar.csv \
  --catalog-radius-arcsec 120
```

The completed survey has 159,505 actually searched stars, 495 insufficient
series, and 319,496 saved peaks. Current TIC rows were obtained for all 3,407
strict-eligible target TICs plus three additional historical-vetting TICs:
3,410 exact matches, 1,541 evolved/giant contexts, and 270 missing radii.
Earlier tables and the old-header queue were preserved before replacement.

The refreshed queue uses current values or explicit missing values, preserves
old header columns, separates evolved/large-radius warnings, and applies an
explicit photosphere-orbit inconsistency warning when current mass is
available. Nearby recognition covers 120 arcsec without assigning the source
or blanket-discarding known hosts. The retained unrecognized queue is 200
signals on 187 distinct stars, with 22 signal rows lacking current radii and
28 lacking masses. This queue remains exploratory: its scores have no
false-alarm-probability interpretation, and every surviving lead requires
independent photometric, source, and literature validation.
