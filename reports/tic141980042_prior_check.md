# TIC 141980042: recognition and blend warning

Checked 2026-10-13 for the approximately 2.321817-day transit-like lead.
**The signal is unrecognized by the checked catalogs, but the unblended planet
interpretation has substantial counterevidence from Gaia image/astrometric
quality.** A photometric dip at this TIC is not yet attributable to a single
well-characterized star. No discovery or credible planetary claim is made.

## Recognition scope

- The seven prepared public planet/prior-signal/EB sources have no TIC or
  positional match within 120 arcsec of RA 96.36705, Dec −72.649775 degrees.
  This includes confirmed planets of all discovery methods, TOI, CTOI,
  DIAmante v2, Prša 2022 EBs, Olmschenk 2021, and Yu 2019. Empty match evidence
  and snapshot hashes are retained in `tic141980042_cached_matches.json`.
- Current live ExoFOP JSON has no confirmed planets, TOIs, CTOIs, or planet
  parameter records. Its contamination ratio is the string `none`, which
  does **not** mean a measured contamination ratio of zero.
- Current Villanova exact TIC detail states, "This TESS ID is not in the
  catalog." Membership absence does not rule out an unreported EB.
- SIMBAD's 120-arcsec cone has no object rows, so no target-associated
  bibliography was available through that positional query.
- An arXiv metadata query of numeric TIC, 2MASS component
  `J06252809-7238591`, Gaia identifier `5265261465421496832`, UCAC4 component
  `087-011411`, and APASS `33617578` returns zero results. This does not search
  all paper full text or all supplementary data and does not establish novelty.
- VSX Version 2026-10-05 and the archived ASAS-SN `II/366` catalog cone queries
  have schema output but no variable-star data rows within 120 arcsec. Gaia
  DR3 EB and variable-classifier queries return no target rows. These are
  incomplete variability/recognition checks, not false-positive vetoes.

## Stellar assumptions fail a Gaia quality check

The old QLP header gives Rstar 1.152 solar radii and Teff 5826 K. In the current
TIC v8.2 data supplied by ExoFOP, radius, mass, effective temperature, and log g
are all absent. Gaia DR3 astrophysical parameters also supply no temperature,
gravity, radius, or mass for this source. A radius-derived planetary size
using the old header is therefore unsupported by a current independent
stellar characterization.

Gaia source 5265261465421496832 has:

| Quantity | Value | Interpretation limit |
|---|---:|---|
| Astrometric parameters solved | 3 | Position-only solution; no DR3 parallax or proper-motion solution. |
| RUWE | absent | Cannot be interpreted as a normal astrometric fit. |
| Astrometric excess noise | 7.09764 mas | Large image/position residuals. |
| Excess-noise significance | 35,648.156 | Catalog evidence of severe astrometric residuals. |
| IPD fraction with multiple peaks | 71% | Strong evidence of multiple-image/crowding problems in Gaia windows; not a companion orbit. |
| IPD goodness-of-fit harmonic amplitude | 0.149989 | Additional image-fit diagnostic, not an orbital measurement. |
| BP/RP flux-excess factor | 1.61562 | Crowding/photometry warning; interpretation is color dependent. |
| Duplicated source | false | This flag alone does not negate the other quality problems. |
| Median Gaia G | 13.44133 | Old DR2/ExoFOP G is 13.1642; release differences need investigation. |
| G-band flux/error | 146.040 | Catalog mean-flux metric, not transit significance. |
| G observations | 270 | Does not guarantee resolved photometry. |

These flags are consistent with a close blend, unresolved multiple stars, or
crowded-image problems. They do not by themselves prove an EB or exclude a
planet in a multiple system. They **do** invalidate treating missing RUWE as
good astrometry or treating the source as an established isolated dwarf.
The ExoFOP distance 527.273 pc has no listed uncertainty and cannot substitute
for the absent DR3 parallax without investigation of its original derivation.

Gaia median RV is 3.94144 ± 6.57660 km/s from ten transits, without a robust
amplitude or RV-variability probability/fit metric. This low-precision median
does not establish a planetary mass, a binary period, or a mass upper limit.

## Nearby stars can supply the observed dip

A 120-arcsec Gaia cone returns 138 sources. The shorter ExoFOP nearby CSV
contains several sources that have enough nominal TESS flux to produce an
approximately 0.84% blended eclipse if they are deeply eclipsing and the
aperture collects enough of their light:

| Neighbor | Separation arcsec | TESS magnitude | Nominal flux ratio to target |
|---|---:|---:|---:|
| TIC 141980044 | 31.14 | 17.3288 | 0.01433 |
| TIC 141980047 | 34.46 | 17.7247 | 0.009953 |
| TIC 141980027 | 50.81 | 17.5085 | 0.01215 |

Ratios use target TESS magnitude 12.7196 and `10**(-0.4*delta_mag)`.
They are full-star catalog ratios, not measured aperture fractions. They
identify plausible contaminants, without establishing whether any neighbor
eclipses or contributes enough collected flux. Pixel/PSF/difference-image
tests must evaluate this rather than equating catalog flux with aperture flux.

An additional Gaia source 5265261538440949632 at RA 96.38773330,
Dec −72.63158328, separation 69.16 arcsec, has G 15.08205, approximately
22% of the target's Gaia G flux. MAST identifies it as **TIC 141980051**,
TESS magnitude 14.4044, about 21.2% of target nominal TESS flux.
It is outside the shorter ExoFOP list and
can be an important PSF-wing source. Another source, 5266012195647515648 at
RA 96.31438055, Dec −72.62247219, separation 113.42 arcsec, G 15.64560,
has approximately 13% of target G flux. It is **TIC 141980058**, TESS
magnitude 14.9034, about 13.4% of target nominal TESS flux.
G ratios are color-dependent proxies
and must not be taken as TESS-band aperture fractions. The MAST TIC cone is
retained for exact catalog associations.

The high Gaia multiple-peak fraction also motivates sub-arcsec inspection:
even a clean TESS difference image cannot resolve all potential close blends.
Consistent light-curve depths in neighboring QLP apertures and an
inconclusive secondary-eclipse test do not remove this issue.

## Required next tests and evidence

Prioritize target/neighbor difference imaging, exact independent-sector
recurrence, half/double-period tests, and exposure-integrated transit versus
grazing/diluted-binary models. Obtain a reliable stellar characterization
and higher-resolution source/blend information before assigning a planetary
radius or promoting this lead to a credible new candidate. The Gaia quality
flags justify ranking it below an otherwise comparable well-characterized,
spatially localized target.

All public response bodies and URL/hash/time manifests are preserved under
`data/broad_catalogs/tic141980042_*`; nothing in the original source snapshots
was overwritten. They can be rechecked by retrieving manifest URLs with TLS
verification enabled and saving new dated snapshots if the sources change.
Remaining novelty scope includes complete TCE/DV reports, journal
supplementary tables, Planet Hunters, non-arXiv literature, and all aliases.
No paid external API requests were used.
