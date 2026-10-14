# TIC 139148478: plausible dwarf, major rotational-neighbor contamination

Checked 2026-10-14 for the approximately 0.701686-day exploratory box feature.
The target's current stellar parameters support a dwarf, but the photometry
has major source-contamination concerns. **No credible new planetary claim
follows from this catalog check.** A known rotational variable lies 51.5 arcsec
away and is substantially brighter. Its measured period differs from the
box feature; the catalog evidence does not identify the source of that feature.

## Current stellar and astrometric context

| Parameter | Current TIC v8.2 / ExoFOP | Gaia DR3 |
|---|---:|---:|
| Stellar radius solar units | 0.717793 ± 0.043900 | GSP-Phot 0.7316; FLAME 0.736508 |
| Stellar mass solar units | 0.87 ± 0.107199 | FLAME mass not supplied |
| Temperature K | 5149 ± 127.634 | GSP-Phot 5148.56 |
| log g | 4.66559 ± 0.083133 | GSP-Phot 4.5969 |
| TESS magnitude | 12.1213 | — |
| Gaia G | ExoFOP 12.6663 | 12.65315 |
| Distance pc | 215.504 ± 3.322 | Parallax 4.60345 ± 0.01373 mas |

The old QLP header radius is −1, an invalid sentinel. It must not be used as
a physical radius. Current independent estimates are consistent with a dwarf.
Gaia source 6540946402503401216 has RUWE 0.94105, five-parameter astrometry,
IPD multiple-peak fraction 0%, excess noise 0.04618 mas, and no NSS flag.
It has no returned Gaia rotation, variable-classifier, or EB solution row;
`phot_variable_flag=NOT_AVAILABLE` is not a photometric-quietness measurement.

Gaia target median RV is −8.46759 ± 0.93428 km/s over 30 transits. Its robust
RV amplitude is 17.2748 km/s, but the reported constant-RV chi-square
probability is 0.57417 and renormalized goodness of fit −0.67907. The
amplitude alone is therefore not supported evidence of RV variability or an
orbital companion. No planetary mass or useful binary exclusion follows.

## Known bright rotational neighbor

TIC 139148479 = HD 217027 = HIP 113422 is approximately 51.5 arcsec away,
at Gaia position RA 344.56635157, Dec −44.47286090 degrees. It has TESS
magnitude 9.3069, compared with target 12.1213: approximately **13.36 times
the target's full catalog TESS flux**, before aperture/PSF coupling.
Its distance is 59.5951 ± 0.21615 pc and it is a separate foreground star.

Gaia identifies the neighbor as variable with rotation modulation, and its
classifier reports `SOLAR_LIKE` with score 0.49717. The Gaia rotation solution
gives best period **9.50244 ± 1.33408 days**, with segment periods 13.2050,
8.0863, 9.6854, 16.1865, and 10.7356 days. VSX Version 2026-10-05 lists the
same source as `ROT` with period 9.5024 days and approximately G 9.840–9.880.
SIMBAD associates HD 217027 with Gaia's solar-like rotation publication,
2023A&A...674A..20D, [doi:10.1051/0004-6361/202244178](https://doi.org/10.1051/0004-6361/202244178).
Its other references include RAVE, nearby-star, and acceleration catalogs;
their presence does not itself establish a short-period companion.

These rotational periods do **not** directly recognize the target's
0.701686-day box feature. Quasi-periodic neighbor variability is an empirical
nuisance hypothesis, rather than a demonstrated identification of every dip.
The independent neighbor light curves were retained by the parent analysis
with retrieval information in `data/broad_followup_queries/tic139148479_qlp.json`.
Direct waveform/correlation/source-localization tests should determine its
actual contribution instead of assigning the feature to it from proximity.

Current ExoFOP gives the target contamination ratio **0.995876**, approximately
equal predicted contaminating and target flux for the catalog aperture model.
This is not a measured QLP aperture fraction and is not a source-localization
test. Nevertheless, interpreting this as an isolated dwarf is unsafe. The
neighbor's full catalog flux is likewise not all collected in the aperture.
Dilution, aperture choices, and detector geometry can change measured depths.

## Finite prior-recognition scope

Current ExoFOP lists no target TOI, CTOI, confirmed planet, planet parameters,
or uploaded follow-up files. The current live Villanova exact-TIC detail states
that the target is not in its EB catalog. No target or neighboring prior
planet/candidate/EB source appears within 120 arcsec in the seven prepared
catalogs at the examined period. No Gaia EB solution was returned for either
target or bright neighbor. No dedicated complex-periodic-variable literature
catalog was comprehensively searched; absence of a Gaia class does not rule
out that phenomenology.

SIMBAD's cone contains HD 217027 and background galaxies, without a
target-associated star entry. The target ExoFOP aliases are TIC 139148478,
2MASS J22581648-4427311, APASS 22080247, Gaia DR2/DR3
6540946402503401216, UCAC4 228-201127, and WISE J225816.48-442731.2.
Finite arXiv metadata searches of these aliases return zero entries.
The ASAS-SN archived cone has no source rows; VSX returns the already-known
bright rotational neighbor. These are not complete full-text, supplementary
table, TCE, or all-alias novelty checks. Nonrecognition is not discovery.

## Photometric interpretation and next tests

The final detailed analysis reported roughly 29 sampled local events,
formal raw/processed S/N about 8.7, and depth near 3800 ppm. The raw sector-1
variability spans about 10%, and measured QLP aperture depths vary from
approximately 2387 to 3800 ppm. A tiny QLP y-centroid offset has formal
significance near 6.8 sigma. These are screening/diagnostic statistics without
a global trials or robust centroid-systematics calibration.

The parent multi-sector analysis reports weak sector-28 recurrence: local
raw S/N around 4, joint/local statistic around 3, and depth about 1454 ppm
versus approximately 3472 ppm in the first sector. That photometric check is
recorded in `candidates/broad_s1/tic139148478_multisector/`; it is independent
of the recognition checks here. Aperture/coupling and ephemeris uncertainty
must be considered before treating weak recurrence as an exact exclusion,
but it fails to strengthen the planetary interpretation.

Required work before promotion would include source-localized independent
recurrence, neighbor-template contamination tests, phase/event robustness,
rotation/complex-waveform hypotheses, and a trials/red-noise assessment.
Current evidence supports keeping this as a weak unexplained, contaminated
photometric lead rather than a credible new candidate.

All public responses and URL/UTC/hash manifests are retained under
`data/broad_catalogs/tic139148478_*`. No frozen ranking, survey database,
raw observation, or current-star table was changed by these checks. No paid
external API calls were used.
