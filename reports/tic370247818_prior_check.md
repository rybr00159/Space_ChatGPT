# TIC 370247818: recognition, stellar context, and contamination checks

Checked 2026-10-13. Classification from these checks: **unrecognized
transit-like signal with incomplete novelty and false-positive assessment**.
There is no confirmed discovery claim. The approximately 3.523-day photometric
lead requires independent observations and spatial localization.

## Recognition checks

The seven prepared catalog sources find no entry at the target TIC or within
120 arcsec of its position (RA 55.45564053, Dec −79.11323158 degrees), including
no exact or 2–5 integer harmonic/subharmonic match to 3.5236501678 days.
Sources: NASA confirmed planets of all discovery methods, TOIs, ExoFOP CTOIs,
DIAmante v2, Prša 2022 published TESS EBs, Olmschenk 2021, and Yu 2019.
Exact snapshot metadata and the empty-match evidence are saved in
`data/broad_catalogs/tic370247818_cached_matches.json`. Catalog absence is not
proof that a signal has never been reported.

The current public ExoFOP target JSON has empty TOI, CTOI, confirmed-planet,
planet-parameter, spectroscopy, time-series, imaging, and uploaded-file lists.
The current live Villanova exact TIC detail states, "This TESS ID is not in
the catalog." Its exact search returns zero objects. This is a genuine
membership check rather than treating a light-curve link as EB membership.
An unlisted or unresolved eclipsing binary remains possible.

SIMBAD's two-arcsec cone returns no object; its 120-arcsec cone returns only
the brighter neighboring star CD-79 151. Consequently no target-associated
SIMBAD bibliography was available. The ExoFOP aliases are:

- TIC 370247818
- 2MASS J03414936-7906478
- APASS 34591811
- Gaia DR2 / DR3 4625697305555206912
- UCAC4 055-002273
- WISE J034149.35-790647.6

An arXiv metadata query using the numeric TIC, the 2MASS designation component,
Gaia source identifier, UCAC4 designation component, and APASS number returns
zero results. This searches indexed metadata, not all paper full text,
supplementary tables, survey databases, or every alias. It is not a complete
literature novelty search. The archived VSX Version 2026-10-05 cone query
within 120 arcsec returned schema information but no variable-star data rows.
The ASAS-SN CDS catalog `II/366/catv2021` (table version 2021, CDS metadata
updated 2024-09-11) also has no variable-star data rows in the 120-arcsec cone.
The successful query selected the catalog collection `II/366`; an initial
incorrect `II/366/catalog` path returned a table-path error and is explicitly
marked unusable in its manifest. Neither the error nor the older ASAS-SN table
is a complete current variability veto.

All requests were read-only public requests; no paid API was used.

## Updated stellar parameters

The old QLP header radius (1.31128 solar radii) should not be the sole stellar
radius used in interpretation. Current ExoFOP/TIC v8.2 and Gaia estimates are
larger and mutually consistent.

| Parameter | TIC v8.2 through ExoFOP | Gaia DR3 |
|---|---:|---:|
| Effective temperature K | 6383 ± 137 | GSP-Phot 6344.53 |
| Stellar radius solar radii | 1.68401 ± 0.08391 | GSP-Phot 1.7453; FLAME 1.73505 |
| Stellar mass solar masses | 1.28 ± 0.20331 | FLAME 1.33971 |
| log g | 4.09259 ± 0.09212 | GSP-Phot 4.0569 |
| Distance pc | 760.545 ± 14.0085 | Parallax 1.314279 ± 0.010765 mas; no distance inversion used here |
| TESS magnitude | 12.2996 | — |
| Gaia G | ExoFOP 12.6819 | 12.67596 |

Gaia target RUWE is 1.06193, `non_single_star=0`, and photometric variability
is `NOT_AVAILABLE`. No row is present in the Gaia DR3 variable classifier or
eclipsing-binary solution table. These nulls do not exclude a binary or planet;
they indicate absent Gaia classifications/solutions in the queried tables.
The source has a five-parameter astrometric solution (`astrometric_params_solved=31`),
IPD multiple-peak fraction 0%, excess noise 0.07182 mas (catalog significance
5.854), BP/RP excess factor 1.19708, and `duplicated_source=false`. Its G-band
flux/error is 3822.15 from 348 observations. No severe Gaia multiple-image
quality warning comparable to the second lead TIC 141980042 appears in these
fields; this still cannot spatially rule out a TESS contaminant or an
unresolved companion below Gaia's effective detection limits.
Gaia supplies target median radial velocity 13.6498 ± 3.6838 km/s from 15
transits, but no robust amplitude or RV variability probability/fit statistic.
The reported median and its imprecise error do not supply a planetary mass or
exclude an eclipsing companion. For bright CD-79 151, Gaia reports median
−12.7460 ± 0.2227 km/s, robust amplitude 1.4274 km/s over 12 transits,
constant-RV chi-square probability 4.37e−5, and renormalized goodness of fit
3.7597. This indicates RV variability in that catalog metric, without an
orbital period, and does not associate it with the target's 3.52-day dips.
No Gaia median RV is supplied for TIC 370247802.

The simple radius proxy `Rstar * sqrt(depth)` is approximately 14.3–15.5 Earth
radii for the refined/screened 0.608–0.716% depths and TIC v8.2 radius. It
assumes negligible dilution and a sufficiently nongrazing transit. With mass
1.28 and radius 1.684 solar units, the approximately 3.5225-day circular
central-transit duration is roughly 4.6 hours, including the proxy companion
radius. The refined BLS box duration is only 1.76 hours. A rough circular
geometry implies impact parameter near 1 and possible grazing, where
`sqrt(depth)` can underestimate companion size. A larger stellar companion
can remain consistent with a shallow grazing event. Exposure-integrated
transit fitting and explicit grazing/binary alternatives are essential;
30-minute exposures provide weak ingress/egress resolution.

## Contamination: extend beyond the ExoFOP nearby list

The ExoFOP nearby CSV includes five neighbors at approximately 31–51 arcsec,
all TESS magnitudes 17.83–19.49. Its finite cone omitted the much brighter
star at 86.64 arcsec. A fresh Gaia 120-arcsec cone returns 27 sources, and a
MAST TIC cone identifies their TIC counterparts. Important neighbors are:

| Neighbor | Separation arcsec | TESS magnitude | Relevant context |
|---|---:|---:|---|
| **TIC 370247831 = CD-79 151** | **86.64** | **9.0714** | Gaia 4625697168116255360, G 9.71409, RUWE 3.14106. About 19.6 times the target's catalog TESS flux. |
| TIC 370247802 | 109.05 | 13.8027 | Gaia 4625697717872066432, G 14.73775, RUWE 54.5749. About 25% of the target's catalog TESS flux. Very poor astrometric fit merits scrutiny. |
| TIC 370247833 | 116.22 | 16.202 | Gaia 4625696480921489920. About 2.7% of target catalog TESS flux; possible off-aperture contaminant. |
| TIC 370247816 | 30.86 | 17.9571 | Faint neighbor near the central aperture boundary. |
| TIC 370247820 | 44.79 | 17.8266 | Brightest of the five neighbors in the shorter ExoFOP list. |

Relative fluxes above use catalog TESS magnitudes, without PSF/aperture
coupling. Gaia G ratios are also saved as clearly labeled G-band proxies and
must not be substituted for TESS aperture fractions. Gaia astrometric RUWE is
not an EB diagnosis.

TIC v8.2's contamination ratio is 0.007605, comparable to the signal depth.
This catalog ratio is a model prediction for an assumed aperture, not a
measurement of this particular QLP aperture or a spatial veto. Bright-star
PSF wings, unmodeled faint sources, eclipsing neighbors, and background
subtraction must be considered. Conversely the full catalog flux of a distant
neighbor is not all collected in the target aperture. A time-series/difference
image test is required to establish or reject actual contamination.

CD-79 151 has no match in the prepared planet/EB/prior-candidate index. Its
SIMBAD bibliography points to Hon et al. 2021, ApJ 919, 131,
[doi:10.3847/1538-4357/ac14b1](https://doi.org/10.3847/1538-4357/ac14b1), the
158,505 oscillating-red-giant QLP catalog. This is evidence of prior giant-star
analysis, not evidence of a 3.52-day eclipse. Aliases include
2MASS J03413878-7908090, TYC 9371-522-1, SAO 256017, CPD-79 115,
CPC 21.3 410, GSC 09371-00522, and PPM 368145. Independent neighbor
photometry at the candidate ephemeris is an immediate falsification test.
The current stellar/recognition lookup for both important neighbors is retained
in the candidate-prefixed JSON source files.
Current ExoFOP reports radius 9.84042 solar radii and temperature 4880 ± 122 K
for CD-79 151, consistent with the oscillating-giant bibliography. Its live
Villanova detail also states it is not in the EB catalog, and current ExoFOP
TOI/CTOI lists are empty for it and TIC 370247802. No stellar radius or mass
is supplied for the high-RUWE TIC 370247802. These lookups do not establish
whether either neighbor varies at the candidate's photometric period.

## Evidence and reproduction

Public responses, exact retrieval URLs or MAST service/parameters, UTC
retrieval times, and SHA-256 manifests are under
`data/broad_catalogs/tic370247818_*`. Existing source catalog snapshots were
not overwritten. The small derived Gaia neighbor table includes separations
and band-labeled flux ratios. The candidate photometric diagnostic source is
`candidates/broad_s1/tic370247818_p3.523650/vetting.json`.

For rechecking, use each public `url` in the corresponding manifest with TLS
verification enabled; the MAST TIC cone source records endpoint, service, and
parameters. Current catalog responses can change, so retain a new dated
snapshot rather than replacing the recorded source/hash. Recompute the
cached broad match with `CatalogIndex(position_radius_arcsec=120)` and the
recorded TIC, period, RA, and Dec. A current match report cannot be substituted
for the original dated snapshot.

Outstanding recognition checks include full SPOC/QLP TCE/DV listings, recent
journal supplementary tables, Planet Hunters reports, all alias variants,
and non-arXiv full-text literature. Outstanding scientific checks include
independent sectors, common-mode artifacts, exposure-integrated transit/binary
fits, centroid/difference images, neighbor light curves, and a search-null
significance assessment accounting for the broad number of stars/periods.
At this stage no credible novelty or planetary claim is supported by catalog
absence alone.
