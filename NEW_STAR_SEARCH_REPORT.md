# Fresh search for uncatalogued stars and unrecognized companions

Research checkpoint: evening 14 October 2026, America/Halifax. Archive manifests
record UTC, including 15 October. **No genuinely uncatalogued star is established.
Two possible previously unreported wide stellar associations are the strongest
new leads. They involve existing Gaia sources; neither binding nor publication
priority is proven.** Two additional sparse radial-velocity leads remain less
secure. These are new investigations, separate from the two earlier TESS systems.

The previous 159505-star TESS period search did not exhaust public astronomical
data. This phase deliberately changes the observable: direct image sources,
infrared proper motion, common stellar astrometry, and repeated spectroscopy.
It also distinguishes a new sky source from an existing source with a newly
suspected companion. Neither a catalogue omission nor a fitted motion alone
qualifies as a new star.

## Actual search scope

| Track | Data actually examined | Outcome |
|---|---|---|
| Direct imaging | 12-field pilot, then 64 distinct fields spanning G 10–15; original/model/residual g/r/z images, 1035 recurring peak detections |495 Gaia-missing peaks already in deeper catalogues. No validated uncatalogued star. |
| Infrared motion |1946 quality/motion-selected CatWISE rows in 10 successful 5-degree-radius fields, approximately 785 square degrees |250 Gaia-unmatched cases; independent older/optical epoch tests reject or leave unverified the apparent new-source interpretations. |
| Gaia associations |4103 primaries at approximately 1–2 kpc;1086 neighbor rows, 1081 distinct neighbors |11 unique screened pairs;5 known, 6 unrecognized in bounded checks, 2 cleaner leads. |
| Repeat spectroscopy |3430 repeat-RV stars, 13176 quality visits;1799-star stricter subset |39/40 strongest changes already binary-flagged; remaining strong case fails reduction quality. Two weaker sparse-RV cases remain provisional. |

The tracks overlap and count different units; their counts must not be added
into an independent-star total. The 356869 single/duplicate quality catalogue
visits returned by an aggregate count are not 356869 analyzed stellar systems.
The 1799-star strict sample is a subset, not an additional sample. Image peaks
include galaxies and artifacts. The 76 cutouts cover only about 0.026 square
degrees before core/edge cuts, targeted near moving stars. No all-sky search
or global completeness claim is made.

The [machine-readable summary](results/new_stars/search_summary.json),
[ranked leads](results/new_stars/ranked_leads.csv), and complete per-track tables
separate examined data, controls, catalogue metadata and candidate assessments.

## Strongest new association leads

| Internal pair | Gaia DR3 sources | Angular / projected separation | Why retained |
|---|---|---|---|
|1|4977459112829640832 /4977459112829654400|4.034 arcsec /7134 AU at nominal 1.77 kpc|Compatible astrometry, clean Gaia quality, two stellar images on a 1993 plate |
|3|4923799509418222080 /4923799509418221952|4.763 arcsec /4889 AU at nominal 1.03 kpc|Compatible astrometry, clean photometric/astrometric diagnostics, distinct 1996 plate images |

Their colors and absolute magnitudes are broadly compatible with ordinary
stellar components at a common distance. That is not an isochrone age/mass
measurement. Their nominal transverse relative speeds are only 0.035 and
0.242 km/s, but proper-motion errors are material: propagated minimum-binding-mass
16/50/84 percentiles are 0.28/1.14/3.10 and 0.24/0.97/2.65 solar masses.
Those transformed measurement draws are not binding probabilities or mass
measurements. Unknown radial motion, true line-of-sight separation, angular
Gaia error correlations and selection effects prevent a secure bound-orbit claim.
Neither component has a public Gaia RV or NSS solution in the retained queries.

Five other pairs are already Tian 2020 pairs and are rejected as new. Four
of the six remaining associations have compromised photometry, astrometric
warnings or marginal dynamics and rank lower. One DR3/DR2 identifier changes;
recognition was rerun through its actual alias instead of assuming IDs persist.

Three shifted-position samples total 12309 primary trials with zero passing
matches. Under a homogeneous independent Poisson model this gives an approximate
95% upper expectation of one accidental directed pass in the native 4103-primary
search. It is not zero contamination or an individual binding probability;
stellar associations and correlated Gaia systematics can violate the model.
Known-pair retrieval and measured-separation controls succeeded.

**Catalogue scope is crucial:** all six unrecognized pairs were ineligible for
Tian 2020, fall below SUPERWIDE's proper-motion selection, and lie outside
El-Badry 2021's distance selection. Those nulls do not support novelty. WDS and
SIMBAD cones show no binary recognition, but bibliography/catalogue supplements
are not exhausted. The correct label is *possible unreported association among
catalogued stars*, not newly discovered stars.

Read [the pair dossier](candidates/new_stars/companions/README.md),
[quality and selection checks](candidates/new_stars/companions/quality_and_catalog_scope.json),
and [historical images](candidates/new_stars/companions/historical_images.png).

## Why the direct-source searches did not establish a new star

The bright 12-field pilot produced three promising-looking peaks on saturated
detector trails. A less-saturated 64-field expansion found 27 residual peaks
after initial cuts; independent-band width and centroid checks reduced these
to one. The final compact-looking peak at RA 202.550980,Dec−8.576843 lies inside
the native saturation masks in all g/r/z bands (mask value 18462). It cannot
support a reliable new point-source detection. This does not exclude a real
faint object behind the mask.

A separate faint pilot peak at RA 216.428278,Dec+23.612155 remains unidentified,
with descriptive g/r magnitudes≈22.76/22.09 and widths broader than the local
PSF. Independent PS1 images show no convincing detection at either a fixed
position or the bright primary's back-propagated common-motion position.
Deep PS1-stack and SDSS queries also return no counterpart in the checked cone.
Faintness, extension, variable emission or an image artifact remain alternatives;
it is **an unverified image peak, not a credible new-star candidate**.

3832 flat-AB, static-Gaussian injection trials test only the detection stage.
At 4 arcsec, basic recovery falls from 61% at magnitude 17 to 2.8% at 21; at 16 arcsec
the corresponding rates are approximately 99.6%/98.7%. This is not full-pipeline
completeness: PSF subtraction, real stellar colors, motion and later artifact
cuts are omitted. The compact-band test accepts 200/277 detected catalogue-PSF
controls. These data are weak for precisely the close, faint companions most
likely to remain hidden. More low-resolution cutouts alone would not resolve
that limitation. See [image analysis](reports/new_stars/imaging.md).

The infrared track independently tests 250 Gaia-missing sources and 10 matched
controls against 2MASS positions.27 no-Gaia cases have consistent older positions,
but all 27 are already classified STAR in TIC; several are known brown dwarfs.
The 28 optical follow-ups of sources without 2MASS counterparts give only one
apparent motion confirmation. Additional SDSS observations falsify that optical
interpretation: the proposed counterpart is a galaxy moving only 0.126 arcsec
over 6.87 years, compared with 2.64 arcsec predicted from CatWISE. A different blended
infrared mover is not excluded, but no independently confirmed new stellar source
remains. See [infrared analysis](reports/new_stars/infrared.md).

## Spectroscopic leads and attempted falsification

Quality-controlled repeat GALAH DR4 velocities initially yield 40 large-amplitude
leads. The all-star catalogue already flags 39 as binaries; the remaining case
has too few reduction-clean epochs. Requiring an unflagged all-star entry and
clean individual reductions leaves 1799 repeated stars. Seven have catalogue
ranges above 3 km/s: three have published Gaia NSS solutions and two others have
strong published Gaia RV variability. The remaining sources are:

- Gaia 4690783166942834816, UCAC 4 093-001240 /TIC 51939231: catalogue range 7.75 km/s;
  already a VSX rotational variable with reported 20.766-day modulation.
- Gaia 5437472031746145536, TYC 7162-1658-1 /TIC 24731098: range 5.80 km/s; a hotter
  broadened-line star, with pulsation/profile variability among the alternatives.

Most clean observations are clustered in 2014, with one later clean measurement
approximately nine years afterward. Same-field old/new controls have median
offsets−0.105/+0.141 km/s and robust scatters 0.449/0.502 km/s, so a uniform field
zero-point shift does not explain everything. Nevertheless, four downloaded
observed-flux spectra give region-dependent shifts: approximately−7.5..−8.8 km/s
for the first source and−7.2..−10.1 km/s for the second, depending on masking.
These exceed the small formal errors and weaken a simple precise RV interpretation.

The public fit products already use rest-frame wavelengths. The documented
coordinate shift was inverted using matching fit products before the regional
test. This is a same-data consistency check with frame assumptions, not an
independent wavelength-calibration confirmation. Synthetic shifts verify sign
and units only. Three QLP light curves provide a separate pulsation/activity
check, but no common dominant short-period frequency is established across
the hotter star's two sectors. Neither a nondetection nor these detrenders
exclude all pulsations or slow rotation. No orbit, mass function or secondary
stellar mass can be measured from these sparse data.

Actual SB9, full Traven 2020, VSX, WDS, ASAS-SN, Gaia alias/NSS and SIMBAD checks
are retained with positive controls. An erroneous ASAS-SN table name originally
returned a service error; it was rejected, corrected to II/366/catv2021 and
validated on an actual positive source. A nearby known SB2 is 111 arcsec from
the second target and is not assigned to it. No first-companion claim follows.
Read [the spectroscopy report](reports/new_stars/spectroscopy.md).

## Preservation and the next scientifically useful steps

Raw images, spectra, QLP data and all 218 intermediate image plots remain in
`/workspace/tess-research-data/new_stars/`. The 186 MB bulk plots were moved intact
outside ordinary Git and verified before/after; four useful image diagnostics
remain tracked. Original URLs, processing identifiers, dates, failed responses,
raw hashes, model outputs and method snapshots are preserved. Small public-query
responses are additionally archived for exact offline restoration. External
storage durability across new cloud tasks remains unverified.

The final checkpoint verifies 10 method files, 479 scientific artifacts and
5,955 cached files (377.9 MB). An offline restoration test recovered all 5,466
archived query files with identical sizes and SHA-256 hashes. The archive is
1.36 MB compressed; larger observations retain retrieval URLs and checksums.

Run `python scripts/new_star_checkpoint.py verify`; offline query restoration
is `python scripts/new_star_checkpoint.py restore-queries`. Missing static
observations can be retrieved with `download --max-downloads N`. Byte identity
does not establish astrophysical validity or publication priority. Public
archive/worker failures are documented rather than treated as empty sky. The
parallel workers hit a usage limit; the primary agent continued local spectral,
photometric, mask, catalogue, SDSS and preservation checks. No paid APIs,
external messages, submissions or observing requests were used.

The next work should prioritize differential RVs and expanded source-specific
literature checks for pairs 1/3; genuinely reduced spectra or additional epochs
for the RV leads; and independent-night, higher-resolution images for the faint
unidentified peak. A new blind-source phase should exploit higher angular
resolution or motion-aligned multi-epoch images rather than equating Gaia
absence with a new star. The completed sample does not exhaust the archives,
and the investigation remains open.
