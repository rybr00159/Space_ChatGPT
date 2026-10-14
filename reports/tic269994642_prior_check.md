# TIC 269994642: recognition scope and source context

Checked **2026-10-14 UTC**. This is a recurring, spatially localized
**transit/eclipse signal of unresolved origin**, not an accepted new planet.
No previous planet or variable identification was found in the finite checks
below. Gaia astrometry and the grazing shape leave stellar-binary alternatives
open; the absence of a catalog entry does not establish novelty.

## Current host measurements and validity

Live [ExoFOP](https://exofop.ipac.caltech.edu/tess/target.php?id=269994642)
identifies TIC 269994642 with Gaia DR2/DR3 **6806240691223695744** and aliases
2MASS J20455077−2438595, APASS 12371401, UCAC4 327−194685, and
WISE J204550.76−243859.9. Its response cites TIC v8.2 parameters dated
2019-04-15: **R=0.835257±0.0872364 solar radii**, **M=0.631±0.0830436 solar
masses**, **Teff=4035±126.582 K**, logg 4.39446±0.117306, and distance
202.6±2.0125 pc. The current TESS magnitude is 12.9234 and the modeled TIC
contamination ratio is 0.003582. That ratio is a catalog estimate, not a measured
aperture fraction or exclusion of an unresolved companion. The old sector-1 QLP
header radius 0.629247 should not be substituted for the refreshed radius.

Gaia DR3 provides a five-parameter solution with G=13.745143,
BP−RP=1.73653 and parallax **5.01109±0.04304 mas**. Its proper motion is
−14.1168/−46.7558 mas/yr. **RUWE=2.081322**, astrometric excess noise
0.218268 mas and significance 35.5775 indicate an imperfect single-source
astrometric fit. These warrant a multiplicity/source-quality warning, but they
neither prove a binary nor identify the source of the eclipse. IPD multiple-peak
fraction is 0%, IPD harmonic amplitude 0.02685 and BP/RP excess factor 1.3050.

`non_single_star=0` means no published NSS processing flag here. It does not
exclude a binary. `phot_variable_flag=NOT_AVAILABLE` and empty Gaia
variability-classifier, rotation-modulation and eclipsing-binary tables likewise
do not establish constant flux. The retrieved astrophysical-parameter row has
**no GSP-Phot/FLAME temperature, radius or mass**, so Gaia does not supply an
independent stellar-radius estimate in this check.

Gaia median RV is **−39.9797±2.0989 km/s**, based on 26 transits. Robust RV
amplitude, variability p-value and renormalized goodness of fit are missing.
The median is not an orbital amplitude or mass constraint, and its existence
cannot rule out a stellar companion.

## Spatial neighborhood and contamination alternatives

The preserved Gaia 120-arcsec cone contains **44 rows**. The expanded MAST
TIC cone is `COMPLETE` and contains **45 rows**. The short ExoFOP nearby list
does not include the entire 120-arcsec neighborhood, so the expanded cone is
used for the following context. Separations below use the Gaia target position;
they can differ slightly from TIC/ExoFOP positions because the reference epochs
and proper-motion handling differ.

| Neighbor | Separation (arcsec) | TESS mag | Full TESS flux / target | Context |
|---|---:|---:|---:|---|
| TIC 269994643 / Gaia 6806240691224982016 | 20.85 | 18.7516 | 0.00467 | DR2 and DR3 source IDs differ; faint closest ordinary neighbor |
| TIC 269994635 / Gaia 6806287729705524352 | 50.91 | 14.9032 | 0.1615 | Gaia RUWE 1.055, current TIC radius 1.3711 |
| TIC 269994641 / Gaia 6806287592266573312 | 71.26 | 14.3331 | 0.2730 | Gaia RUWE 0.952, current TIC radius 1.1751 |
| TIC 269994646 / Gaia 6806287587969025920 | 80.85 | 15.5202 | 0.0915 | Gaia RUWE 1.618 and IPD multiple-peak fraction 90% |
| TIC 269994651 / Gaia 6806240656863954944 | 81.86 | 14.8847 | 0.1642 | Gaia RUWE 0.943, current TIC radius 1.0470 |
| TIC 269994655 / Gaia 6806240278906832768 | 117.63 | 14.3533 | 0.2680 | Gaia RUWE 0.969, current TIC radius 1.0547 |

These are full catalog flux ratios, **not measured aperture contributions**.
The closest faint source alone has less catalog TESS flux than the observed
approximately 2–3% event; the more distant sources have enough full flux to
warrant spatial tests. Flux ratios alone cannot establish or exclude actual
blend contamination, given passband, aperture, background and variability
uncertainties.

The independent retained pixel analysis in
`candidates/broad_s1/tic269994642_pixels/sector1/` and `sector92/` recovers
the narrow signal at the target neighborhood in distinct sectors/cameras.
Sector 92 provides 200-second rather than sector-1 1800-second sampling.
The source-localization uncertainties and differential-image modeling must
be consulted before assigning a subpixel offset to a companion. This report
does not independently reprocess the pixel data; it supplies the catalog
counterhypotheses for that analysis.

SIMBAD returns only **1RXS J204550.6−243846**, an X-ray source roughly
14 arcsec from the target, within the 120-arcsec cone. It has no linked
references in the saved SIMBAD response. ROSAT angular resolution and
the nearby optical sources prevent a unique association here. It is possible
activity context, not proof that the target is active, rotating or binary.

## Finite prior-recognition checks

The saved `tic269994642_cached_matches.json` checks exact TIC identity and
positions within 120 arcsec against the dated project snapshots of NASA
Exoplanet Archive confirmed planets (including RV-only hosts), TOIs, ExoFOP
CTOIs, the actual Prša TESS EB membership table, DIAmante v2 candidates and
published Olmschenk/Yu candidate lists. It finds **no host, nearby, exact-period
or harmonic-period matches** at the tested period 1.74649452 days. The matching
scope, tolerance, original snapshot dates and hashes are embedded in that file.

Additional live checks on 2026-10-14 are preserved as public response bodies:

- Target ExoFOP has empty confirmed-planet, TOI, CTOI and planetary-parameter
  fields, and no uploaded files, imaging, spectroscopy or time-series entries.
  ExoFOP pages for bright neighbors TIC 269994635, 269994641, 269994651 and
  269994655 also have empty planet/TOI/CTOI fields and no files.
- The live Villanova TESS EB endpoint does not identify TIC 269994642 as a
  catalog member; this is consistent with the actual EB membership snapshot.
- VizieR **VSX version 2026-10-05** and **ASAS-SN II/366** within 120 arcsec
  return no entries. The ASAS-SN mirror cites a 2024 database retrieval and is
  not an exhaustive contemporary variability search.
- SIMBAD's 120-arcsec cone supplies only the ROSAT source above, with zero
  linked references; its bibliography query is empty.
- An arXiv API query for the exact TIC string and all listed ExoFOP aliases
  returns **zero results**. This tests indexed alias strings, not every paper,
  table, supplementary file, survey detection or differently named object.
- The exact-target MAST CAOM query (`obs_collection=TESS`, bare
  `target_name=269994642`) returned `COMPLETE` with zero observation rows.
  Its public response and request are retained in
  `tic269994642_mast_target_observations.json` and its UTC/hash manifest.
  This limited indexing query cannot exclude public FFI-based TCEs, reports,
  or products indexed under different names; it is not an exhaustive TCE search.

These checks do not exhaust ADS, all public SPOC TCE reports, survey candidate
databases or literature tables. A surviving planetary interpretation would
require a broader literature/TCE and alias investigation before a novelty
claim. Catalog absence here is recorded as **no recognition found in the
checked scope**, not “previously undiscovered.”

## Competing photometric interpretations

The sector-1 local raw-SAP review recovered P≈1.74649452 days, a fitted width
about 1.73 hours and raw local S/N≈18.7. Subsequent independently extracted
sector-92 photometry shows recurrence near the same photometric period.
These are strong evidence for a repeatable eclipse/transit feature; they are
not a calibrated search-wide false-alarm rate or a planet classification.

The retained joint physical analysis in
`candidates/broad_s1/tic269994642_physical/` explicitly tests grazing, dilution,
density and binary alternatives. Its single-star dark-companion fit is
grazing and permits a broad companion-radius profile extending into stellar
scales. An **equal-light grazing binary at twice the photometric period**
(about 3.493 days) fits comparably: equal alternating eclipses evade both a
half-phase secondary test at the shorter period and an odd/even test. Equal
unresolved stellar flux would also make the catalog combined-light radius
different from each component radius. An unequal-brightness binary improves
the adopted diagonal residual statistic in that analysis. These fits are
conditional illustrations, not calibrated model odds or a binary diagnosis;
correlated errors, baseline choices and stellar assumptions remain material.

The retained event-level parity review in `parity_robustness.json` separately
fits fixed-shape amplitudes and constants to 24 events. The even/odd amplitude
difference is 1.28 conditional standard errors in sector 1 and 3.09 in sector
92, in the same direction. A conservative pooled sector comparison gives
3.33 conditional standard errors; 3,000 within-sector event bootstrap draws
and event-scatter/error inflation were also evaluated. This **targeted,
same-photon diagnostic favors a stellar-binary interpretation** but is not a
search-wide false-alarm probability. Fixed-shape mismatch, activity and baseline
systematics remain alternatives, and two sectors cannot estimate unobserved
long-term systematics. These diagnostics must not be multiplied as independent
probabilities.

The useful next tests are parity consistency in additional independent data
at the doubled period, all-phase variability on unnormalized photometry,
high-resolution imaging/stellar characterization, and phase-resolved RV
observations. The current evidence supports retaining this as a strong
**unclassified eclipsing source**, with the stellar-binary explanation favored
by conditional alternating-depth evidence and planetary alternatives not
excluded. It does not support accepting a new exoplanet.

## Preservation and reproduction

All source bodies and request URL/ADQL/UTC/hash manifests are retained as
`data/broad_catalogs/tic269994642_*`. The initial combined Gaia query timed
out; the successful split-query cone has HTTP 200, 44 parsed rows and a
matching SHA-256. Its adjacent `gaia_nearby.csv.error.json` records a **later
display-parser failure on an empty G value**, not failure of the saved cone.
Missing fields are preserved, not converted to zero.

`python scripts/tic269994642_source_context.py` reconstructs the derived
neighborhood CSV and verifies every retained response hash without changing
original response bodies. The verified inventory is
`reports/tic269994642_prior_check_sources.json`. Original observations,
native FFI identifiers, extraction diagnostics and conditional fit products
remain in their separate candidate folders. No paid model/API calls were used.
