# Current stellar and finite recognition checks for four secondary QLP leads

Checked 2026-10-14 UTC. **None of these four leads currently supports a credible
new planetary interpretation.** Three target stars are evolved giants, and
the fourth has large Gaia radial-velocity variability plus a bright 2.47-arcsec
companion. These are falsification and prioritization results; they do not
identify the exact origin of every dip or establish confirmed binaries.

## Summary

| TIC | Refined period (days) | Current TIC radius (solar) | Gaia radius (solar) | Principal finding |
|---|---:|---:|---:|---|
| 79354493 | 6.521071 | 9.78722 | GSP-Phot 9.4753; FLAME 9.943618 | Giant, independently associated with published oscillating-giant work |
| 231825824 | 6.246274 | 9.71328 | GSP-Phot 11.1553; FLAME 10.0826235 | Giant, independently associated with published oscillating-giant work |
| 259901286 | 3.789038 | 9.57773 | No stellar solution in queried Gaia fields | TIC giant estimate and published oscillating-giant association; header radius missing |
| 308600649 | 1.625939 | 2.18313 ± 0.100543 | GSP-Phot 2.1927; FLAME 2.181306 | Gaia robust RV amplitude 41.743 km/s; unresolved TESS neighbor 2.47 arcsec away |

The sector-1 QLP header radii are respectively 0.78682, 0.825081, -1,
and 1.60697 solar radii. Those early header values must not be used as current
host-radius evidence. ExoFOP is a current response, but its TIC stellar rows
explicitly cite TIC v8.2, dated 2019-04-15; this is not a newly measured spectrum.

Using the existing local-quadratic KSPSAP depths of 0.00119898, 0.00120324,
0.00128338, and 0.00418405, the illustrative expression
`Rstar * sqrt(depth)` gives 0.3389, 0.3369, 0.3431, and 0.1412 solar radii
(about 37.0, 36.8, 37.4, and 15.4 Earth radii). This is a depth-size proxy,
not a transit fit or radius posterior. It assumes an event on the target,
negligible dilution, and suitable nongrazing geometry. Its first three values
are incompatible with an ordinary planet-sized occultor on the catalog giant;
oscillations/granulation, detrending artifacts, binary eclipses, and a blended
source must be tested. A background planet around a different star is not
excluded by the giant target radius, but that interpretation needs independent
source localization. The fourth proxy is potentially giant-planet sized,
but RV and contamination warnings substantially weaken that interpretation.

## Target details

### TIC 79354493

Current TIC temperature is 4963 ± 122 K and radius 9.78722 solar radii; Gaia
GSP-Phot gives 4921.3535 K, log g 2.6633, radius 9.4753 (reported bounds
9.3343–9.6220), and FLAME mass 2.8014076 solar masses. Gaia DR3 source is
6466365360560070912. RUWE is 1.171754, IPD multi-peak fraction 1%, parallax
1.125026 ± 0.020441 mas, `non_single_star=0`, and
`duplicated_source=false`. Astrometric excess-noise significance is 18.70;
RUWE and the other catalog fields do not exclude an unresolved companion.

Median Gaia RV is 8.315189 ± 0.236914 km/s from 22 transits, robust amplitude
3.514428 km/s, RV chi-square probability 0.00163408, and renormalized RV
goodness of fit 2.02570. This suggests catalog RV variability but does not
provide an orbit or planetary/binary mass. The small-amplitude variability
could include stellar variability and catalog/systematic effects.

The 120-arcsec Gaia cone has 44 rows including the target. The nearest
catalog neighbor is at 11.15 arcsec, G 19.7951 versus target G 10.2914;
its full G-band flux ratio is 0.000158. No neighbor has G-band flux ratio
above 5% in this finite cone. This is not a pixel-level contamination veto.

SIMBAD associates CD-52 9823 within 0.0676 arcsec, type `*`, with five
references. One is Hon et al. 2021, ApJ 919, 131, DOI
10.3847/1538-4357/ac14b1, the QLP oscillating-red-giant catalog. The other
references include RAVE and APASS catalog papers. The association strongly
supports prior giant-star analysis but does not establish that the 6.52-day
screened dips themselves were previously characterized.

### TIC 231825824

Gaia DR3 source 4712315388482342528 differs from its recorded DR2 source
4712315384186079488; the ExoFOP DR3 binding was used explicitly. Current TIC
temperature is 4944.65 ± 45.96 K and radius 9.71328 solar radii; Gaia GSP-Phot
gives 4919.3574 K, log g 2.5778, radius 11.1553 (reported bounds
10.9614–11.9237), and FLAME mass 2.801597 solar masses. Radius estimates
differ at a scientifically relevant level, but all place this star in the
evolved-giant regime rather than the 0.825-solar-radius old header regime.

RUWE is 0.9011672, IPD multiple-peak fraction 1%, parallax
0.866397 ± 0.010935 mas, `non_single_star=0`, and
`duplicated_source=false`. Its astrometric solution is six-parameter
(`astrometric_params_solved=95`). Median Gaia RV is 14.498685 ± 0.248335
km/s over 24 transits, robust amplitude 4.067407 km/s, chi-square probability
0.00269767, and RV renormalized goodness of fit 2.31865. These metrics
indicate variability with the same limitations described for TIC 79354493.

There are 16 Gaia rows in the 120-arcsec cone including the target. Neighbor
Gaia 4712315384186079360 = TIC 231825825 is 38.42 arcsec away, G 12.3484,
full G-flux ratio 0.2626; ExoFOP gives TESS magnitude 11.9252 versus
target 10.3004, full TESS-flux ratio about 0.224. Another G 12.7202 neighbor
is 115.67 arcsec away with G-flux ratio 0.1865. These are aperture-dependent
blend warnings; source attribution remains required.

SIMBAD associates TYC 8852-667-1 within 0.0980 arcsec, type `*`, with five
references including Hon et al. 2021 and RAVE/APASS. The stellar context is
already published, although the specific screened period was not found in
the finite candidate catalogs below.

### TIC 259901286

Gaia DR3 source is 4646596960014607488. The ExoFOP TIC v8.2 radius is
9.57773 solar radii and temperature 4821 ± 122 K. All queried Gaia
GSP-Phot/FLAME stellar values are empty for this source: Gaia did not
independently deliver a radius in this query. Its parallax is
0.630385 ± 0.010755 mas, RUWE 0.9519212, IPD multi-peak fraction 0%,
`non_single_star=0`, and `duplicated_source=false`.
The median RV is -9.219384 ± 0.375494 km/s over 21 transits, robust amplitude
4.313019 km/s, chi-square probability 0.185586, and RV renormalized goodness
of fit 0.650693. Unlike the first two, these catalog variability statistics
do not strongly reject constant RV.

The 120-arcsec Gaia cone has 30 rows including the target; G 10.9738 source
4646608814124343040 at 108.34 arcsec has 1.856 times the target's full
G-band flux (target G 11.6452). Its coupling into the actual aperture has not
been measured. SIMBAD associates TYC 9148-821-1 within 0.0493 arcsec, with
one reference, Hon et al. 2021. This corroborates the evolved/oscillating
star context despite the missing Gaia stellar solution. A direct extraction
of the Hon supplementary TIC row would be stronger than bibliography alone.

The two coarse peaks (1.958448 and 3.795368 days) are preserved separately
in the cached signal matches; the photometric refined lead is 3.789038 days.
These are not interchangeable orbital ephemerides.

### TIC 308600649

Gaia DR3 source 5277483675334676480 has RUWE 0.88822484, IPD multi-peak
fraction 0%, parallax 1.697192 ± 0.012698 mas,
`non_single_star=0`, and `duplicated_source=false`. GSP-Phot temperature
is 6476.50 K, log g 3.9344, radius 2.1927 (reported bounds 2.0857–2.2529),
and FLAME mass 1.5170454 solar masses, broadly consistent with the TIC
subgiant radius.

**Its Gaia robust RV amplitude is 41.742954 km/s across 16 transits**, with
RV chi-square probability recorded as 0.0 (finite catalog rounding, not
mathematical probability zero) and renormalized goodness of fit 10.970441.
Median RV is 19.02023 ± 3.91555 km/s. This is a strong nonplanet warning:
spectroscopic binarity, variable/composite spectra, or blended-source RV
bias require assessment. Gaia has no orbit/phase linkage to the 1.625939-day
photometric signal here, so no companion mass is inferred from the amplitude.

The 120-arcsec cone has 92 Gaia rows. **TIC 806435656 / Gaia
5277483671040252928 lies only 2.4733 arcsec away**, much smaller than a
TESS pixel. Its G is 13.82823, full G-flux ratio 0.12570, RUWE 1.46663.
ExoFOP nearby sources give target TESS magnitude 11.1906 and this neighbor
13.3655: full TESS-band flux ratio about 0.1349. If both stars entered the
aperture fully, an eclipse on the neighbor at only roughly 3.5% intrinsic
depth could produce a combined 0.418% dimming. This is an illustrative
blend feasibility estimate, not an aperture model or observed neighbor eclipse.
Its TIC distance 590.106 ± 7.391 pc and target 588.496 ± 8.769 pc and
similar proper motions also warrant assessment as a possible physical companion.
Distance agreement alone does not prove binding.
An additional Gaia query measures parallax 1.688591 ± 0.019250 mas for the
neighbor versus 1.697192 ± 0.012698 mas for the target; its median RV is
18.886608 ± 4.124171 km/s over 14 transits, without a published robust RV
amplitude or RV variability statistic. Similar catalog median RVs do not
resolve the target's large RV variability or establish a bound orbit.
Current ExoFOP estimates neighbor temperature 6137 ± 151 K, radius
0.889741 solar radii, and mass 1.16 solar masses, without radius/mass errors.
Gaia supplies no GSP-Phot/FLAME stellar values for the neighbor in the queried
fields. Its ExoFOP TOI/CTOI/confirmed-planet lists are empty and the live
Villanova page has no associated EB entries. These nulls do not exclude an
eclipsing companion at the candidate ephemeris.

Current TIC contamination ratio is 0.176996, a model prediction for an
assumed aperture, not the measured contamination of this QLP light curve.
Several other relatively bright Gaia neighbors are 78–116 arcsec away.
A planet interpretation requires deblended spatial photometry and a
velocity/source model. Normal target RUWE does not negate either warning.

SIMBAD associates TYC 8928-2416-1 within 0.1195 arcsec, type `*`, with
zero references. Live Villanova returns **“No associated EB entries found.”**
The same detail includes HLSP product links, which are not EB membership.

## Finite prior-recognition scope

- Cached `CatalogIndex(position_radius_arcsec=120)` checks use the existing
  frozen NASA confirmed-planet table for all discovery methods (including RV),
  TOIs, ExoFOP CTOIs, DIAmante v2, Prsa2022 TESS EBs, Olmschenk2021,
  and Yu2019. Exact TIC and <=2-arcsec host associations, nearby positions,
  and exact/2–5 integer harmonic/subharmonic matches within 1% of both
  coarse periods were evaluated. **No cached catalog matches** were found
  for any of the four targets within 120 arcsec. Source snapshot dates,
  exact hashes, missing-optional-source state, and per-period results are
  retained in `tic{ID}_cached_matches_current.json`; “current” here is the
  computation date, not a claim that all underlying catalog snapshots are live.
- Current ExoFOP target responses have empty TOI, CTOI, confirmed-planet and
  planet-parameter lists for all four targets. Empty lists do not exclude
  unpublished or recent reports.
- Live Villanova explicitly states “This TESS ID is not in the catalog” for
  the first three; TIC 308600649 states “No associated EB entries found.”
  Membership is not inferred from a product link.
- Gaia variable classifier and eclipsing-binary tables contain no target
  entries for all four. Target `phot_variable_flag=NOT_AVAILABLE` indicates
  absent published Gaia classification, not photometric constancy.
- VizieR's VSX mirror, explicitly **Version 2026-10-05**, returns no
  variable-star rows in each 120-arcsec cone. This is a finite mirrored
  catalog check rather than a guarantee about the continuously current AAVSO
  database.
- Finite arXiv metadata queries using the TIC identifier and available
  Gaia DR3, Tycho, and 2MASS tokens return zero indexed-metadata results for
  each target. They do not search full paper text or every supplementary
  table; the SIMBAD Hon references demonstrate why metadata nulls cannot
  establish novelty.

The exact previous recognition of any screened dips remains incomplete.
Outstanding checks include full TCE/DV tables, recent non-arXiv papers and
supplements, Planet Hunters, all variable/binary surveys, and all alias variants.
No discovery or confirmation claim is supported by these nulls.

## Preserved evidence and reproduction

Public bodies, exact URLs/ADQL, actual UTC retrieval times, sizes, HTTP status,
and SHA-256 are under `data/broad_catalogs/tic{ID}_*` with
`.manifest.json` sidecars. Relevant suffixes are
`exofop_current.json`, `exofop_nearby_current.csv`,
`gaia_target_current.csv`, `gaia_stellar_current.csv`,
`gaia_neighbors_current.csv`, `gaia_neighbors_derived_current.csv`,
`gaia_class_current.csv`, `gaia_eb_current.csv`,
`simbad_current.csv`, `simbad_current_aliases.csv`,
`simbad_current_references.csv`, `tesseb_live.html`,
`vsx_current.tsv`, and `arxiv_current.xml`. Current response files were
newly added; no global catalog source or earlier scientific response was replaced.
Gaia G flux ratios in derived tables are explicitly band-labeled and are
not TESS aperture fractions.
Additional close-neighbor responses for TIC 308600649 use suffixes
`gaia_close_neighbor_current.csv`,
`close_neighbor806435656_exofop_current.json`, and
`close_neighbor806435656_tesseb_current.html`, with the same manifests.

Run `/workspace/openrouter-check-venv/bin/python scripts/secondary_leads_catalog_check.py`
to produce new stellar/Gaia/cached recognition snapshots for these four IDs.
Existing names receive a UTC timestamp; the script refuses replacement of a
name collision. Other read-only checks can be replayed with the exact public
URLs recorded in their manifests, with TLS verification enabled, and saved
as new snapshots. The archived database used to extract RA/Dec and periods
was opened with SQLite `mode=ro`; no broad-prioritization tables were modified.

Photometry comes from the existing sector-1 candidate `vetting.json` files.
These qualitative contextual results must be integrated into ranking and
subsequent physical validation; statistical false-alarm probabilities still
require correlated-noise and search-trials calibration. No paid requests,
credentials, or API key inspection were used.

## Next scientific actions

Prioritize independent sectors and source localization for other compact-star
leads before spending more time on these three giant-host interpretations.
For giant signals, directly join Hon et al.'s public supplement and compare
oscillation/granulation timescales with the screened event morphology; do not
blindly treat a stellar variability catalog as the exact cause of every dip.
For TIC 308600649, first obtain spatially resolved source photometry and
assess the close neighbor and Gaia RV data, then test eclipse/ellipsoidal/
binary models alongside any planetary model. This is a high-value
false-positive study rather than a credible new-planet claim at this point.
