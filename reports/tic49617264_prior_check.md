# TIC 49617264: current stellar context, prior recognition, and nearby sources

Checked 2026-10-14 UTC. **This is an already publicly followed planetary
candidate, not a newly unreported signal.** It is a public SPOC TCE with a
4.8319416046-day period in the current investigation, and public ground
follow-up from 2024 explicitly identifies it as a neighbor planetary
candidate to TOI-2417.01. Empty TOI fields on this particular TIC do not
establish novelty.

## Positive prior recognition: public ground follow-up

The public ExoFOP target JSON lists three completed ground light-curve
observations (2024-09-15 TRAPPIST-South in zp; 2024-10-05 Brierfield in Rc;
2024-10-15 El Sauce in B), plus 35 associated files. The Brierfield public
notes file, retrieved without authentication at
https://exofop.ipac.caltech.edu/tess/get_file.php?id=1004492, states:

> “TIC 49617264.11 (TOI 2417.11)”

> “[NPC to original candidate TIC 49617263 (TOI-2417.01)] SAFFIs constrain
> event to this star relative to Gaia DR2 stars bright enough to cause the
> TESS detection.”

It reports a full, on-time, on-target 13.1-ppt Rc event in a 6.6-arcsec aperture
on 2024-10-05, and quotes a 12.8-ppt zp event measured by TRAPPIST on
2024-09-15 in a 4.5-arcsec aperture. The reported revised period is
4.83261 days, compared with the original 4.8326600-day period. It quotes an
approximately 2h46m duration and describes VPC comments and follow-up
requests to test chromaticity. It also reports 4/6 nearby sources cleared
in the local NEB check, with T2/T3 not cleared because they contaminate each
other.

These are **existing observers' reported conclusions**, not our independent
reanalysis or planetary confirmation. The term “TOI 2417.11” in this note
must not be promoted to official current TOI table membership on TIC49617264;
the current ExoFOP TOI and CTOI lists for this TIC are empty. Nonetheless the
positive follow-up history conclusively prevents interpreting a blank host
TOI row as evidence of a previously unrecognized photometric candidate.

The full raw note and its exact URL, actual UTC retrieval, size, HTTP200,
and SHA-256 are preserved in
`data/broad_catalogs/tic49617264_followup20241005_notes.txt` and its manifest.
The note has small date/telescope transcription inconsistencies; the dated
ExoFOP observation metadata provides the observing log. This research has
not independently verified its fitted parameters or source-clearing claims.

## Nearby matching catalog signal and contamination

The cached 120-arcsec all-method/position/period crossmatch finds one nearby
TOI entry, **TIC49617263 / TOI-2417.01**, period **4.83266 days**,
duration 2.308 hours, epoch BTJD2139.2966, **FP disposition**. The relative
period difference versus 4.8319416046 is 0.000148654, or **0.0148654%**.
This is an exact-frequency nearby match, not a new extra-planet period or a
harmonic. Its cached separation is 43.18 arcsec using the stored catalog
coordinates; current Gaia/TIC positions give roughly 43.42 arcsec. These
differences do not represent distinct neighboring stars.

The live ExoFOP response for TIC49617263 likewise identifies TOI2417.01 with
TFOP disposition FP, TESS disposition PC, and QLP sectors3,30. Its note says
“Depth-aperture. Potential visual binary 2arcsec away in Gaia with similar
distance and proper motion but 3 mags fainter.” This is consistent with a
crowded-source origin requiring careful source assignment.

Initially a matching bright neighbor could be suspected of contaminating
the target. **The independent prior ground follow-up above specifically
reports the opposite source assignment: TIC49617264 is the transit source
that contaminated the original brighter TOI target.** Thus the correct
working interpretation is a previously recognized neighbor candidate,
pending independent pixel/source checks, rather than assuming that the
brighter neighbor itself eclipses. An FP disposition on the original TOI
does not, by itself, establish that the relocated signal is nonplanetary.

The 120-arcsec Gaia cone contains 11 sources including the target:

| Source | Separation arcsec | Gaia G | Full G-flux ratio to target | Relevant TIC/TESS context |
|---|---:|---:|---:|---|
| TIC49617264 / Gaia4952839776170905472 | target | 13.51845 | 1 | T12.9915 |
| TIC49617263 / Gaia4952839741811566720 | 43.42 | 10.70592 | 13.3357 | TOI2417.01 FP; T10.2589, about12.39 times target full TESS flux |
| TIC631878531 / Gaia4952839741810218880 | 43.96 | 13.79513 | 0.7750 | T13.4007; RUWE1.57498, near the bright TOI |
| Gaia4952839879250119552 | 70.70 | 14.36117 | 0.4602 | Finite-cone contaminant to consider |
| Gaia4952838951537184896 | 116.44 | 15.70800 | 0.1331 | Distant aperture-wing possibility |

Full catalog flux ratios are not aperture fractions. Gaia G ratios are not
TESS-passband ratios. The current TIC contamination ratio **0.316295** is a
model estimate for an assumed aperture, not the measured dilution of the
specific SPOC/QLP aperture. Neighbor pixel attribution and time-series
photometry are required; the bright TOI's extremely close secondary merits
an eclipse/blend test alongside the already-reported target transit.

SIMBAD has no object within2arcsec of TIC49617264. It associates the bright
neighbor with TYC7555-691-1, aliases TIC49617263 and TOI-2417, at43.61arcsec.
Its bibliography includes the SOAR TESS companion survey,
2021AJ....162..192Z, DOI10.3847/1538-3881/ac17f6, as well as RAVE/APASS and
red-clump-calibration work. This bibliography is about the neighboring star,
not automatically the target. A catalog radio source NVSSJ023134-400013 is
6.128arcsec from the target and is referenced in SPECFINDV3.0; it is not
a sufficiently secure target association to infer stellar radio activity.

## Current stellar and Gaia target quality

Current TIC v8.2 parameters through ExoFOP are temperature5201±127.96K,
radius1.03812±0.0606231 solar radii, mass0.88±0.114027 solar masses,
logg4.35005±0.08610, and distance472.188±3.345pc. The row is dated
2019-04-15; a current retrieval does not make this a new spectroscopic fit.

Gaia DR3 source4952839776170905472 gives:

- GSP-Phot temperature5194.238K (reported lower/upper5183.026–5211.444).
- GSP-Phot logg4.4432; radius0.9941 solar radii (bounds0.9748–1.0026).
- FLAME radius1.0069487 solar radii; no supplied FLAME mass.
- Parallax2.108891±0.011401mas; RUWE1.0264244; IPD multi-peak fraction0%.
- Five-parameter astrometry (`astrometric_params_solved=31`), excess
  astrometric noise0.0, duplicated_sourcefalse, non_single_star0.
- Median RV23.265215±1.347732km/s across34transits, but **no robust RV
  amplitude, RV chi-square probability, or renormalized RV fit statistic**.
  This median cannot exclude a binary or yield a planetary mass.

The Gaia radius estimates support a roughly solar-radius host, with no
severe target astrometric warning in these fields. That is useful contextual
evidence, not a planetary validation. Gaia variable classifier and
eclipsing-binary tables have no target row; phot_variable_flagNOT_AVAILABLE
means missing published classification, not constant flux.

## Finite other recognition checks and limitations

The cached catalog match includes frozen NASA confirmed planets of all
discovery methods, TOIs, CTOIs, DIAmantev2, Prsa2022 EBs, Olmschenk2021,
and Yu2019, plus positional and 2–5 integer harmonic/subharmonic checks.
There is no exact host match; the nearby TOI is the decisive positive
match. Dates and hashes of every underlying catalog are preserved in
`tic49617264_cached_matches_current.json`.

Current target ExoFOP lists empty TOIs, CTOIs, confirmed planets and
planet_parameters. **These nulls are superseded for novelty by its positive
uploaded follow-up record.** Target Villanova explicitly states “This TESS
ID is not in the catalog”; bright-neighbor Villanova has “No associated EB
entries found,” despite photometry links. VSX's VizieR mirror explicitly
Version2026-10-05 has no variable-star row within120arcsec. Finite arXiv
metadata queries for target and bright-neighbor identifier/alias tokens
return zero entries. These finite nulls neither establish novelty nor
exclude a binary; they are preserved as scope documentation.

The reported transit's planetary nature remains unresolved without
independent transit fitting, chromaticity, deblending, blend/EB modeling,
ground-photometry reanalysis and preferably spectroscopy. This target may
still be scientifically useful for candidate validation, but **it cannot
serve as an undiscovered-signal success for this research objective**.

## Preservation and reproduction

All public response bodies and SHA-256+URL+UTC manifests use
`data/broad_catalogs/tic49617264_*`; earlier response files and global source
catalogs were not overwritten. The Gaia neighbor-derived CSV records band
labels and hashed sources. Primary suffixes are
`exofop_current.json`, `exofop_nearby_current.csv`,
`nearby49617263_exofop_current.json`, `followup20241005_notes.txt`,
`gaia_target_current.csv`, `gaia_stellar_current.csv`,
`gaia_neighbors_current.csv`, `gaia_neighbors_derived_current.csv`,
`gaia_class_current.csv`, `gaia_eb_current.csv`, `tesseb_live.html`,
`tesseb_live_brightneighbor.html`, `simbad_current.csv`,
`simbad_current_brightneighbor_{aliases,references}.csv`,
`simbad_current_radio_{aliases,references}.csv`, `vsx_current.tsv`,
`arxiv_current.xml`, `arxiv_current_brightneighbor.xml`, and
`target_page_current.html`.

From the repository root run
`PYTHONPATH=. /workspace/openrouter-check-venv/bin/python scripts/tic49617264_catalog_check.py`
to save new ExoFOP/Gaia/cached catalog snapshots. Replay other checks using
the exact public URLs recorded in their manifests, with TLS verification
enabled and a fresh dated output. Do not substitute a failed/changed public
response for the preserved snapshot. No paid APIs, credentials, or API-key
inspection were used. Independent pixel vetting of sectors106/107 is
performed elsewhere in this research; this report does not claim its results.
