# TIC 214311908: independent catalogue and source audit

Checked 2026-10-17 using public services, independently of the other model's
findings. These checks identify the target and potential blends; they do not
establish that any measured TESS dip belongs to the target or validate a planet.

## Target

TIC 214311908 is **Gaia DR3 5925373272850784768 / 2MASS
J17235292-5119271**, a catalogued nearby cool dwarf. ExoFOP coordinates are
RA 260.96969903043°, Dec −51.32383313824°; Gaia DR3 positions refer to J2016.0
and differ slightly. Its substantial proper motion must be propagated when
comparing different epochs.

| Quantity | Public catalogue value |
|---|---:|
| TIC TESS magnitude | 13.0135 ± 0.0075 |
| TIC effective temperature | 3301 ± 157 K |
| TIC radius | 0.380658 ± 0.011660 solar radii |
| TIC mass | 0.369492 ± 0.020483 solar masses |
| TIC density | 9.445379 ± 0.344157 g cm⁻³ |
| TIC distance | 64.1712 ± 0.3479 pc |
| Gaia DR3 G magnitude | 14.281043 |
| Gaia DR3 parallax | 15.497306 ± 0.025122 mas |
| Gaia DR3 proper motion (RA cos Dec, Dec) | (−113.6712 ± 0.0311, +85.5436 ± 0.0224) mas yr⁻¹ |
| Gaia DR3 RUWE | 1.07620 |
| Gaia DR3 radial velocity | 16.2837 ± 2.2174 km s⁻¹; 21 transits |

The stellar estimates are catalogue inputs, not new measurements. Gaia DR3
astrophysical-parameter fields for this target are empty. Its Gaia variability
flag is `NOT_AVAILABLE`, and exact source-ID queries return no Gaia eclipsing
binary, rotation, classifier, astrometric-orbit, acceleration, or nonlinear
spectroscopic solution. These nonentries do not exclude binarity or variability;
the catalogue supplies neither an RV semi-amplitude nor a planet mass.

## Prior recognition

- Fresh ExoFOP target JSON: no TOI, CTOI, confirmed planet, planet parameter,
  uploaded file, or photometric follow-up entry.
- Fresh NASA confirmed-planet query: no exact TIC or source within 120″.
- Fresh NASA TOI query: no exact TIC or source in a conservative bounding box
  enclosing the 120″ cone. An initial spherical query failed with the archive's
  `ORA-00904 HTM20` error; that failure and successful replacement are retained.
- Fresh complete ExoFOP CTOI download: 5,138 rows, no exact TIC entry.
- RAVEN Zenodo release 19661443: no exact TIC in any of its four overlapping
  tables (14,815 / 2,170 / 143 / 207 rows). All downloaded files match the
  release's MD5 values as well as their retained SHA-256 manifests.
- Live TESS-EB catalogue: "This TESS ID is not in the catalog."
- The project's existing frozen catalogue index, with the core snapshot from
  2026-10-16 and older versioned support tables, finds no target association in
  NASA confirmed planets, TOIs, CTOIs, DIAmante, Prša 2022 EBs, Yu 2019,
  Olmschenk 2021, RAVEN, or the supplied manuscript.
- SIMBAD's proper-motion-star entry has one reference: the Gaia Catalogue of
  Nearby Stars (2021A&A...649A...6G). Finite arXiv metadata searches for the TIC,
  Gaia ID, 2MASS name and UCAC4 alias return zero results. This does not search
  every paper's full text, supplements, unpublished reports or all aliases.
- WDS and El-Badry et al. 2021 binary-catalogue queries within 30″ return no rows.

This is **no prior report found in the checked scope**, not a claim of novelty.
The source is already catalogued as a star. Public pipeline threshold-crossing
events and additional literature can exist without any TOI/CTOI entry.

## Crowding and plausible alternative sources

The 120″ Gaia cone contains **1,031 sources**. The TIC contamination ratio is
**2.229999** (catalogue estimated contaminating flux divided by target flux).
That implies approximately 31% target flux under the TIC estimate; it is not a
measurement of any particular QLP/SPOC aperture. Nearby sources are already
catalogued objects and must not be assigned to TIC 214311908.

| Neighbor | Separation | Relevant public information |
|---|---:|---|
| Gaia 5925373272827260544 / TIC 1515709327 | 3.15″ | G 17.48, T 16.85; unresolved within a TESS pixel |
| Gaia 5925373272827258752 / TIC 1515709326 | 4.27″ | G 18.60, T 17.82 |
| Gaia 5925373272827253632 / TIC 1515709322 | 11.29″ | G 15.74; incomplete astrometry, substantial excess-noise significance |
| Gaia 5925373272850782208 / TIC 214311910 | 11.67″ | G 14.65, T 13.37; near-comoving neighbor; RUWE 13.72 |
| Gaia 5925373268528684288 / TIC 214311906 | 17.78″ | Bright known Mira-like variable ASASSN-V J172354.64−511925.3; T 11.26; VSX period 294 d, Gaia LPV period about 326 d |
| Gaia 5925373272827278208 / TIC 1515709337 | 26.63″ | Known Gaia eclipsing binary, period 0.3160829620 ± 0.0000010405 d |
| Gaia 5925373204107739648 | 62.12″ | Known Gaia eclipsing binary, period 0.2609258493 ± 0.0000048504 d |

Gaia G flux ratios are not TESS flux ratios. The very red Mira is approximately
5.0 times the target's flux using the TIC T magnitudes, before aperture coupling;
its strongly varying brightness can influence sector normalization and detrending.
The nearby Gaia EBs have primary/secondary eclipse depths about 0.43/0.42 mag
and 0.23/0.18 mag respectively, before dilution. Their periods and harmonics
provide concrete contamination tests if a similar TESS period appears.

The 11.67″ neighbor has parallax 16.4126 ± 0.3929 mas and proper motion
(−112.7960 ± 0.4412, +87.2053 ± 0.3229) mas yr⁻¹, suggestive of a physical wide
pair with the target. Its high RUWE and 20% multiple-peak image fraction make
its astrometry uncertain. Association is not proved here. A projected separation
would be roughly 750 AU if both objects are at the target's distance. This would
not identify which star produces any transit-like signal. ExoFOP includes
**duplicate TIC 1515709367 for this same Gaia source**; do not count it as a
second separate contaminant. ExoFOP leaves the Mira's Gaia DR3 field blank but
retains the correct identical Gaia DR2 ID, so a DR3-only join would miss it.

The neighboring EB model reference epochs and phases are preserved in
[`known_neighbor_eb_ephemerides.json`](known_neighbor_eb_ephemerides.json).
Gaia's model times are barycentric JD **in TCB**, offset by 2455197.5 days;
exact comparison with TESS requires conversion to TDB and propagation of
period/phase uncertainty. An apparent phase mismatch years later is not by
itself a disproof of a blend.

## Reproduction and limitations

Run from the repository root:

```bash
python scripts/tic214311908_catalogs.py
```

The retrieval script uses only the Python standard library, verifies cached
responses, and retains failure bodies. Successful public snapshots include full
URLs or ADQL queries, UTC retrieval times, sizes and SHA-256 checksums under
`data/tic214311908/catalogs/`. Derived neighbor distances, full source IDs and
G-band flux ratios are in
[`gaia_neighbors_derived.csv`](gaia_neighbors_derived.csv); the machine-readable
audit is [`catalog_check.json`](catalog_check.json). The previous frozen-index
comparison and its exact input metadata are separately retained as
[`frozen_catalog_index.json`](frozen_catalog_index.json).

Pixel localization, recurrence, transit shape and false-positive diagnostics
must determine whether there is a credible signal and which source hosts it.
Catalogue nonmatches alone supply none of those tests.
