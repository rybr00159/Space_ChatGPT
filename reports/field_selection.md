# Independent cool-dwarf field selection

Selection completed at 2026-10-13T16:39:36.887888+00:00.

The field is centered at RA 100°, Dec -60° with radius 3°. TIC cuts: Tmag 9–11.5, Teff 3200–4800 K, radius 0.2–0.8 solar radii, finite contamination ≤ 0.05. Exclude all matches in current saved confirmed-planet and TOI catalogs, including dispositions other than planet candidate. Sort remaining stars by magnitude, contamination, then TIC ID. Check up to 20 stars, accept the first 6 with at least 4 actual SPOC `-s_lc.fits` observation products. No flux or BLS result influenced selection.

| TIC | Tmag | Teff (K) | Radius (solar) | TIC contamination | Sector count |
| --- | --- | --- | --- | --- | --- |
| 348839788 | 9.300 | 3797 | 0.556 | 0.0192 | 44 |
| 150359500 | 9.748 | 3286 | 0.365 | 0.0077 | 42 |
| 278727563 | 10.396 | 4736 | 0.663 | 0.0054 | 33 |
| 150360085 | 10.674 | 4160 | 0.658 | 0.0059 | 33 |
| 348839628 | 10.927 | 3626 | 0.482 | 0.0390 | 33 |
| 260476305 | 11.091 | 4245 | 0.709 | 0.0150 | 21 |

Selected 6 targets after checking 7 of 12 eligible catalog targets. Product URIs, sectors, exact request, catalog hashes, exclusion reasons, and observation-query statuses are saved in `data/field_sample.json`. Complete raw JSON responses remain in `/workspace/tess-research-data/field_queries`. FITS downloads are a separate pipeline step.

This is a comparison pilot; these targets are not proven to be unstudied. Before any novelty claim, check CTOIs, TCEs, aliases, and literature. Incomplete catalog coverage remains a limitation. TIC temperatures/radii/contamination are catalog estimates and require verification for survivors.

Reproduce with `python src/select_field_sample.py`. Existing exact-query cache entries avoid repeated archive traffic. Apply the same catalog snapshots or record updated hashes. A larger/different field must be logged as a strategy change.
