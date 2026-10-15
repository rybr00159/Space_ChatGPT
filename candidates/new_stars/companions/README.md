# Candidate associations among catalogued Gaia stars

**Pairs 1 and 3 are plausible wide stellar associations with no recognition in
the finite checked scope. Neither gravitational binding nor first publication
is established. Both component sky sources were already catalogued.**

| Pair | Gaia DR3 primary / neighbor | RA / Dec of primary, degrees | Separation | Nominal projected separation |
|---|---|---|---|---|
|1|4977459112829640832 /4977459112829654400|4.52642935 /−47.54760148|4.03386 arcsec|7134 AU |
|3|4923799509418222080 /4923799509418221952|2.46128905 /−54.22809809|4.76273 arcsec|4889 AU |

The first pair has G 14.76/16.21 and BP−RP 0.728/0.820 at a nominal 1.77 kpc.
The second has G 16.33/17.83 and colors 1.262/1.845 at 1.03 kpc. Both have
reasonable stellar color/luminosity combinations and clean Gaia diagnostics.
They have separate images in 1993/1996 photographic plates. Plate fluxes are
nonlinear and the images do not independently prove matching distances/motions.

Nominal relative transverse velocities 0.035/0.242 km/s are smaller than their
uncertainty scale. Propagated minimum-binding-mass 16/50/84 percentiles are
0.28/1.14/3.10 and 0.24/0.97/2.65 solar masses. These are measurement-error
transforms, not measured masses or binding probabilities. Radial velocities
are unavailable. An unbound comoving association remains possible.

[Historical images](historical_images.png) show both components;
[projected dynamics](projected_dynamics.png) is a conditional geometry display.
The [reviewed pairs](reviewed_pairs.csv), [full measurements](leads_full.json),
[quality/scope audit](quality_and_catalog_scope.json), [summary](summary.json)
and [raw-response provenance](provenance.json) preserve all eleven cases.

Pairs 2/5 have compromised photometry or astrometry;9/11 have weak dynamics.
Five other pairs are already known. Exact DR2 mapping, catalogue positive
controls and three shifted fields were checked. Null Tian/El-Badry/SUPERWIDE
results supply no novelty evidence because these candidates fail those
catalogues' selections. WDS/SIMBAD absence is finite, not exhaustive literature
clearance. See [the full report](../../../NEW_STAR_SEARCH_REPORT.md).

Reproduce from retained data:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python src/new_companion_search.py
# --download restores missing public responses; failed query variants remain.
```

The program verifies raw hashes, known-pair recovery, alias completeness and
catalogue-selection failures. Next useful evidence is differential spectroscopy
and broader prior-source review. No submissions or observing requests were sent.
