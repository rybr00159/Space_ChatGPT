# Independent prior-report context for pooled-sector leaders

Retrieved 16 October 2026 UTC. These are finite recognition checks, not discoveries
or independent confirmations of the photometric signals. The ten initial targets
were selected from the already-completed pooled search; TIC 238182368 was added
after independent southern-sector vetting. This target-based follow-up has no
unbiased occurrence-rate or all-sky completeness interpretation.

## Findings that change interpretation

* **TIC 375090561:** already reported by Yu et al. (2019, AJ 158, 25) at
  5.423940 days. Actual public combined SPOC DV XML fits 5.42428305 days,
  approximately 4,661 ppm, and also contains a 2.7123614-day TCE. The new pooled
  8.136209-day peak is almost exactly **3/2 of the published period**. The frozen
  integer-only harmonic matcher did not label this rational alias as a harmonic;
  it did retain the previously reported host. This lead must not be interpreted
  as a new planetary signal without separating the already-reported event train.
  The target also has a Gaia `Acceleration9` solution and TIC contamination ratio
  1.308345. These facts do not establish which source produces the TESS events.
* **TIC 382190230:** explicitly classified **“RUWE Binary”** by Avallone et al.
  (2022, ApJ 930, 7), table 1, with published v sin i 75.60 km/s and radius
  1.35 solar radii. This is prior multiplicity recognition, not proof of an orbit.
  Its Gaia RUWE is 3.04 and radial-velocity variability is significant. Public
  SPOC DV has three TCEs near 0.84260–0.84279 days. Their existence is prior
  detection context; it does not identify the pooled 1.28239-day peak as the same
  signal. No rotation-period row for this target occurs in Avallone table 2.
* **TIC 353933809:** has a published Gaia **AstroSpectroSB1** orbit at
  1124.79934 ± 38.11657 days and RUWE 10.75. A nearby faint source has similar
  parallax and proper motion. The 18.9784-day pooled peak is not identified with
  that long-period orbit by these catalog checks. Phase extrapolation would need
  the orbital covariance and correct Gaia time reference.
* **TIC 349972227:** already a same-position VSX **RS** variable (0.38 arcsec
  separation), independently classified RS in Gaia DR3 with score 0.6723. Its
  Gaia robust RV range is 150.29 km/s with published chi-square p-value zero to
  catalog precision. This is strong stellar-variability context, not a new star
  discovery or a companion mass measurement.
* **TIC 359679661, 366490373 and 293268396:** have significant Gaia RV
  variability and/or poor astrometric fits. Particularly, RUWE is 17.33 for
  366490373, and robust RV range is 81.13 km/s for 293268396. Such diagnostics
  are warnings, not automatic eclipse-source assignments or planet exclusions.

The other targets' catalog nulls do not strengthen a novelty claim. All ten
initial public ExoFOP JSON pages had empty TOI, CTOI and file lists, despite the
real prior TESS report for TIC 375090561. The pages therefore cannot replace
prior literature or SPOC-product checks. The JSON schema does not constitute an
inventory of every comment, proprietary note or pipeline TCE.

## Evidence and scope

* [target_context.csv](target_context.csv): searchable ten-target identities,
  stellar properties, warning flags, and individual assessments.
* [review.json](review.json): actual target JSON, SIMBAD aliases and source-tagged
  bibliography, Gaia 60-arcsec neighbors and target NSS/variability checks,
  VSX 30-arcsec cones, ASAS-SN 2021 60-arcsec cones, and exact-TIC arXiv checks.
* [alias_resolution.json](alias_resolution.json): official Gaia DR3-to-DR2
  neighborhood rows agree with all ten ExoFOP identifier pairs; 2MASS aliases
  retained. Angular-distance values keep their source catalog definition and are
  not treated as an eclipse localization.
* [rotation_literature.json](rotation_literature.json): actual published
  Avallone classification row, from complete 4,216-row table 1 and 169-row table 2.
* [spoc_dv_inventory.json](spoc_dv_inventory.json): 136 exact-target SPOC CAOM
  observations including positive control TIC 55525572; 344 product rows,
  including 90 DV product rows. The initial targets with returned DV products
  were 375090561 and 382190230. This finite product query is not a complete
  inventory of every team's TCEs, every processing version or every publication.
* [spoc_dv_xml.json](spoc_dv_xml.json): numerical fits and diagnostics from
  **11 real public XML reports**, including a positive TOI-813 b recovery near
  83.8963 days. TCEs are unvalidated pipeline detections. The original XML bodies
  and all HTTP manifests remain in the referenced timestamped data directory.
* [extra_tic238182368/](extra_tic238182368/): additional mixed photometric lead.
  Gaia RUWE 0.961, no target NSS solution or strong RV variability, and valid
  empty finite VSX/ASAS-SN cones. TIC radius 1.1067 solar radii, Teff 5346 K and
  contamination ratio 0.203952. A pair of G≈14 stars lies about 39 arcsec away.
  The lack of recognition does not validate its mixed-sector photometry.

Nearby SIMBAD bibliography is explicitly tagged by source and separation.
For example, red-giant publications near TIC 287086419 and a Gaia paper near
TIC 356821129 must not be assigned to the selected target. No searched target's
RS/EB/variable classification alone proves the origin of a TESS dip.

The first VSX request for TIC 287086419 returned HTTP 503; the failure body is
retained. A new valid response returned no rows. ASAS-SN cones all returned valid
schema with no rows. The remote positive-control check recovered the actual
published source J030022.84+002836.8. An initial attempted ASAS-SN control used a
VSX-named ASAS-SN variable absent from this particular published 2021 table;
that unsuccessful control and its response remain preserved. It was replaced
with a verified member of the table, not reinterpreted as an empty service.
The arXiv positive request recovered RAVEN arXiv:2603.22597v2. Exact-TIC searches
returned no entries; this narrow text scope plainly misses real alias-based
literature and is not an exhaustive publication search.

## Reproduction and preservation

All HTTP bodies, including failures, have URL/request, status and SHA-256
manifests. Existing source and older results were not modified. The reusable
driver is [scripts/pooled_catalog_review.py](../../../scripts/pooled_catalog_review.py).
Its calls use the frozen `src/multisector_recognition.py`. No paid requests or
external submissions were made.

Re-derive the current retained ten-target table and compact DV parameters:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/pooled_catalog_review.py summarize
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/pooled_catalog_review.py report
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/pooled_catalog_review.py dvsummary
```

Commands `fetch`, `supplement`, `dvxml`, `rotation`, `aliases`, `remote-controls`,
`extra238`, and `tce238` retrieve new dated evidence rather than changing old HTTP
files. The `controls` command checks retained recognition and parser controls.
The redundant full XML parse tree was preserved outside ordinary Git; its
location and checksum are in `full_xml_parse_preservation.json`. Compact derived
fits and exact original XML are sufficient to reproduce the scientific summary.
