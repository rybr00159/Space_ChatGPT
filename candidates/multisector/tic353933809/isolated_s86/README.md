# TIC 353933809: isolated Sector 86 event

**The selected QLP dip does not survive this phase as a robust, source-localized transit.** Native aperture shape, apparent depth, and timing depend strongly on background treatment and local baseline. A weaker real event remains possible; these checks do not establish its absence or an exact instrumental cause.

This is an isolated-event test around BTJD **3637.21848**, without assuming the rejected 18.9784-day recurrence period. Raw native pixels share photons with QLP, so these are processing and source diagnostics, not independent detection significance.

The native 9×9 TESSCut product contains 11,207 cadences, with 10,601 passing quality and finite-flux/error cuts. Full raw FITS, the downloaded ZIP, original FFI identifiers, timing, WCS, FITS-integrity results, byte hashes and science-array hashes are preserved. The exact URLs and raw paths are in [native_manifest.json](native_manifest.json) and [retrieval_manifest.json](retrieval_manifest.json). Gaia J2016 proper motions were propagated to the event date. The primary and catalogued companion are separated by only about 0.61 native pixels; their flux is strongly blended.

The first extraction preserved here used unsubtracted native pixels. Its shape fits hit depth/time boundaries, while a huge difference signal was an almost uniform additive sheet of roughly 480 electrons/s per pixel. Positive Gaussian widths hit bounds and source amplitude was negligible relative to the sheet. A high summed-pixel score in that extraction therefore supplies no localized transit evidence. These failed fits are retained in the JSON and diagnostic arrays directly in this directory.

Two subsequent sky treatments used the same explicitly defined peripheral pixels: radius above 3.1 pixels from the primary and out-of-event reference flux below its 45th percentile. One subtracts their median; the other fits a spatial plane. Both retain original raw bytes. Their diagonal error propagation and empirical residual floors still omit full sky covariance.

| Background | Half baseline window | Primary aperture shape | Selected-template difference S/N |
|---|---:|---|---:|
| Peripheral median | 0.8 d | Brightening: −336 ± 88 ppm | 4.19 |
| Peripheral median | 1.2 d | Dip: 347 ± 97 ppm, midpoint 3637.45185, duration 0.312 ± 0.066 d | 3.97 |
| Peripheral plane | 0.8 d | Dip: 613 ± 1422 ppm; ingress parameter near boundary | 5.04 |
| Peripheral plane | 1.2 d | Dip: 557 ± 480 ppm, midpoint 3637.38637, duration 0.795 ± 0.143 d | 5.27 |

These are selected, approximate exposure-integrated trapezoids with quadratic baselines. Both positive and negative depth signs were allowed. Other apertures often hit the duration or timing bounds, and their signs can reverse between baseline widths. Conditional 100-draw residual-block summaries are retained but do not account for sky-policy uncertainty, long instrumental correlations or search trials.

For each sky treatment, both difference-map windows use the template selected from that treatment's 1.2-day-window primary-aperture fit. They are not separately optimized maps. All 216 native cadences within 0.25 days of the initial epoch have quality zero; clean flags do not exclude background or pointing errors.

The remaining difference maps contain positive and negative lobes around the bright stellar image. A single empirical Gaussian fits poorly: diagonal χ² approximately 316–596 for 31 degrees of freedom. The summed-aperture score is below the predeclared threshold of 7 for attempting a localization bootstrap. **No Gaussian centroid, host assignment, or positional confidence is claimed.**

A separate diagnostic fits the difference image as a scaled out-of-event reference image, its x/y gradients, and a spatial background plane. The gradients allow pointing-like redistribution without requiring stellar flux loss. After this decomposition the nominal fractional flux decrements are:

| Background | Half window | Reference-image decrement |
|---|---:|---:|
| Median | 0.8 d | 468 ± 374 ppm |
| Median | 1.2 d | 573 ± 478 ppm |
| Plane | 0.8 d | 824 ± 694 ppm |
| Plane | 1.2 d | 895 ± 725 ppm |

The quoted errors include a residual-scaled diagonal covariance. The decomposition itself remains poor (χ² about 281–501 for 33 degrees of freedom), so it does not identify a complete instrumental model. It demonstrates that a strong, well-measured stellar flux decrement is not robust to plausible image redistribution. A blended genuine event can also project onto gradients; this is not a formal exclusion test.

The published Gaia AstroSpectroSB1 orbit (1124.8 days) and the nearby Gaia companion are prior context. The nominal photocentric inclination inferred from the published Thiele–Innes coefficients is about 29.4°, with no covariance propagated here. That outer orbit does not exclude an inner planet, a different eclipsing source, or other astrophysical variability, and this isolated dip has not been assigned to it.

Inspect [median-sky shapes and maps](median/isolated_pixel_validation.json), [plane-sky shapes and maps](rim_plane/isolated_pixel_validation.json), and the corresponding `difference_redistribution_test.json` files. The plots retain failed and alternative fits as well as fitted values.

To reproduce from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/pooled_north_isolated_pixels.py --background median
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/pooled_north_isolated_pixels.py --background rim_plane
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/pooled_north_isolated_pixels.py --background median --decompose-only
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/pooled_north_isolated_pixels.py --background rim_plane --decompose-only
```

The script verifies retained raw hashes and preserves the native archive. Source snapshots associated with each result remain under `method_snapshots/`. No physical planet fit, companion radius, false-positive probability, mass or novelty claim is justified.

Regenerated TESSCut ZIP/FITS headers can change. If a missing archive is retrieved and differs from its saved expected hash, the script preserves the unexpected response separately and stops without overwriting the original manifest. Compare the recorded science-array hash explicitly before scientific reuse; do not substitute a new expected byte hash simply to make verification pass.
