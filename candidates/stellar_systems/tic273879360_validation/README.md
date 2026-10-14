# Independent physical validation of TIC 273879360

Checked 2026-10-14 UTC. **The eclipses support an eccentric luminous binary
near the catalogued background source. Prior public neural-network eclipse
recognition exists; this is independent characterization, not first discovery.**

The source is Gaia DR3 4692910966756908288 at RA 32.252462825 degrees,
Dec−71.323170344 degrees. Both eclipse families were independently localized
near it in [the earlier source dossier](../tic273879360/README.md). The bright
foreground star TIC 273879358 is 65.59 arcseconds away; its radial velocities
must not be assigned to this source.

This new extraction uses the four original native FITS cubes, QUALITY=0,
an independent source-core aperture and blank-sky estimate. Twenty-two
coverage-qualified event windows contain 6234 cadences. Per-event simultaneous
raw baselines and two BATMAN occultations represent two luminous spheres in
one eccentric orbit, with sector-specific dilution and native exposure
integration. No earlier flattened light curve is used in the physical fit.

Main conditional solution: P 7.586890837 days, primary BTJD 1327.900292,
secondary phase 0.691727, e 0.32145, inclination 88.48 degrees, R2/R1 ≈ 0.680,
a/R1 ≈ 20.63, TESS L2/L1 ≈ 0.145. These are dimensionless eclipse-model parameters;
catalog single-star radii/temperatures do not measure either component.
Dilution, grazing geometry, approximate limb darkening and source calibration
prevent absolute radii/masses or a secure mass ratio. The conjunction
approximation differs from exact minimum projected separation by ≈ 0.68 minutes.

| Artifact | Purpose |
|---|---|
| [physical_fit.json](physical_fit.json) | Eight initial geometries, best solution, dark/circular alternatives, exact raw manifests and exclusions |
| [physical_fit.png](physical_fit.png) | Raw local eclipse shapes, physical model and empirical-error residuals |
| [physical_controls.json](physical_controls.json) | Aperture, polynomial and limb tests; exposure convergence; failed epoch holdouts; 24 residual-block draws |
| [alias_and_period_tests.json](alias_and_period_tests.json) | Missing-eclipse windows for half-period hypothesis and broadened early-period likelihood profile |
| [alias_and_period_tests.png](alias_and_period_tests.png) | Early period profile and signed extra-eclipse template depths |
| raw_model_arrays.npz | Compact exact timestamps, raw flux/errors, model, baseline and residuals |

The bootstrap resamples 30-minute residual blocks within each existing event,
not whole sectors/events; it omits shared calibration and broad model
uncertainty. Conditional e 16th–84th percentile range 0.3197–0.3222 is not a
complete error bar. The radius-ratio interval 0.664–0.738 is wider than changes
in the six explicit extraction/limb choices. Dark and circular penalties
mostly express the same missing-secondary failure and must not be added as
independent evidence. No survey-wide false-positive probability is supplied.

The early-only prediction failed when an unconverged restricted period fit
was treated as precise. Broadened, converged profile fits quantify the weak
early timing constraint; any residual mismatch remains documented. Half-period
controls test covered phases where extra eclipses should appear, rather than
only comparing fits inside selected real events. The initial coarse-profile
outputs `alias_and_period_tests_v1.*`, a second plot with an incorrect period
offset axis (`alias_and_period_tests_v2.*`), and exploratory `epoch_early_fit.json` /
`epoch_later_fit.json` are retained as superseded diagnostics; use the current
profile and fit-status fields for interpretation. The version2 method is
retained in `methods/` with its original source identity. Its numerical
profile is unchanged in the current corrected plot.

Archived [ASAS-SN photometry](../deep_context/ground_photometry/README.md) provides
nine independent primary events, eight before TESS. Its V-band depth,
field-center HJD timing and unresolved blending are handled separately from
the TESS physical orbit. [Context and novelty audit](../../../reports/stellar_systems/deep_context.md)
explains the prior Kostov 2025 Table 2 eclipse flag and the limitations of
Gaia assumed-binary fits. [The consolidated report](../../../STELLAR_VALIDATION_REPORT.md)
compares both systems and follow-up needs.

Reproduce from the repository using the prepared astronomy environment:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/workspace/tess-research-data/matplotlib
python src/validate_background_binary.py
python src/background_binary_controls.py
python src/background_binary_alias_tests.py
python scripts/stellar_background_verify.py
```

Raw FITS remain outside Git. Their original URLs, archive identifiers and
hashes are included in `physical_fit.json` and the earlier sector manifests.
Restore missing cubes using the previous pixel configurations listed in
PROGRESS.md, preserving mismatches for investigation. Source RA/Dec, exposure
times, normalization, measured noise floors and event rejections are recorded.
The cloud overlay is not an independently verified durable volume. No paid
API calls or external submissions were made.
