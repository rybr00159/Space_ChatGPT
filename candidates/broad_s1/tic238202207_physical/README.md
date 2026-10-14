# TIC 238202207: physical transit/eclipse shape assessment

This is a conditional shape analysis, **not planet validation or a mass measurement**. Finite exposure, grazing geometry, stellar assumptions and dilution materially affect the inferred companion radius.

TIC/FLAME R≈1.63 Rsun versus Gaia GSP-Phot R=1.1929 Rsun. TIC radius/mass uncertainties are absent, so 10% radius/15% mass widths are explicit scenario assumptions. Gaia RUWE=4.128, IPDmulti_peak=22% undermine a clean single-star interpretation.

4490 unbinned local samples from sectors [1, 2, 5, 6, 9, 10, 12, 13, 27, 28, 29, 30, 32, 33, 35, 36, 39] enter the diagonal likelihood. Exact retained FITS TIMEDEL exposure values are integrated with 31 BATMAN subsamples; 61 subsamples provide a numerical convergence check. Baselines were fitted only outside protected events. Approximate error floors and baseline variance are included, but residual correlations and baseline covariance between samples are not fully modeled. Local normalization can absorb slow orbital light variations; their absence here cannot exclude a binary.

Circular, no-dilution, negligible-companion-mass best fit: P=2.321825581 d, epoch BTJD 1327.599821, k=0.3024, b=1.1219, a/R=5.701, T14=2.104 h, T23=0.000 h; grazing=True. Conditional radius=0.492 Rsun (4.79 Rjup). Data χ²=3857.5/4466 nominal dof, prior penalty=4.30.

The sampled radius profile at Δ penalized χ²≤3.84 admits k=0.150–0.800, or 2.38–12.67 Rjup for the adopted radius. Upper grid boundary reached: `True`. This upper limit is artificial when the boundary is reached. All other physical/timing parameters are reoptimized at each k. The plotted contours include broad approximate limb-darkening and log(a/R) priors. They are not calibrated confidence intervals or Bayesian posterior probabilities; conditional Gaussian covariance is unreliable along grazing degeneracies. Stellar radius uncertainty has not been marginalized into the profile radius range. Exact numerical contours are in `fit_results.json`.

| Scenario | k | b | Conditional Rcomp/Rjup | Data chi² | Prior chi² |
|---|---:|---:|---:|---:|---:|
| free_density | 0.1631 | 0.9216 | 2.58 | 3857.1 | 0.01 |
| dilution20percent | 0.3118 | 1.1051 | 4.94 | 3857.6 | 6.65 |
| dilution50percent | 0.3279 | 1.0494 | 5.19 | 3858.0 | 13.21 |
| massratio0.3 | 0.2863 | 1.1010 | 4.53 | 3857.5 | 1.86 |
| eccentric0.3_periastron | 0.1643 | 0.9245 | 2.60 | 3857.1 | 0.04 |
| eccentric0.5_periastron | 0.8000 | 1.6532 | 12.67 | 3857.8 | 0.15 |
| gaia_radius1.1929 | 0.1625 | 0.9204 | 1.89 | 3857.1 | 0.01 |
| gaia_radius_dilution50percent | 0.2460 | 0.9295 | 2.86 | 3857.5 | 1.39 |

An explicit **equal-light binary** has physical orbit twice the photometric period (4.64365114 d) and two identical grazing eclipses. Individual radii are the combined-light inferred radius divided by √2, with two equal nominal stellar masses in Kepler's law. It fits b=1.785, a/R=15.26, T14=2.112 h, data χ²=3857.7, prior χ²=2.54. Equal primary/secondary eclipses at alternating photometric epochs evade both odd/even differences and a secondary test at photometric half-phase. Their combined photocenter can be astrometrically balanced. This alternative is conditional on equal Teff/brightness and illustrative mass/radius assumptions; it does not require detectable planetary-period secondary light. Compare its absolute data residuals and priors before any classification. No population odds or Bayes factors are inferred.

Allowing unequal surface brightness at fixed equal radii gives secondary/primary TESS brightness ratio J=0.8930±0.0194 (conditional diagonal covariance), data χ²=3830.1. The data improvement over equal brightness is Δχ²=27.55. This explicitly models alternating primary/secondary depths; it is not a calibrated population probability. The individual radius combined-light/√2 approximation becomes approximate for unequal brightness, and neither luminosity nor depth proves a stellar mass.

A separate fixed-physical-shape amplitude+constant fit to 140 individual events checks whether diagonal sample errors overstate this alternation. Error inflation includes per-event residual scatter, conservative even/odd event scatter, sector-to-sector scatter, and 3,000 event-cluster bootstrap draws within observed sectors. The conservative sector comparison is 4.39 conditional standard errors from zero; this is a **targeted diagnostic, not search-wide significance**. Sector and event measurements are in `parity_robustness.json` and `event_amplitude_measurements.csv`. With only two sectors, the TIC 269994642 bootstrap cannot estimate unobserved long-term systematics. No independent probabilities are obtained by multiplying tests that use the same photons. The alternating depths favor a binary interpretation while activity and baseline/shape mismatch remain limitations.

Dilution is a fraction of total aperture light; 0.2 and 0.5 are illustrative unresolved-source scenarios, not measured contaminations. A faint stellar companion can reproduce a grazing primary while its optical secondary remains weak. The total-system mass affects Kepler a/R; the mass-ratio 0.3 scenario tests that density effect without imposing a stellar mass-radius relation. Eccentric scenarios change transit velocity with omega=90 degrees and fit the primary alone; they are not orbit/secondary solutions or evidence of eccentricity. Both primary luminosity and unknown extra light make a dark-transit radius conditional. Radii of planets, brown dwarfs and faint stars overlap; mass requires other information.

The physical fits leave unresolved transit-like leads for further falsification. Follow source-specific odd/even and all-phase secondary tests, pixel localization, astrometric multiplicity, catalogs/literature and RV/high-resolution imaging before classification. **No new exoplanet is accepted on this shape analysis.**

Numerical integration check: max|31 vs 61 subsample model|=2.29e-05; Δ data χ²=0.07404. Exact source/input hashes and model arrays are preserved in this directory. Reproduce: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python src/fit_broad_leads.py --tic 238202207`.

The doubled-period binary implementation also passes an independent native-time full-orbit test, adding separately calculated primary and secondary occultations: max flux discrepancy 1.1e-13. That verification checks the physical period/phase mapping; it does not validate a binary classification.
