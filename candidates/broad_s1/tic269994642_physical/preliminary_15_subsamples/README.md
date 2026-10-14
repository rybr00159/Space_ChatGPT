# TIC 269994642: physical transit/eclipse shape assessment

This is a conditional shape analysis, **not planet validation or a mass measurement**. Finite exposure, grazing geometry, stellar assumptions and dilution materially affect the inferred companion radius.

CurrentTIC R=0.835257±0.0872364Rsun,M=0.631±0.0830436Msun,Teff4035K. Catalog contamination ratio0.003582 is only a prior estimate; unresolved sources are not excluded.

1838unbinned local samples fromsectors[1, 92] enter the diagonal likelihood. Exact retained FITS TIMEDEL exposure values are integrated with15BATMAN subsamples;31subsamples provide a numerical convergence check. Baselines were fitted only outside protected events. Approximate error floors and baseline variance are included, but residual correlations and baseline covariance between samples are not fully modeled. Local normalization can absorb slow orbital light variations; their absence here cannot exclude a binary.

Circular, no-dilution, negligible-companion-mass bestfit: P=1.746372123d, epochBTJD1326.350049, k=0.2128, b=0.9174, a/R=6.113, T14=1.756h, T23=0.000h; grazing=True. Conditional radius=0.178Rsun (1.73Rjup). Dataχ²=1930.7/1829nominaldof, priorpenalty=0.05.

Sampled radius profile: `{'minimum_penalized_chi2': 1930.7302623523658, 'grid_min': 0.1, 'grid_max': 0.8, 'delta_chi2_le_1': {'sampled_k_range': [0.18, 0.8], 'upper_bound_reached': True, 'conditional_radius_rjup_range': [1.463038836935473, 6.502394830824325]}, 'delta_chi2_le_3.84': {'sampled_k_range': [0.17541970440805751, 0.8], 'upper_bound_reached': True, 'conditional_radius_rjup_range': [1.425810223959605, 6.502394830824325]}, 'warning': 'Conditional sampled penalized likelihood profile; no calibrated coverage, covariance or mass classification. Fixed stellar-density/limb priors and diagonal errors; grazing upper bounds are artificial.'}`. All other physical/timing parameters are reoptimized at each k. The plotted contours include broad approximate limb-darkening and log(a/R)priors. They are not calibrated confidence intervals or Bayesian posterior probabilities; conditional Gaussian covariance is unreliable along grazing degeneracies. Stellar radius uncertainty has not been marginalized into the profile radius range.

| Scenario | k | b | Conditional Rcomp/Rjup | Data chi² | Prior chi² |
|---|---:|---:|---:|---:|---:|
| free_density | 0.2290 | 0.9445 | 1.86 | 1930.7 | 0.00 |
| dilution20percent | 0.2493 | 0.9318 | 2.03 | 1930.7 | 0.02 |
| dilution50percent | 0.3135 | 0.9205 | 2.55 | 1930.8 | 1.04 |
| massratio0.3 | 0.1941 | 0.8814 | 1.58 | 1930.8 | 0.62 |
| eccentric0.3_periastron | 0.8000 | 1.5947 | 6.50 | 1930.8 | 3.15 |
| eccentric0.5_periastron | 0.8000 | 1.5955 | 6.50 | 1930.9 | 7.51 |
| catalog_dilution | 0.2133 | 0.9176 | 1.73 | 1930.7 | 0.04 |

An explicit **equal-light binary** has physical orbit twice the photometric period (3.49274423d) and two identical grazing eclipses. Individual radii are the combined-light inferred radius divided by√2, with two equal nominal stellar masses in Kepler's law. It fits b=1.704,a/R=15.89,T14=1.769h, dataχ²=1930.8,priorχ²=0.88. Equal primary/secondary eclipses at alternating photometric epochs evade both odd/even differences and a secondary test at photometric half-phase. Their combined photocenter can be astrometrically balanced. This alternative is conditional on equal Teff/brightness and illustrative mass/radius assumptions; it does not require detectable planetary-period secondary light. Compare its absolute data residuals and priors before any classification. No population odds or Bayes factors are inferred.

Dilution is a fraction of total aperture light;0.2and0.5 are illustrative unresolved-source scenarios, not measured contaminations. A faint stellar companion can reproduce a grazing primary while its optical secondary remains weak. The total-system mass affects Kepler a/R; the mass-ratio0.3 scenario tests that density effect without imposing a stellar mass-radius relation. Eccentric scenarios change transit velocity with omega90degrees and fit the primary alone; they are not orbit/secondary solutions or evidence of eccentricity. Both primary luminosity and unknown extra light make a dark-transit radius conditional. Radii of planets, brown dwarfs and faint stars overlap; mass requires other information.

The physical fits leave unresolved transit-like leads for further falsification. Follow source-specific odd/even and all-phase secondary tests, pixel localization, astrometric multiplicity, catalogs/literature and RV/high-resolution imaging before classification. **No new exoplanet is accepted on this shape analysis.**

Numerical integration check: max|15vs31subsamplemodel|=9.88e-05; Δdataχ²=0.5574. Exact source/input hashes and model arrays are preserved in this directory. Reproduce: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python src/fit_broad_leads.py --tic 269994642`.
