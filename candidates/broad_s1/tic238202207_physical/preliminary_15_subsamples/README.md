# TIC 238202207: physical transit/eclipse shape assessment

This is a conditional shape analysis, **not planet validation or a mass measurement**. Finite exposure, grazing geometry, stellar assumptions and dilution materially affect the inferred companion radius.

TIC/FLAME R≈1.63Rsun versus GaiaGSP-PhotR=1.1929Rsun. TIC radius/mass uncertainties are absent, so 10%radius/15%mass widths are explicit scenario assumptions. GaiaRUWE4.128,IPDmulti_peak22% undermine a clean single-star interpretation.

4490unbinned local samples fromsectors[1, 2, 5, 6, 9, 10, 12, 13, 27, 28, 29, 30, 32, 33, 35, 36, 39] enter the diagonal likelihood. Exact retained FITS TIMEDEL exposure values are integrated with15BATMAN subsamples;31subsamples provide a numerical convergence check. Baselines were fitted only outside protected events. Approximate error floors and baseline variance are included, but residual correlations and baseline covariance between samples are not fully modeled. Local normalization can absorb slow orbital light variations; their absence here cannot exclude a binary.

Circular, no-dilution, negligible-companion-mass bestfit: P=2.321825569d, epochBTJD1327.599826, k=0.3021, b=1.1215, a/R=5.708, T14=2.101h, T23=0.000h; grazing=True. Conditional radius=0.492Rsun (4.78Rjup). Dataχ²=3857.4/4466nominaldof, priorpenalty=4.36.

Sampled radius profile: `{'minimum_penalized_chi2': 3861.712898127833, 'grid_min': 0.1, 'grid_max': 0.8, 'delta_chi2_le_1': {'sampled_k_range': [0.18556078980696014, 0.5398027896111792], 'upper_bound_reached': False, 'conditional_radius_rjup_range': [2.9387485630736103, 8.548921752075222]}, 'delta_chi2_le_3.84': {'sampled_k_range': [0.14820227227359112, 0.8], 'upper_bound_reached': True, 'conditional_radius_rjup_range': [2.347097224264583, 12.669696291466774]}, 'warning': 'Conditional sampled penalized likelihood profile; no calibrated coverage, covariance or mass classification. Fixed stellar-density/limb priors and diagonal errors; grazing upper bounds are artificial.'}`. All other physical/timing parameters are reoptimized at each k. The plotted contours include broad approximate limb-darkening and log(a/R)priors. They are not calibrated confidence intervals or Bayesian posterior probabilities; conditional Gaussian covariance is unreliable along grazing degeneracies. Stellar radius uncertainty has not been marginalized into the profile radius range.

| Scenario | k | b | Conditional Rcomp/Rjup | Data chi² | Prior chi² |
|---|---:|---:|---:|---:|---:|
| free_density | 0.1620 | 0.9194 | 2.57 | 3856.9 | 0.01 |
| dilution20percent | 0.3113 | 1.1046 | 4.93 | 3857.5 | 6.72 |
| dilution50percent | 0.3278 | 1.0494 | 5.19 | 3857.8 | 13.31 |
| massratio0.3 | 0.2851 | 1.0996 | 4.51 | 3857.3 | 1.90 |
| eccentric0.3_periastron | 0.1635 | 0.9230 | 2.59 | 3856.9 | 0.05 |
| eccentric0.5_periastron | 0.8000 | 1.6532 | 12.67 | 3857.7 | 0.15 |
| gaia_radius1.1929 | 0.1618 | 0.9191 | 1.88 | 3856.9 | 0.01 |
| gaia_radius_dilution50percent | 0.2444 | 0.9267 | 2.84 | 3857.3 | 1.44 |

An explicit **equal-light binary** has physical orbit twice the photometric period (4.64365113d) and two identical grazing eclipses. Individual radii are the combined-light inferred radius divided by√2, with two equal nominal stellar masses in Kepler's law. It fits b=1.785,a/R=15.28,T14=2.109h, dataχ²=3857.5,priorχ²=2.59. Equal primary/secondary eclipses at alternating photometric epochs evade both odd/even differences and a secondary test at photometric half-phase. Their combined photocenter can be astrometrically balanced. This alternative is conditional on equal Teff/brightness and illustrative mass/radius assumptions; it does not require detectable planetary-period secondary light. Compare its absolute data residuals and priors before any classification. No population odds or Bayes factors are inferred.

Dilution is a fraction of total aperture light;0.2and0.5 are illustrative unresolved-source scenarios, not measured contaminations. A faint stellar companion can reproduce a grazing primary while its optical secondary remains weak. The total-system mass affects Kepler a/R; the mass-ratio0.3 scenario tests that density effect without imposing a stellar mass-radius relation. Eccentric scenarios change transit velocity with omega90degrees and fit the primary alone; they are not orbit/secondary solutions or evidence of eccentricity. Both primary luminosity and unknown extra light make a dark-transit radius conditional. Radii of planets, brown dwarfs and faint stars overlap; mass requires other information.

The physical fits leave unresolved transit-like leads for further falsification. Follow source-specific odd/even and all-phase secondary tests, pixel localization, astrometric multiplicity, catalogs/literature and RV/high-resolution imaging before classification. **No new exoplanet is accepted on this shape analysis.**

Numerical integration check: max|15vs31subsamplemodel|=5.2e-05; Δdataχ²=0.1426. Exact source/input hashes and model arrays are preserved in this directory. Reproduce: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python src/fit_broad_leads.py --tic 238202207`.
