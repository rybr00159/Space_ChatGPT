# TIC 269994642: physical transit/eclipse shape assessment

This is a conditional shape analysis, **not planet validation or a mass measurement**. Finite exposure, grazing geometry, stellar assumptions and dilution materially affect the inferred companion radius.

Current TIC R=0.835257±0.0872364 Rsun, M=0.631±0.0830436 Msun, Teff=4035 K. Catalog contamination ratio 0.003582 is only a prior estimate; unresolved sources are not excluded.

1838 unbinned local samples from sectors [1, 92] enter the diagonal likelihood. Exact retained FITS TIMEDEL exposure values are integrated with 31 BATMAN subsamples; 61 subsamples provide a numerical convergence check. Baselines were fitted only outside protected events. Approximate error floors and baseline variance are included, but residual correlations and baseline covariance between samples are not fully modeled. Local normalization can absorb slow orbital light variations; their absence here cannot exclude a binary.

Circular, no-dilution, negligible-companion-mass best fit: P=1.746372124 d, epoch BTJD 1326.350048, k=0.2124, b=0.9168, a/R=6.106, T14=1.759 h, T23=0.000 h; grazing=True. Conditional radius=0.177 Rsun (1.73 Rjup). Data χ²=1931.2/1829 nominal dof, prior penalty=0.05.

The sampled radius profile at Δ penalized χ²≤3.84 admits k=0.175–0.800, or 1.43–6.50 Rjup for the adopted radius. Upper grid boundary reached: `True`. This upper limit is artificial when the boundary is reached. All other physical/timing parameters are reoptimized at each k. The plotted contours include broad approximate limb-darkening and log(a/R) priors. They are not calibrated confidence intervals or Bayesian posterior probabilities; conditional Gaussian covariance is unreliable along grazing degeneracies. Stellar radius uncertainty has not been marginalized into the profile radius range. Exact numerical contours are in `fit_results.json`.

| Scenario | k | b | Conditional Rcomp/Rjup | Data chi² | Prior chi² |
|---|---:|---:|---:|---:|---:|
| free_density | 0.2302 | 0.9464 | 1.87 | 1931.2 | 0.00 |
| dilution20percent | 0.2492 | 0.9315 | 2.03 | 1931.2 | 0.02 |
| dilution50percent | 0.3138 | 0.9210 | 2.55 | 1931.3 | 1.01 |
| massratio0.3 | 0.1938 | 0.8808 | 1.58 | 1931.4 | 0.63 |
| eccentric0.3_periastron | 0.8000 | 1.5946 | 6.50 | 1931.3 | 3.18 |
| eccentric0.5_periastron | 0.8000 | 1.5954 | 6.50 | 1931.4 | 7.55 |
| catalog_dilution | 0.2130 | 0.9170 | 1.73 | 1931.2 | 0.05 |

An explicit **equal-light binary** has physical orbit twice the photometric period (3.49274423 d) and two identical grazing eclipses. Individual radii are the combined-light inferred radius divided by √2, with two equal nominal stellar masses in Kepler's law. It fits b=1.704, a/R=15.87, T14=1.772 h, data χ²=1931.3, prior χ²=0.91. Equal primary/secondary eclipses at alternating photometric epochs evade both odd/even differences and a secondary test at photometric half-phase. Their combined photocenter can be astrometrically balanced. This alternative is conditional on equal Teff/brightness and illustrative mass/radius assumptions; it does not require detectable planetary-period secondary light. Compare its absolute data residuals and priors before any classification. No population odds or Bayes factors are inferred.

Allowing unequal surface brightness at fixed equal radii gives secondary/primary TESS brightness ratio J=0.8582±0.0230 (conditional diagonal covariance), data χ²=1898.2. The data improvement over equal brightness is Δχ²=33.13. This explicitly models alternating primary/secondary depths; it is not a calibrated population probability. The individual radius combined-light/√2 approximation becomes approximate for unequal brightness, and neither luminosity nor depth proves a stellar mass.

A separate fixed-physical-shape amplitude+constant fit to 24 individual events checks whether diagonal sample errors overstate this alternation. Error inflation includes per-event residual scatter, conservative even/odd event scatter, sector-to-sector scatter, and 3,000 event-cluster bootstrap draws within observed sectors. The conservative sector comparison is 3.33 conditional standard errors from zero; this is a **targeted diagnostic, not search-wide significance**. Sector and event measurements are in `parity_robustness.json` and `event_amplitude_measurements.csv`. With only two sectors, the TIC 269994642 bootstrap cannot estimate unobserved long-term systematics. No independent probabilities are obtained by multiplying tests that use the same photons. The alternating depths favor a binary interpretation while activity and baseline/shape mismatch remain limitations.

Dilution is a fraction of total aperture light; 0.2 and 0.5 are illustrative unresolved-source scenarios, not measured contaminations. A faint stellar companion can reproduce a grazing primary while its optical secondary remains weak. The total-system mass affects Kepler a/R; the mass-ratio 0.3 scenario tests that density effect without imposing a stellar mass-radius relation. Eccentric scenarios change transit velocity with omega=90 degrees and fit the primary alone; they are not orbit/secondary solutions or evidence of eccentricity. Both primary luminosity and unknown extra light make a dark-transit radius conditional. Radii of planets, brown dwarfs and faint stars overlap; mass requires other information.

The physical fits leave unresolved transit-like leads for further falsification. Follow source-specific odd/even and all-phase secondary tests, pixel localization, astrometric multiplicity, catalogs/literature and RV/high-resolution imaging before classification. **No new exoplanet is accepted on this shape analysis.**

Numerical integration check: max|31 vs 61 subsample model|=4.35e-05; Δ data χ²=0.2608. Exact source/input hashes and model arrays are preserved in this directory. Reproduce: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python src/fit_broad_leads.py --tic 269994642`.

The two observed eras initially imply 1,419 photometric cycles. Explicit 1,418/1,419/1,420-cycle fits on the same native photons give data χ² 1418: 2245.9, 1419: 1931.2, 1420: 2192.9. Adjacent aliases are therefore disfavored by the selected two-era shape fit. The formal period covariance conditions on this cycle mapping and stellar/shape assumptions; no general timing-variation model is fitted.

The doubled-period binary implementation also passes an independent native-time full-orbit test, adding separately calculated primary and secondary occultations: max flux discrepancy 6.7e-13. That verification checks the physical period/phase mapping; it does not validate a binary classification.
