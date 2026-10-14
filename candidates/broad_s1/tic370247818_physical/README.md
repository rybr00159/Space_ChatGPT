# TIC 370247818: physical transit/eclipse assessment

This is a conditional light-curve fit, **not a mass measurement or planet
validation**. A grazing stellar companion remains an allowed explanation.
No physical radius is inferred by taking the square root of observed depth.

Eight independent TESS sectors contribute 2099 locally normalized
samples around the recurring 3.522-day signal. Exact raw FITS TIMEDEL values
are used: 30 minutes for sectors 1/5/12/13 and 10 minutes for 27/28/32/39.
BATMAN quadratic limb darkening is numerically integrated across each exposure
with 25 sub-samples. Current TIC R=1.68401±0.08391 solar radii and
M=1.28±0.20331 solar masses anchor a conservative Gaussian log(a/R) prior;
the reported TIC density uncertainty is used if larger than independent M/R
propagation. Free, physically permitted limb-darkening coefficients have broad
approximate priors u1=0.30±0.15 and u2=0.25±0.20. These are approximate TESS
priors, not a matched atmosphere-grid calculation.

Fitted parameters include k=Rcomp/Rstar over 0.03--0.5, b/(1+k) over
0--0.9999, a/R, both limb-darkening coefficients, period, epoch, and profiled
per-sector constant baselines. The prior event-timing ephemeris initializes
the fit; it is **not reused as an independent prior** because it uses the same
photons. No transit phase binning is used in the likelihood. Errors inherited
from local fits include observational scatter and approximate baseline variance;
their remaining correlations are not represented by the diagonal likelihood.

Circular, undiluted, negligible-companion-mass best fit:

| Parameter | Value |
|---|---:|
| Period d | 3.522195916 |
| Epoch BTJD | 1326.916389 |
| Radius ratio k | 0.27782 |
| Impact parameter b | 1.15807 |
| a/Rstar | 6.4899 |
| Inclination deg | 79.7210 |
| T14 h | 2.2785 |
| T23 h | 0.0000 |
| u1, u2 | 0.2698, 0.2007 |
| Companion radius, conditional Rsun | 0.4679 |
| Companion radius, conditional Rjup | 4.553 |
| Data chi² / nominal dof | 1522.6 / 2084 |
| Prior penalty chi² | 0.263 |

Grazing status is `True`: b>1-k removes the fully interior phase
and makes depth a degenerate function of radius, limb darkening, and geometry.
The best-fit radius above assumes the target star is the transited source,
the stated stellar radius, no dilution, a dark companion, and circular orbit.
It is not evidence of planetary mass; planets, brown dwarfs, and faint stars
overlap in radius. The giant-planet/stellar-size k=0.12 boundary is a useful
shape comparison, not a mass or classification boundary.

Radius-ratio profile diagnostics (reoptimizing all other variables):
`{'conditional_k_ranges_by_delta_penalized_chi2': {'1.0': [0.1297640593881414, 0.5], '3.84': [0.1057777893565491, 0.5], '9.0': [0.09747438681974015, 0.5]}, 'minimum_penalized_chi2': 1522.8927709132483, 'best_k_below_0_12': {'vector': [0.11957776955863655, 0.8544113485975474, 1.9541925445229418, 0.18176063906350698, 0.2878093203771764, 4.3525394466787936e-05, 0.10765731804292229], 'radius_ratio': 0.11957776955863655, 'impact_parameter': 0.9565799519484289, 'grazing': True, 'a_over_rstar': 7.058217525758671, 'inclination_degrees': 82.21089774253637, 'period_days': 3.5221960224942808, 'epoch_btjd': 1326.916339655214, 'duration_T14_hours': 2.2409317818126535, 'duration_T23_hours': 0.0, 'u1': 0.24540577025171498, 'u2': 0.18092818225911542, 'conditional_companion_radius_rsun': 0.20137015971443956, 'conditional_companion_radius_rjup': 1.959558982459052, 'data_chi2': 1522.0532462411513, 'prior_chi2': 2.31093050915797, 'penalized_chi2': 1524.3641767503093, 'n_points': 2099, 'nominal_dof': 2084, 'assumptions': {'eccentricity': 0.0, 'omega_degrees': 90.0, 'dilution_fraction': 0.0, 'companion_mass_ratio': 0.0, 'density_prior': True, 'contact_duration_approximation': 'circular exact geometry'}}, 'best_k_atleast_0_12': {'vector': [0.2778214315527998, 0.9062832213205897, 1.8702431966274653, 0.2213752573461963, 0.2867247014782047, 9.276941138616746e-05, 0.09696674105443306], 'radius_ratio': 0.2778214315527998, 'impact_parameter': 1.1580681232601588, 'grazing': True, 'a_over_rstar': 6.489874522978415, 'inclination_degrees': 79.72095880016215, 'period_days': 3.5221959155885108, 'epoch_btjd': 1326.916388899231, 'duration_T14_hours': 2.2785423126796585, 'duration_T23_hours': 0.0, 'u1': 0.26981099557180777, 'u2': 0.2006943257011788, 'conditional_companion_radius_rsun': 0.46785406894923043, 'conditional_companion_radius_rjup': 4.552748255201098, 'data_chi2': 1522.629735111469, 'prior_chi2': 0.2630358017793242, 'penalized_chi2': 1522.8927709132483, 'n_points': 2099, 'nominal_dof': 2084, 'assumptions': {'eccentricity': 0.0, 'omega_degrees': 90.0, 'dilution_fraction': 0.0, 'companion_mass_ratio': 0.0, 'density_prior': True, 'contact_duration_approximation': 'circular exact geometry'}}, 'high_minus_low_penalized_chi2': -1.4714058370609564, 'note': 'Conditional penalized likelihood profile; not probability of planet versus star.'}`. Approximate Δpenalized-chi²≤1 and ≤3.84 ranges are conditional
profile contours, not calibrated Bayesian credible intervals or global
false-alarm probabilities. Reported penalties include stellar-density and
limb-darkening priors. Compare data chi² and prior penalties separately in the
JSON. Eight-sector bootstrap (100 draws) is `{'draws': 100, 'converged_draws': 100, 'radius_ratio_p16_median_p84': [0.20017674925451737, 0.2749393921593299, 0.3003446335987328], 'impact_parameter_p16_median_p84': [1.070972051673143, 1.155489143694115, 1.1813725593253466], 'a_over_rstar_p16_median_p84': [6.449805392393976, 6.495891453397066, 6.5824842458737765], 'duration_T14_hours_p16_median_p84': [2.251664318452434, 2.2802083316888098, 2.31475900631092], 'period_days_p16_median_p84': [3.5221921144086408, 3.5221960061580644, 3.5222030374848066], 'epoch_btjd_p16_median_p84': [1326.9151192287077, 1326.9162728787219, 1326.91710409585], 'conditional_companion_radius_rjup_p16_median_p84': [3.2803601248708176, 4.505519358039695, 4.921843138353318]}`;
it probes sector sensitivity but only eight independent groups are available.
Bootstrap refits condition on the same stellar priors and use the best-fit
initialization; their narrow range is not a radius posterior and does not
replace the broad, upper-bound-reaching radius profile. Gaussian covariance
is particularly unreliable along the grazing degeneracy.

Alternative illustrative fits vary eccentricity (omega=90 deg), dilution,
stellar mass ratio, density-prior removal, and Gaia stellar density. Their
likelihoods and shapes remain in `fit_results.json`; no Bayes factors or
eccentricity posterior are claimed. A circular central transit would last
roughly twice the measured duration, supporting a high-impact interpretation
under circular geometry; eccentric velocity can also shorten a transit.
An unresolved luminous companion or contaminating star changes both depth and
stellar interpretation. Stellar-companion total mass changes Kepler a/R and
is explored separately without assuming a mass-radius relation.
Eccentric alternatives shown here fit the primary alone with omega=90 degrees;
they are illustrative velocity changes and do not jointly fit the observed
secondary phase. Exact transit/occultation timing and companion luminosity
must be included before interpreting an eccentric binary orbit.

An independently evaluated eclipse occurs at phase
0.47546039, about 2.07 hours earlier than exact half-phase. It was missed by
the original half-phase-only check. Raw SAP depth is
1326.5±232.2 ppm
(conditional S/N 5.71), and small-aperture depth is
1136.8±129.7 ppm
(S/N 8.77), consistent across all eight sectors.
The initial selected-phase statistic is not a blind-search false-alarm probability.
The small-e timing approximation implies e cos(omega)≈-0.03855 if the same
binary produces both eclipses; precise inference requires a joint orbit fit.
The detected secondary also requires geometric overlap. With the main
grazing primary geometry, the usual conjunction approximation requires
e sin(omega) ≲ 0.049.
The illustrative e=0.3 and e=0.6 periastron planet-like primary alternatives
predict no secondary overlap and therefore cannot describe both eclipses.
They demonstrate primary-only degeneracy, not allowed joint orbit solutions.

For a generous *illustrative* companion radius ≤2 Rjup, the sampled circular
primary profile gives a coupled Lambert reflected/thermal maximum of
193 ppm. An even more generous comparison sums
independent reflected (geometric albedo ≤1) and zero-albedo, uniform-substellar
temperature blackbody thermal extremes at the reddest assumed wavelength,
1.1 micron: 713 ppm. These extremes
need not coexist physically; this comparison grants full occultation despite
the grazing geometry. The sampled maximum is at k=0.120,
slightly below the exact illustrative 2-Rjup cutoff, and is not a likelihood
profile jointly marginalized over orbital/stellar/atmospheric uncertainties.

The observed secondary therefore strains an ordinary passive, planet-sized
interpretation under the adopted geometry. It is readily compatible with
larger, intrinsically luminous grazing companions: illustrative uniform-disc
secondary inversions require optical brightness temperature near 3650 K over
k≈0.2--0.4. At the primary best-fit k≈0.278, only ≈14% of the companion disc
is eclipsed, with companion/star total light ratio ≈0.0093. This is a plausible
faint-star geometry, although a full luminous binary fit is still needed.

The budget assumes blackbody spectra with a 0.6--1.1-micron top-hat photon
weight, rather than an actual TESS throughput/atmosphere calculation. It couples
Bond and geometric albedo using Lambert scattering, assumes steady passive
irradiation, and uses circular separation. Internal luminosity, eccentricity,
non-Lambert scattering, unusual atmospheres, target misidentification, and
blending can change it. Two Rjup is a comparison radius, not an absolute planet
limit. Eclipse brightness and radius never establish companion mass.
The measured secondary is box-averaged: unresolved ingress and exposure
averaging can make its true peak deeper. The illustrative temperature/light
ratios are approximate and should be replaced by a joint luminous eclipse fit.
Nonetheless, the repeated off-phase secondary and grazing shape favor an
**eclipsing-binary explanation**, and no credible planetary candidate is promoted.


Gaia median RV 13.65±3.68 km/s is a systemic-velocity estimate. The target has
no published amplitude in the cached table, so it cannot exclude a stellar
companion. Spectroscopy and source-localized higher-cadence transit photometry
remain needed. All fit uncertainties are conditional; no discovery claimed.

Reproduce/resume: `/workspace/openrouter-check-venv/bin/python src/fit_370247818.py`.
Checkpointed profile and bootstrap results, plots, model arrays, exact input
hashes, raw retrieval identifiers, and software versions are preserved here.
No raw datasets are deleted or downloaded; no paid calls are made.
