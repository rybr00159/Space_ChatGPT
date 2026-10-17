# Independent physical-model checks after reading the external report

17 October 2026. The external shape, star, background-source and brown-dwarf reports were reviewed at public commit `99d4f66fca99c32c3245396f367b8a73d5d13fa5`. Their source-specific stellar constraints provide information beyond our original free-density/free-flux stellar grid.

The supplied planet geometry agrees with our previous independent work: radius ratio approximately 0.119 versus 0.121–0.124, conditional radius 4.95 versus approximately 5.1 Earth radii, and duration approximately 1.65–1.70 hours. Agreement is descriptive because the reductions share TESS photons.

## New tests on our original observations

`shape_checks.py` uses our preserved 13 original QLP event windows, 1,961 samples, the same event baselines and errors, exposure-integrated batman models and fixed limb coefficients [0.125, 0.533]. No external Python code or external extracted light curve is executed. The small parameter choices prompted by the public report are explicitly encoded in our script; these are subsequent tests, not an additional blind discovery.

| Model | Light-curve χ² | Relevant result |
|---|---:|---|
| Planet-compatible, fixed nominal TIC density | 2284.27 | k 0.1219, b 0.565 |
| Planet-compatible, free circular effective density | 2284.10 | k 0.1197, density 12.99 g/cm³ |
| Large body k ≥ 0.197, fixed density | 2371.33 | Circular nominal-density geometry fits poorly |
| Large body k ≥ 0.197, free effective density | 2299.57 | Grazing solution at density 2.06 g/cm³, Δχ² +15.47 relative to the free-density small-body model |
| Specified K + M background pair at P | 2358.99 | Δχ² +74.73 relative to the fixed-density planet; flux-prior penalty is additional |
| Specified equal-K background stars at 2P | 2386.60 | Δχ² +102.34 relative to the fixed-density planet |

The last two rows use specified source-compatible stellar masses/radii, total background reference flux approximately 0.031 of target flux with a 10% Gaussian nuisance term, and circular orbits. They fit the same original QLP samples and baselines. They do not survey every binary configuration, include new multiband data or calculate FPP. Their failure nevertheless independently supports the supplied report's physical objection to an ordinary binary on the faint neighbor; matching the flux budget alone is insufficient.

The large-body radius cutoff is an illustrative model choice, not a universal planet/brown-dwarf mass boundary. Freeing its effective density substantially improves the fit. For the TIC density and ignoring companion mass, the low-density solution requires minimum apoastron eccentricity approximately 0.468. The supplied deblended-light-curve brown-dwarf models similarly require roughly e ≥ 0.5. Thus a circular-size argument alone does not exclude an eccentric grazing brown dwarf.

## Photo-eccentric sign and statistical limits

Numerical batman contact-duration checks at e=0.5 independently confirm `rho_circular = rho_star * g^3`, with `g=(1+e*sin(omega))/sqrt(1-e^2)`. A periastron transit gives approximate density ratio 5.185 versus analytic 5.196; apoastron gives 0.19237 versus 0.19245. This supports the sign correction adopted in the external final report. Older opposite-sign results must not be used as its current conclusion.

The supplied ordinary stellar/background fits include colours, parallax, reddening and source luminosity, but depend on chosen stellar relations, error floors and allowed model classes. Their best 3.1-arcsecond pre-ELM alternative is penalized rather than mathematically impossible. Eccentricity and observational systematics are not exhaustively covered for every background class.

Older subsidiary text called 0/25 V-shaped injections approximately 5σ. That wording is not adopted in the final main report and is not a valid simulation tail measurement: the one-sided 95% upper bound is approximately 0.113. The final report correctly presents finite counts and remains unvalidated. Finite injections and custom model likelihoods should not be interpreted as independently measured tiny FPPs.

The combination of our original recovery, the supplied improved source assignment, source-specific stellar constraints and the new bounded checks supports **strong, unvalidated planet candidate** as a reasonable follow-up classification. No measured mass, complete binary exclusion, calibrated global FAP or independently validated FPP follows.

The machine-readable [results](shape_checks.json) retain exact method/input hashes. Software/forward-model checks are distinct from confirmation of a planetary system.
