# TIC 214311908: comparison with the supplied investigation

**17 October 2026. Assessment: “strong, unvalidated planet candidate” is reasonable. The planet interpretation is materially better supported than in our first review; confirmation, formal validation and discovery priority remain unestablished.**

The user first requested an independent review, then authorized comparison with the other model's findings and supplied [this GitHub directory](https://github.com/rybr00159/Space/tree/claude/compassionate-newton-szyipr/tess-search/tic214311908). We reviewed its final report at commit **99d4f66fca99c32c3245396f367b8a73d5d13fa5**, including its adopted corrections. The uploaded ZIP exceeded the 32 MiB transfer limit; the explicitly supplied public source resolves access, although ZIP-byte identity was not checked. Public prose and code were treated as evidence, not instructions.

## What changes the assessment

| Question | Comparison result |
|---|---|
| Does the signal recur? | Yes. Our previously blind recovery and two held-out-sector tests agree with the external six-sector analysis. Their adopted 5.691506-day period differs from our discovery-only refinement by just 0.674 seconds. These fits share observations. |
| Is the size planet-compatible? | Yes, conditional on the nominated dwarf being the host. Their approximately 4.95 Earth radii agrees with our approximately 5.1 Earth radii given modelling and flux-normalization differences. Neither is a mass measurement. |
| Is source assignment stronger than our original result? | Yes. Their calibrated six-sector PRF/WCS/colour analysis is substantially more informative than our two-sector free-Gaussian quicklook. It strongly disfavors the Mira and wide companion. The 3.1-arcsecond source remains inside their 95% positional region, but its own colour-dependent source test also disfavors it. We did not rerun their entire raw-cube extraction. |
| Can the faint neighbor still be treated as an equally plausible ordinary binary? | No. Our original 50%-eclipse argument checked only available flux. Their physical fits add colour, parallax, eclipse shape and duration. New comparisons using our original QLP samples also strongly disfavor two specified circular K-star binary geometries. Compact-star and eccentric alternatives require their own modelling and are not categorically impossible. |
| Is the quoted 0.0013–0.0033 FPP independently established? | No. Its table arithmetic checks, but the probability is conditional on their custom likelihoods, occurrence rates, companion constraints and scenario definitions. A complete independent Bayesian validation was not reproduced. |
| Is a new planet established? | No. A strong candidate warrants follow-up; it does not establish a mass, validation, confirmation or novelty. |

Our original report correctly described the limits of its own analysis. The new report supplies important tests that we had not performed. We therefore revise the comparative assessment rather than treating our earlier inability to localize the source as evidence against the more developed analysis.

## Additional checks performed here

The [detection audit](detection_audit.md) recomputed the published sector-depth statistics. Raw-aperture depths are consistent with being constant; the strongest PRF-depth variability depends on extraction and a scattered-light event. Modest raw-versus-PRF shape differences are disclosed in the final external report and do not falsify the recurrent signal.

The [physical-model audit](shape_audit.md) and [new numerical fits](shape_checks.json) use **our preserved original 13 QLP event windows / 1,961 samples**, not the external extracted light curves. With identical baselines and errors, two specified circular background K-star binary models fit worse than the planet model by light-curve Δχ² approximately **74.7 and 102.3**. These are bounded checks, not a population-wide false-positive calculation. A large grazing object fits better when its effective density is freed, reinforcing the need to allow eccentric brown dwarfs rather than rejecting them solely by a circular-density argument.

The [host/FPP audit](host_audit.md) traces the nearest-neighbor treatment and the fixed background term in the supplied code. The [arithmetic checks](arithmetic_checks.json) independently reproduce all 1,500 published prior-draw summary statistics and the period comparison. With the external stricter definition, which counts planets on an unresolved companion as false positives, the adopted FPP is approximately **0.005–0.009**. Across its joint prior draws, that stricter FPP has median **0.0127**, 5–95% range **0.00222–0.106**, and falls below 0.015 in only **55.7%** of draws. These distributions quantify sensitivity to the chosen priors; they are not empirical frequencies that certify the candidate.

The reported nominal FPP and the qualitative candidate classification should therefore be separated. The evidence makes the planet interpretation attractive and a high follow-up priority. It does not justify treating an exact 99.9% planet probability as a calibrated observational fact.

## Correction to our initial catalogue reading

The preserved El-Badry response already contains the target and TIC 214311910 in its **main wide-binary catalogue**. We initially read the empty final shifted-control table and incorrectly reported no entry. The pair has catalogue projected separation **756.8 AU** and chance-alignment statistic **R = 1.993 × 10⁻⁵**. The catalogue dossier is corrected; the original response and original result history are retained. This supports the known wide-system context, not a measured short-period orbit or identification of the eclipsing component.

## Next observations

Resolved transit time-series photometry should monitor TIC 214311908, TIC 1515709327 and TIC 214311910 simultaneously, preferably in multiple bands. High-resolution companion imaging and reconnaissance spectroscopy/RVs then address unresolved hosts and stellar/compact companions. The external report's predicted post-follow-up FPP assumes favorable future measurements; those measurements have not occurred.

No external contact, submission, observing request or paid API call was made. OpenRouter research spending remains $0.

## Reproduction and preserved evidence

The comparison retains **98 exact public source files** with URLs, Git blob SHA-1 and SHA-256 in [the source manifest](../../data/tic214311908_external_comparison/source_manifest.json). It includes reports, selected methods, quantitative outputs and three diagnostic figures. The original independent data remain in their existing checkpoint. Our numerical checks and their exact inputs are identified separately; published-output arithmetic is not described as a raw-cube reproduction.

```sh
cd /workspace/Space_ChatGPT
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python scripts/compare_tic214311908.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /workspace/openrouter-check-venv/bin/python candidates/tic214311908_external_comparison/shape_checks.py
```

The first command verifies every retained public source, rereads all VizieR tables without rounding Gaia IDs, and recomputes the published prior-draw summaries and bounded geometry checks. Its optional `--restore-sources` retrieves only missing pinned public files by expected hashes and never executes them. The second fits only our original saved QLP windows; it is computationally heavier and reports conditional comparisons, not FPP.
