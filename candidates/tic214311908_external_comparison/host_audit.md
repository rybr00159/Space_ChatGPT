# Source assignment and false-positive probability audit

Reviewed 17 October 2026, pinned public commit `99d4f66fca99c32c3245396f367b8a73d5d13fa5`. This reviews the supplied final report and selected saved outputs/code; it is not a rerun of the full six-sector raw-cube PRF analysis.

## Localization is materially stronger

Our initial two-sector analysis fitted a free Gaussian under several local sky/baseline treatments. It could not discriminate the target from close sources. The external work uses all six sectors, calibrated detector PRFs, proper-motion propagation, field-star WCS calibration, colour-dependent image shifts, injected-source tests and a separate implementation. S39/S104 provide information absent from our S66/S93 quicklook. Our earlier inconclusive result does not refute this stronger method.

The retained `referee_loc/out/fit_base.json` and its README report alternative-minus-target Δχ² approximately +126 for pair A, +146 for pair B, +285 for the Mira and +35.8 for the faint 3.1-arcsecond source. The latter nevertheless lies inside the 95% free-position region in the M-dwarf-colour coordinate frame. Its own bluer colour-dependent source hypothesis differs from that free-position hypothesis. These two statements are therefore not inherently contradictory.

Tests with 0/1,000 injected alternatives reaching the observed statistic support rejection under the adopted model, but do not directly measure tiny posterior host probabilities. Zero independent exceedances in 1,000 trials has a one-sided 95% binomial upper bound 0.00299. The Gaussian/Student-t likelihood ratios used in later Bayesian work extrapolate beyond sampled tails. PRF mismatch, colour calibration, WCS errors and shared photons remain material assumptions. The final report acknowledges these and does not claim validation.

## The 3.1-arcsecond neighbor was explicitly investigated

The original flux-budget objection is real but incomplete. Its nominal T-band flux is 2.93% of the target's. A 1.58% target-normalized deficit would require a roughly 54% eclipse on it. A uniform-brightness, circular small-angle check for an ordinary dark occultor requires radius ratio at least 0.735 and gives central flat/total duration ratio at most approximately 0.153. Limb darkening changes these approximate bounds; the physical fits below supply the more relevant test. A luminous compact component being occulted can avoid the simple argument, so shape alone is not an exhaustive exclusion.

The supplied `shape/scripts/s12_bg_eb.py` adds source-specific griz photometry, Gaia parallax and stellar relations. It allows E(B−V) from 0 to 0.30, 0.05-mag photometric error floors, a broad distance range and a 10% reference-flux nuisance term. It tests ordinary main-sequence binaries and several compact-star classes. The closest physical pre-ELM alternative has total Δχ² approximately +19.6, including approximately +10.1 from the light curve. A compact light-curve-only solution can fit the transit but fails the measured colours/parallax badly. The best pre-ELM solution reaches the specified mass and reddening bounds, which is a reason for caution about exhaustive exclusion, not evidence for a planet mass.

The new fits to our original QLP samples independently support rejecting two specified circular K-star binary geometries; see [the shape audit](shape_audit.md). The source-specific background models are substantially stronger than our earlier free-flux/free-density existence grid. However, the ordinary background stellar fits are circular, and not every eccentric/blended/compact-star population is exhaustively integrated. “Substantially disfavored” is better supported than “physically impossible.”

## What the FPP means

The initial `fpp/scripts/f04_bayes.py` explicitly includes the nearest background source and its physical scenario fits. It was not simply overlooked. Later `referee2_fpp/scripts/r05_fpp.py` carries a fixed aggregate background posterior of approximately `3e-5` from that earlier stage, rather than re-deriving it after every later modification. The final main report discloses this carry-over. The later photo-eccentric sign correction and the final prior sensitivity tables must take precedence over superseded phase-2 numbers.

We recomputed the final published prior-draw summary from all 1,500 CSV rows; the arithmetic matches. We did not independently regenerate the likelihoods, integrate a new population model or validate the adopted companion constraints. The quoted 0.0013–0.0033 FPP is therefore a conditional result of this framework, not a separately certified probability. Its stricter companion-host definition, wider prior sensitivity, unstable standard-tool output and absence of new imaging all matter. An FPP threshold developed for one code is not a universal validation rule for a different custom likelihood.

These qualifications do not invalidate the qualitative description **strong, unvalidated planet candidate**. The final supplied report already distinguishes that from statistical validation. The decisive practical tests remain resolved multi-band transit photometry, companion imaging and appropriate spectroscopy/RVs.

## Correction to our own catalogue audit

The original `data/tic214311908/catalogs/elbadry_binary_cone.tsv` has a main-catalogue row joining Gaia `5925373272850784768` and `5925373272850782208`: projected separation 756.78493 AU, `R=1.99306465e-5`. The final shifted-chance-alignment table is empty. Reading that empty control table as the whole response caused our incorrect no-entry statement. `scripts/compare_tic214311908.py` now checks all tables and exact string IDs. This corrects the wide-system context without establishing which component hosts the short-period signal.
