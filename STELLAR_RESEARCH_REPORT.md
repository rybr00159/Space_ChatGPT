# Search for unrecognized stellar companions and hidden eclipsing systems

Updated 2026-10-14. **Two promising, potentially unreported eclipsing-system signals survive this phase. The stronger one is a background star, TIC 273879360, whose deep eclipses contaminated a brighter foreground target. The other is a short-period system blended with giant TIC 260352037. No newly uncatalogued star position, confirmed component masses or proven first discovery is established.**

The strategy favored unrecognized companions over new isolated sky sources: existing Gaia/TESS catalog positions do not resolve every close binary, whereas primary/secondary eclipses and source localization can reveal multiplicity. This makes effective use of the preserved astronomical data and avoids treating incomplete planet catalogs as stellar novelty evidence.

## Ranked findings

| Rank | Source | New evidence and interpretation | Remaining uncertainty |
|---:|---|---|---|
| 1 | **TIC 273879360 / Gaia 4692910966756908288** | Background source 65.59″ from the original QLP target. Both eclipse phases localize there in four usable sectors. Direct small-aperture primary/secondary depths approximately 21–27% / 7–9%; period 7.586891 d, secondary phase 0.6918. Strong luminous-companion evidence. | Already a catalogued star; no binary classification found in documented checks. Masses, unresolved close blends and exhaustive first-report status remain unconfirmed. |
| 2 | **Unresolved source in the image of TIC 260352037** | Early two-sector ephemeris predicts primary dips in 23 later sectors. Period 0.880788 d; primary/secondary raw window depths 972 ±26 / 445 ±38 ppm. Three pixel sectors place the signal near the giant image; its inferred giant radius is incompatible with an ordinary detached orbit at this period. | Exact eclipsing star is unresolved. Nearby faint TIC 737466840 is a plausible common-motion companion but is not securely selected by the pixels. |
| 3 | TIC 370247818 | Previously preserved eight-sector eclipses, an off-half-phase secondary, target-neighborhood primary localization and luminous-companion modeling. Expanded current binary catalogs still do not identify it. | Previously investigated lead, not a new detection in this phase. Grazing/source/stellar-radius degeneracies and no mass. |

Read the [background-source dossier](candidates/stellar_systems/tic273879360/README.md), [giant-image analysis](reports/stellar_systems/hidden_sources.md), and [recognition supplement](reports/stellar_systems/prior_audit_supplement.md). The source-specific distinction matters: Gaia RV variability of the foreground TIC 273879358 does not confirm the background eclipsing system.

## Actual analysis scope

The existing broad investigation had searched **159,505** real light curves and retained **319,496** period peaks. This phase reconsidered those results using stellar-eclipse criteria, selected **1,253** eligible stars and fixed a **192-star** native SAP pilot before examining its new measurements. It remeasured all 192, finding 153 robust repeated raw leads, 40 with significant conditional half-phase secondaries, and 48 with strong alternating-depth diagnostics. Counts of catalog rows and geometric archive coverage are not new analyzed-star counts. The three giant/source targets contributed 53 retained-sector investigations; seven actual QLP sectors and native pixels were used to falsify another promising blend. The background source then received four usable pixel-sector analyses and one explicitly excluded off-detector product.

The screening used multiple apertures, raw local baselines, empirical event scatter, cadence integration, primary/secondary and alternating events. Synthetic equal/unequal eclipse trains on actual timestamps checked the local estimator. These tests do not measure survey-wide completeness or false-alarm probabilities. All outcomes and failed tests remain in the [searchable raw scout table](results/stellar_systems/scout_raw_measurements.csv), with fixed selection/checkpoints and per-target JSON evidence.

## Falsification and recognition

The initial strongest new scout signals were largely already known. Eleven underwent expanded individual review: ten had a recognized binary host or nearby known EB with a matching period/harmonic; the other is an existing visual pair with an unresolved eclipsing-subsystem question. Six of eight additional deep-source checks had positive VSX/ASAS-SN recognition. Two of the apparently novel deep systems were published GALAH double-lined binaries, including TIC 38825979 and 52285789. Merely checking EB lists would have missed that prior companion recognition.

The actual [TIC 197764088 pixel test](candidates/stellar_systems/tic197764088/README.md) located both eclipses on known neighboring TIC 197764089, already listed at the same period in the **2025 TESS Ten Thousand catalog**. All 200 bootstrap realizations per eclipse phase favored that neighbor. This is a retained negative result and a validation of the catalog/source workflow.

Expanded recognition uses current VSX/live ASAS-SN, published variable and stellar-binary catalogs, actual Gaia NSS tables, SIMBAD references and selected literature aliases. Complete **7,936 + 2,065** vetted TESS EB rows and **12,760** Traven GALAH SB2 entries were checked separately from the approximately 900,000 unvetted target inventory. None of the main two candidate fields matches the vetted catalogs. TIC 273879360 is present in the unvetted inventory, with position/magnitude only; the star itself is already known. All absences are dated, bounded checks, not proof of novelty.

## Scientific and practical limits

The background source's eclipse depth and luminous secondary strongly favor a stellar companion, but a phenomenological trapezoid is not a physical mass/orbit solution. Its aperture depth changes by camera because of contamination and pixel sampling. The giant-image signal is coherent and robust to baseline choices; its component assignment is considerably weaker. TESS pixels cannot exclude a close background blend or resolve the proposed 4.21″ giant neighbor reliably. Empirical Gaussian positions omit calibrated PRF/WCS and some covariance effects.

Source positions already exist in Gaia/TIC. The potentially new finding is unresolved multiplicity/classification, not two brand-new isolated stars. Claiming newly discovered stars, confirmed masses or confirmed higher-order binding would overstate these data. Resolved ground photometry and spectra are the useful next observations.

## Preservation and continuation

All raw light curves, native pixel cubes, failed coverage products and intermediate scientific evidence remain retained outside ordinary Git. Tracked manifests record archive URLs, FFI/product IDs, processing/timing details, raw/science hashes and retrieval/check times. The large original combined catalog response (36.4 MB) also remains cached outside Git; original vetted supplements and all needed retrieval metadata are tracked. Existing planetary search sources/results were unchanged.

The [progress checkpoint](PROGRESS.md) gives exact reproduction commands. A [64-row reserve queue](results/stellar_systems/scout_next_pilot.csv) is ranked from already remeasured pilot curves after expanded known-catalog exclusions; its remaining live recognition tests are pending, so those rows are not promoted candidates. Additional parallel-agent work reached a usage limit; completed results were preserved and the primary analysis continued locally. No paid API calls were made, and no external submissions or messages were sent.

Final verification checked 335 original catalog snapshots, ten raw native pixel files including failed/pilot evidence, six current source reports, all 192 scout method identities, and source-input/attribution controls. Sixty relevant QLP files also passed recorded-hash/FITS verification. Eleven existing utility tests passed. These checks establish reproducibility and numerical/data integrity, not first-discovery status or companion masses.

The next scientific work is independent resolved eclipse photometry and mass-sensitive spectra for TIC 273879360; finer source separation for the giant blend; and current live catalog/literature checks of the reserve queue before further pixel follow-up. This checkpoint supports credible candidate systems while keeping novelty and stellar classification claims bounded by the observations.
