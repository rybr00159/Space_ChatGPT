# Candidate evidence index

## Stellar-system candidates

[Stellar research report](../STELLAR_RESEARCH_REPORT.md): strongest is [background TIC 273879360](stellar_systems/tic273879360/README.md), independently characterized as a strongly supported eccentric eclipsing-binary candidate through source-attributed TESS eclipses and earlier ground photometry. **Novelty correction, 2026-10-14:** it was already publicly neural-network flagged for eclipse-like features in Kostov2025 Table 2. That unvetted row has no period or vetted binary classification; this is not first recognition of an entirely unflagged system. [Giant-image TIC 260352037](../reports/stellar_systems/hidden_sources.md) contains a robust unresolved short-period eclipse signal whose actual source remains uncertain. Existing sky positions were already catalogued; component masses and first validated publication remain unestablished. [TIC 197764088](stellar_systems/tic197764088/README.md) is a preserved known-neighbor falsification.

No object is promoted as a credible new exoplanet candidate. These dossiers preserve investigated signals and their failures, so that attractive periodograms are not mistaken for discoveries and later research can revisit evidence.

## Broad-phase evidence

The broad survey actually searched 159,505 stars. A detailed dossier or an automated flag-free row is a screening lead, not an accepted planet. The original 187-star detailed batch and its version history are preserved; the current-property full-sample queue receives further checks. Read [the broad report](../BROAD_RESEARCH_REPORT.md), [v3 batch review](broad_s1/V3_REVIEW.md), and [all 389 exact signal rows](broad_s1/all_signal_evidence.csv) and [review summary](broad_s1/all_signal_review_summary.json) together with stellar-context and focused assessments.

| Target | Strongest evidence | Present interpretation |
|---|---|---|
| TIC 269994642 | [Independent 2018/2025 pixels](broad_s1/tic269994642_pixels/README.md), [physical fits](broad_s1/tic269994642_physical/README.md), [prior scope](../reports/tic269994642_prior_check.md) | Strong recurring eclipse source; doubled-period binary fits equally well, alternating depths favor binary; not accepted as planet |
| TIC 238202207 | [17-sector measurements](broad_s1/tic238202207_multisector/README.md), [physical fits](broad_s1/tic238202207_physical/README.md), [prior scope](../reports/tic238202207_prior_check.md) | Strong recurrence; unequal doubled-period binary and unresolved-source astrometry favor stellar explanation |
| TIC 260820871 | [Three-sector follow-up](broad_s1/tic260820871_multisector/README.md), [physical blend model](broad_s1/tic260820871_physical/README.md), [pixel source identification](broad_s1/tic260820871_pixels/README.md) | Independent pixels place eclipse phases at a known faint neighboring variable; hot-star planetary interpretation rejected |
| TIC 370247818 | [Eight sectors and event timings](broad_s1/tic370247818_multisector/README.md), [three pixel reductions](broad_s1/tic370247818_pixels/README.md), [physical fit](broad_s1/tic370247818_physical/README.md) | Real eclipse near target; off-phase secondary and grazing/light-budget evidence favor an EB |
| TIC 141980042 | [23-sector assessment](broad_s1/tic141980042_multisector/README.md) | Strong secondary and crowded astrometry favor EB/blend |
| TIC 308600649 | [Raw-variation review](broad_s1/V3_REVIEW.md), [current Gaia/TIC context](../reports/secondary_leads_20261014.md) | Strong ellipsoidal/RV variation and neighboring source; binary/blend favored |
| TIC 260352037 / 441386747 / 406280158 / 271697430 | [57-curve independent-sector QA](broad_s1/giant_survivors/README.md) | Current giant radii invalidate old small-planet estimates; variability/blends unresolved |
| TIC 49617264 | [Native SPOC follow-up](tce107/EXTRA_FOLLOWUP.md), [independent pixels](tce107/tic49617264_pixels/README.md), [2024 prior recognition](../reports/tic49617264_prior_check.md) | Recovered reported candidate TOI-2417.11; planetary nature unconfirmed and no novelty claim |
| Other sector-107 leads | [Initial assessment](tce107/ASSESSMENT.md), [additional follow-up](tce107/EXTRA_FOLLOWUP.md) | Three variability/baseline aliases and one published complex periodic variable |

Raw broad/follow-up light curves and pixels are retained outside Git; their exact URLs and hashes are in the observation manifests. Historical catalog and batch snapshots remain evidence even when a later assessment changes their interpretation.

## Historical pilot evidence

| Target | Signal investigated | Evidence / interpretation |
|---|---|---|
| TOI-1064 | Original sector peaks and revised 8.87/14.85 day leads | [Original diagnostics](known_host_rejected/REPORT.md), [trimmed diagnostics](known_host_rejected/TOI1064_TRIMMED.md); gap/processing dependence and weak native-flux support |
| TOI-4342 | Boundary peaks, 6.67/13.34 day alias, final 19.45 day joint peak | [Alias report](known_host_rejected/TOI4342_TRIMMED_ALIAS.md), [final joint local/raw test](known_host_rejected/TOI4342_JOINT.md); no credible native-flux/recurrence support |
| TIC 348839788 / 150359500 | Broad 0.272/1.033 day variability and residual harmonics | [Variability assessment](field_variability/README.md); no credible additional short-period signal |
| TIC 150299576 | Weak 11.3 day two-event alignment | [Assessment](tic150299576/REPORT.md); method dependent, weak raw flux and coherent support, regional null unexceptional |
| TIC 150360085 | Joint 13.123 day weak residual | [Assessment](tic150360085/ASSESSMENT.md); mixed event depths and local-baseline significance far below screening score |
| TOI-904 | Broad 13.46 day residual | [Assessment](toi904/ASSESSMENT.md); detrending/baseline dependent, no strong independent support |
| TOI-1803 | Broad 4.59 day residual | [Assessment](toi1803/README.md); broad curvature, local/raw models and narrower-event tests weaken it |

Raw observations remain in the external cache and exact products are identified in `data/observations.json`. Compact normalized arrays, event windows, measurement tables and plots in these directories are scientific evidence, including negative findings. Dossier labels do not imply that the physical source is conclusively identified or that further planets are excluded.
## Targeted claim validation

[TIC 71316387 / CTOI 71316387.01](tic71316387/README.md): five-sector alternating eclipses and independent pixel localization support a roughly 1.30367-day eclipsing interpretation of this existing RAVEN PC. Exact quoted depths/significance were not certified; no new planet or confirmed mass classification.

