> Historical pilot report. The investigation has since expanded to 159,505 actual period searches. See [the broad report](BROAD_RESEARCH_REPORT.md) and [current progress](PROGRESS.md) for the latest findings.

# TESS research report — 2026-10-13 pilot

## Scientific result

**No credible, previously unrecognized exoplanet candidate has been established in this pilot.** This conclusion concerns the searched data and methods; it does not show that these stars lack additional planets. Several apparent high-score signals were investigated and failed specific falsification tests. Weak residuals are preserved as unpromoted signals rather than converted into discovery claims.

Completed screening comprises 135 per-sector variants (including three unmasked controls), 16 additional untrimmed joint variants from the incremental archive checkpoint, and 58 gap-edge-excluded variants covering the final pilot. These retain 775 ranked peaks; not all are individually vetted.

The scientific work extends beyond a pipeline demonstration: observations were retrieved and verified, known events recovered, artificial events injected before detrending, residuals searched, stronger leads inspected against raw flux and other sectors, and prior reports checked. The resulting reproducible negative findings are useful evidence about the pilot and its failure modes.

## Strategy and sample

The two initial tracks were (1) residual searches for additional planets around small known transiting systems and (2) bright cool field dwarfs outside current planet catalogs. Known systems offer established stellar properties and transit geometry, but have substantial prior scrutiny and known-event masking costs. Field stars avoid known-planet masks but require more activity and contamination screening; catalog absence is not proof of little prior analysis. The saved Reddit example concerns publicly reported TIC 4206066 signals near 3.18 and 11.13 days, which were explicitly excluded from novelty claims.

The pilot contains eight known two-planet hosts (TOI-776, 1064, 1468, 1470, 2095, 4342, 1803, 904), ten independently selected field dwarfs, and L 98-59 as a positive control: **19 targets and 66 SPOC light-curve products**. Known hosts usually contribute their latest three available sectors; TOI-1803 has two. Field selection used a 3° cone at RA 100°, Dec −60°, Tmag 9–11.5, TIC Teff 3200–4800 K, radius 0.2–0.8 solar, contamination ratio ≤0.05, then magnitude/contamination/TIC ordering and actual SPOC availability. Twelve stars passed catalog cuts and ten had adequate coverage. Four recent available sectors per field star were downloaded. Selection preceded flux analysis. This is a purposeful feasibility sample, not a statistically matched occurrence-rate experiment.

See `reports/strategy_review.md`, `reports/field_selection.md`, and `data/catalogs/pilot_targets.csv` for exact choices. Available but unexamined sectors remain recorded in query/selection caches.

## Methods and integrity

- Download exact public MAST CAOM product URIs, cache raw FITS outside Git, preserve processing/header metadata, and record SHA-256. Raw FITS header-plus-data ones-complement checksums are enforced. Astropy's verifier mishandles these products' absent DATASUM cards; the standard raw-byte sum independently verifies them without changing source files or expected checksums.
- Require finite positive flux/error and QUALITY=0; sort aligned arrays. Require TDB timestamps with an explicit BJD reference equal to 2457000. Catalog ephemerides use BJD and are converted explicitly to BTJD.
- Normalize each segment, splitting at gaps >0.2 days, with 1 and 2 day robust rolling medians. Remove high flare/outlier samples and bin at 10 minutes; leave negative excursions for inspection. Inflate diagonal errors to at least the robust observed scatter. Store intermediate arrays.
- Mask known confirmed and unmatched TOI ephemerides using conservative duration windows. Masks do not identify every possible TTV or ephemeris-error event; known TOI-1803 timing variations motivated widened linear masks in dedicated vetting; measured individual TTV ephemerides were not reconstructed. Prior RV candidates lacking transit epochs are recorded for interpretation rather than blindly masked.
- Run Astropy fast BLS with periods from 0.7 days to half the sector baseline, normally about 13 days; contiguous-sector searches extend to 25 days. Durations 1.08–5.52 hours, frequency_factor=0.7, oversample=10, objective=S/N. Retain ranked peaks and explicit observed event coverage. Separate robustness searches remove 0.25 days at each segment edge.
- For shortlisted leads, inspect individual events, raw SAP and corrected PDCSAP, fixed-epoch predictions in other sectors, local linear/quadratic baselines, event scatter and red noise. Variability vetting additionally uses Lomb–Scargle/harmonic fits, shorter-period BLS, frozen training ephemerides, 5 day detrending, exact slow BLS at shortlisted parameters, and preservation injections. Sector 97/98 heldout parameter fits were already viewed in initial screening and are not pristine blind validation.

BLS S/N and robust periodogram SDE are ranking statistics, not survey false-alarm probabilities. “Supported events” means enough in-transit and nearby out-of-transit data, not proof that each event is a real dip. Odd/even tests are unavailable when one parity lacks coverage; absence was explicitly marked rather than interpreted as agreement. Fast BLS approximations can be misleading for broad high-duty variability, so exact models were used for that comparison.

## Readiness and empirical checks

All three known L 98-59 planets are recovered at catalog ephemerides in sectors 2/5/8. The strongest free-search signal is the known 3.6906764 day planet, with recovered periods 3.6909/3.6905/3.6917 days and diagonal S/N 28.0/30.1/23.2. Local known-event measurements recover smaller signals too, but sparse coverage can reduce significance sharply.

Six blind full-period injection searches (500/1000/2000 ppm × two detrenders), using 1.8 hour boxes at 4.731 days added to raw flux, recover the injected period as rank one and retain about 97–98% of differential depth. Two additional unchanged-data searches provide controls, for eight total runs. The favorable injection epoch was selected using coverage only; this is basic sensitivity validation, not population completeness. Five utility tests exercise checksum corruption detection, timing/masking, aligned indexing, gap handling and local-baseline fits.

Ninety-six additional half-day block-permutation searches (24 per width for each of two stars) demonstrate ordinary period-grid maxima near S/N 5. Paired-width illustrative plus-one tail fractions are 0.24 for TIC 278727563 sector 94 and 0.76 for L 98-59 sector 2. A dedicated TIC 150299576 regional moving-block check gives 30/64 null maxima exceeding its joint observed value (plus-one 0.477). Fixed-flattening/exchangeability assumptions, small tail resolution, selected targets and incomplete all-target trials prevent treating any of these as global false-alarm probabilities. Full outcomes, seeds and resumption checkpoints are retained.

## Investigated leads

| Target / signal family | Evidence | Current interpretation |
|---|---|---|
| TOI-1064 leading ~12 day sector peaks | High diagonal S/N selected one-event/gap-boundary excursions; raw flux and other sectors disagree | Original periodic interpretation rejected; detailed trimmed alternatives assessed separately |
| TOI-4342 original sector-boundary peaks; trimmed aliases including final 19.45 day joint peak | Dips depend on corrected flux; gap-trim and local raw-flux checks weaken them; sampled predictions in sectors 101/103 lack common dips | No credible candidate; two-event period ambiguity and raw baseline drift acknowledged |
| TIC 348839788 ~0.8166/1.089/1.361 day BLS family | Broad fundamental ~0.27228 days, 505–520 ppm peak-to-peak; harmonics explain ranks; residual training peaks lack robust recurrence | Stellar/contaminant variability favored; no additional short-period candidate |
| TIC 150359500 ~1.033 day waveform and ~0.68977 day residual alias | Broad ~4% modulation agrees in raw flux; residual alias changes with detrending and local-baseline/event-scatter significance is ≤3.1 | Variability/harmonic residual, not a credible transit candidate |
| TIC 150299576 ~11.3 day, two central sector-13 dips | Four-sector coherent score weakens, first raw-flux event much shallower, sampled predictions lack common depth, selected-star null tail ~0.48 | Weak unpromoted signal; no credible planetary interpretation |
| TIC 150360085 ~13.123 day joint residual | BLS ~6.9 becomes local-linear S/N ~3.6 and event-scatter ~2.4; quadratic/raw fits weaker; event signs mixed | Weak exploratory residual; no promoted candidate |
| TOI-904 ~13.46 day broad residual | Strongly detrending/local-baseline dependent, raw-flux significance weak, sector 88 does not independently support it; 5.5 hour duration exceeds ~2.77 hour circular central expectation | Weak unpromoted residual; eccentricity can change duration, so duration mismatch is not proof of exclusion |

Detailed candidate directories preserve figures, measurements, intermediate evidence and limitations. Final trimmed TOI-1064 peaks also lack raw-flux support: the 14.85 day PDC score 9.63 contrasts with SAP −1.74, and its last event has negative native SAP local depths. TOI-1803's broad 4.59 day PDC score 8.31 becomes local-quadratic 2.66 and SAP 0.88; a 945-combination narrower-event safeguard reaches only event-scatter S/N 2.36. A final TOI-4342 19.45 day joint maximum has corrected-flux local significance near 9 but native SAP linear/quadratic S/N −1.37/0.62, indicating reduction dependence rather than credible repeated transits. Only two selected cycles have support: the first has two ten-minute bins, versus 31 in the second. Both are the same cycle parity, so odd/even agreement cannot be assessed. The candidate ranking contains no promoted object unless a later report documents a survivor. Remaining lower-ranking peaks are screening entries, not individually verified planetary exclusions.

## Prior recognition and contamination

Dated snapshots contain 6,445 confirmed planets (including RV-only), 8,148 TOIs and 5,138 CTOIs. Twenty transit ephemerides are used for pilot/control masks. Cross-checks include exact TIC and position matches; selected-host abstracts, arXiv metadata, SIMBAD references, ExoFOP, Gaia DR3 and TESS EB object/detail pages provide additional context. Missing or unavailable entries are not treated as proof of novelty.

The current TOI table already includes TOI-2095.03 at 57.6869 days. Published papers report a ~47.5 day RV candidate on TOI-4342 and a ~1.7361 day RV candidate on L 98-59 absent from the queried confirmed catalog. These periods are prior recognition, not discoveries here. Literature coverage is explicitly limited: abstracts/host aliases can miss supplemental and unindexed reports.

Low TIC contamination and high CROWDSAP do not rule out blends. TIC 150299576 has a 7.83″ neighbor with Tmag 16.813 versus 11.223; its ~0.0058 flux ratio could dilute a 10–20% eclipse into the weak observed depths if admitted to the aperture. Other field stars have SIMBAD companions/activity flags, Gaia solar-like variability flags or nearby CTOIs. No source was positively classified as an eclipsing binary just because the EB website returned an object: the relevant detail pages say `In catalog? False`.

No signal has survived sufficiently to warrant a full pixel-localization/physical-transit validation campaign. Pixel files, centroid/difference imaging, secondary-eclipse tests and a population-wide false-positive analysis therefore remain incomplete. This precludes any planetary confirmation or discovery claim.

## Reproducibility, storage, and next investigation

Raw observations and intermediates are retained at `/workspace/tess-research-data`; ordinary Git stores code, compact evidence, plots, searchable `results/targets.csv` and `results/signals.csv`, exact product manifests and catalog snapshots. `python src/retrieve_manifest.py` restores missing exact products and verifies original hashes. Nothing scientifically important was deleted to reduce size. Published cloud snapshots can retain prepared files, but separate durable storage and new-task restoration were not verified; processes are not persistent.

The current negative pilot cannot establish a yield difference between tracks. A justified next phase is to prioritize quiet field dwarfs, improve sub-hour sensitivity and gap-aware detrending, and add heldout archival sectors before a much wider search. Additional-planet searches remain scientifically worthwhile with explicit known-event/TTV masks and independent pixel reductions. Longer periods/single transits and weak radius/period completeness need separate methods and representative injection grids. Review `PROGRESS.md` before resumption and reuse caches/checkpoints.

Research-phase OpenRouter calls/spending: zero. No other paid service used. Installation and research startup instructions are saved in the environment draft; saving is not publication. Review/save and publish through environment settings to snapshot the prepared machine.
