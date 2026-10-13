# Recurring field signals: independent waveform and residual vetting

Assessment date: 2026-10-13 UTC. These are **rejected transit-search interpretations**,
not confirmed classifications of the stars. No credible additional short-period
transit candidate survives this pilot screen. The photometric source has not been
localized with pixels, and a weak planet could remain below this search's sensitivity.

## Main findings

| TIC | Period selected on sectors 95–96 | Sinusoidal peak-to-peak amplitude, sectors 95–98 | Sinusoid fraction of total weighted variance | Best periodic box fraction |
| --- | --- | --- | --- | --- |
| 348839788 | 0.272276 d | 505–520 ppm | 4.4–7.3% | 3.2–5.6% |
| 150359500 | 1.033288 d | 37,537–39,811 ppm | 83.9–94.0% | 66.1–74.8% |

The quoted period precision describes the numerical fit, not an established
physical uncertainty. The waveforms vary, and the baseline contains correlated
noise. The low total explained fraction for TIC 348839788 reflects substantial
other variability/noise; the periodic folded profile itself is smooth and broad.
The original 0.8166, 1.089 and 1.361 d BLS peaks are approximately 3, 4 and 5 times
the 0.2723 d fundamental. Their recurrence is explained by a repetitive smooth
waveform that BLS can represent with broad boxes. The 1.03 d TIC 150359500 peak
similarly follows a nearly sinusoidal, approximately 4% modulation.

Both signals recur in raw `SAP_FLUX` as well as `PDCSAP_FLUX`. Raw-SAP fitted
amplitudes are 486–534 ppm and 37,086–39,458 ppm respectively. A 2 d rather than
5 d normalization window preserves the same approximate fundamentals. This
supports a real aperture variation, while leaving stellar rotation, other stellar
variability, blended sources and instrumental contributions open. It does not
identify a planet or establish which star supplies the flux variation.

No known eclipsing-binary classification is asserted. The other catalog review
found that a TIC detail page can exist with **“In catalog? False”**; page existence
alone is not a positive catalog match. Catalog and literature novelty checks are
recorded separately in `data/catalogs/`, and cannot make these broad shapes into
planetary-transit candidates.

## What was tested

`src/vet_field_variability.py` operates on immutable cached FITS products listed
in `data/observations.json`, selecting only these two TICs and sectors 95–98.
Actual timestamps are retained, so gaps do not become artificial uniformly
spaced samples. `QUALITY == 0`, finite positive flux and uncertainty cuts are
applied. A gap-separated running median normalizes 2 d or 5 d variability;
weighted 5-minute bins and 0.25 d trims on either side of gaps define the pilot.
Positive-only 7-MAD clipping limits large flare leverage; negative points remain.
Photon errors are propagated through binning. A robust residual-scatter floor is
applied **after** Fourier removal, rather than inflated by the star's original
large modulation.

Lomb–Scargle frequencies are computed with the actual observation times.
Sectors 95–96 choose the frequency and 1–4 harmonic complexity, using weighted
least squares and diagonal-error BIC as a descriptive complexity criterion.
For TIC 348839788 the chosen frequency is restricted to the previously noticed
subday family, rather than its stronger long-timescale noise. Sectors 97–98 had
already appeared in the initial pilot screen; they are held out from these
parameter fits, but are not a pristine blind discovery validation set.
Both the frozen training waveform and a separately amplitude-refitted waveform
are examined in sectors 97–98. The latter is exploratory, and can absorb an
additional signal exactly coincident with a harmonic.

The residual BLS search covers 0.1–0.7 d with durations 0.01, 0.02, 0.035 and
0.05 d. All durations are shorter than the minimum trial period. A logarithmic
grid limits period mismatch to approximately one third of the shortest trial
duration over the observed baseline. Fast BLS ranks trials; all shortlisted peaks
are reevaluated at the same period with `method='slow'`. The broad-shape box
comparison also uses the slow method. This distinction matters: the fast method
badly underfit a 45% duty-cycle box in one diagnostic, reporting a 23% explained
fraction where exact evaluation recovered 74%. That misleading comparison was
corrected before the final results were saved.

The top three residual ephemerides from **each** training sector were propagated
unchanged to both held-out sectors. Across two TICs and three preprocessing
variants, none of the 36 training-selected ephemerides has diagonal S/N ≥ 5 in
both held-out sectors. Five is an explicitly uncalibrated screening threshold;
these S/N values are not discovery probabilities or global false-alarm rates.
Many periods and durations were searched, so isolated S/N maxima are expected.
The searchable measurements are in `frozen_ephemeris_checks.csv`. This rejection
uses independent recurrence and waveform shape, rather than merely rejecting
the highest BLS peak.

## A residual lead investigated further

An early 2 d-normalization run, before correcting the positive-flare clipping,
showed a TIC 150359500 residual near **P = 0.689773 d**, approximately two thirds
of the strong 1.033 d waveform. Some held-out global box S/N values remained
large enough to warrant a dedicated falsification experiment even though the
peak later fell below the top three training trials.

Sector 96 with 2 d normalization selected a fixed epoch of BTJD 3910.201273
and duration 0.050 d. Every other sector and normalization evaluated exactly
that period, epoch and duration:

| Processing | Sector 95 global box S/N | Sector 96 | Sector 97 | Sector 98 |
| --- | --- | --- | --- | --- |
| PDC, 2 d median | −0.1 | 8.3 | 6.8 | 7.2 |
| PDC, 5 d median | 1.2 | 1.7 | 4.1 | 2.3 |
| SAP, 5 d median | 1.4 | 1.1 | 4.1 | 1.8 |

The folded residual is a broad positive/negative wave, rather than a localized
flat-bottomed dip. A separate weighted linear baseline around each putative
event, with uncertainty inflated by observed event scatter, gives S/N only
1.6, 3.1 and 2.0 for sectors 96–98 with the 2 d normalization; the 5 d values
are 0.8, 1.6 and 0.05. Robust 60-minute scatter exceeds the independent-error
expectation by approximately 3.2–3.7. Thus the attractive global box statistics
measure correlated residual waveform structure and are not reliable transit
significance.

A nonblind preservation control multiplies the original native flux by a
1000 ppm box at this same fixed ephemeris **before detrending**, then repeats
normalization and harmonic fitting. The recovered measurement difference is
785–849 ppm with the 2 d median and 790–865 ppm with the 5 d median. This
shows appreciable attenuation in both variants, while demonstrating that the
5 d treatment does not remove all ordinary box depth at this ephemeris. The
observed signal's much larger sensitivity to normalization is therefore a
warning about its broad waveform origin. This small fixed-phase control is
not a blind recovery experiment, a completeness curve, or a measured false
positive rate.

Disposition: **detrending-sensitive commensurate residual, rejected as a
credible transit candidate**. More flexible stellar-variability models or a
longer held-out baseline could improve sensitivity, but the present evidence
does not justify planet-candidate language.

## Reproduce and resume

With the project's Python dependencies installed:

```bash
python src/vet_field_variability.py
python src/vet_field_variability.py --width 2
python src/vet_field_variability.py --column SAP_FLUX
python src/vet_field_variability.py --probe-only
```

The current environment's interpreter is
`/workspace/openrouter-check-venv/bin/python`. No network request or paid model
call is made by these commands. Metadata, product identifiers and SHA-256 values
are embedded in every per-target JSON report. Raw observations remain in the
external cache; prepared residual arrays are retained under
`$TESS_DATA_DIR/intermediate/tic*_FLUX_w*_variability_residuals.npz`, defaulting
to `/workspace/tess-research-data/intermediate/`. The project retrieval manifest
reconstructs original observations if this local cache does not persist.

`waveform_summary.csv` has 24 target/sector/variant rows;
`frozen_ephemeris_checks.csv` has 36 residual-ephemeris rows.
JSON files preserve all shortlisted peaks, Fourier coefficients, alternative
frozen-waveform searches, individual-event measurements and injection controls.
PNG files show period spectra, full-phase waveform shapes and residual folds.

Useful next investigation: search less variable field targets, where weak
transits require fewer stellar-model assumptions. If revisiting these objects,
use an older independent sector, compare a controlled stellar-variability
model against protected-transit detrending, and perform multiple-phase blind
injection-recovery plus pixel/source diagnostics for any localized surviving
event. The present analysis makes no completeness statement for periods
outside 0.1–0.7 d, gap-adjacent transits, or depths below its residual noise.
