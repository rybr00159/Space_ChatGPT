# TOI-1803: falsification of the 4.59 d residual lead

Assessment date: 2026-10-13 UTC. **No credible additional planet candidate is
retained from this signal.** The global box score is sensitive to long-timescale
baseline structure and to the SPOC flux reduction. More local baseline models
leave only weak, inconsistent evidence. This is a rejection of the current
candidate interpretation, not proof that TOI-1803 has no additional planets.

TIC 144401492 / TOI-1803 has catalogued planets at 6.293287 and 12.885779 d.
Catalog provenance is preserved in `data/known_planets.csv` and `vetting.json`,
including the Zingales et al. 2025 reference. The present analysis uses their
linear catalog ephemerides; it does not reconstruct individual measured transit
timing variations. Broader masks test one possible timing-leakage explanation.

## The initial lead and its local-baseline test

Sector 22 with 2 d normalization and 0.25 d gap trims produced a BLS peak at
P = 4.590202809 d, epoch BTJD 1903.019175, duration 0.2295 d (5.51 hours),
depth about 588 ppm and diagonal-error S/N 8.31. Six apparent events were
reported. The duration reached the pilot search's maximum trial duration.
All quoted S/N values below are descriptive ranking or fitted-depth statistics;
none is a calibrated false-alarm probability or discovery significance.

Every putative event was refitted with an individual constant-plus-slope or
constant-plus-slope-plus-curvature baseline and a shared box-shaped event
indicator. Points within half a duration define the event; points between one
and three durations define the baseline. The intervening region is omitted.
Each event's uncertainty is conservatively enlarged if its residual chi-squared
exceeds the independent-error expectation. The combined uncertainty is further
enlarged for excess event-to-event depth scatter. Values may therefore differ
from early diagnostic runs using only uninflated photon errors.

| Sector 22 processing | Frozen global box S/N | Local linear event-scatter S/N | Local quadratic depth | Local quadratic event-scatter S/N |
| --- | --- | --- | --- | --- |
| PDC, native normalized flux | 22.86 | 3.42 | 278 ppm | 2.78 |
| PDC, 1 d median | 3.89 | 2.18 | 293 ppm | 2.80 |
| PDC, 2 d median | 8.31 | 2.76 | 290 ppm | 2.66 |
| PDC, 5 d median | 14.94 | 3.67 | 297 ppm | 1.85 |
| SAP, native normalized flux | 26.78 | 0.83 | 134 ppm | 1.05 |
| SAP, 1 d median | 2.60 | 0.62 | 165 ppm | 1.12 |
| SAP, 2 d median | 5.71 | 0.31 | 129 ppm | 0.88 |
| SAP, 5 d median | 12.42 | 0.44 | 203 ppm | 1.04 |

Large native-flux global box values arise from comparing the putative events
against a single baseline across a variable sector; their collapse under local
baseline fits illustrates why they are insufficient evidence. The 1 d and 2 d
PDC quadratic depths agree at roughly 290 ppm but have low significance and
remain substantially weaker in SAP. The first three apparent dips have some
agreement between reductions, while later individual SAP/PDC profiles differ.

The sector's good-data segments span BTJD 1900.15864–1912.50737 and
1915.16985–1926.49749. The predicted event at 1912.19958 lies 0.308 d before
the first segment ends; after the 0.25 d trim it lacks a usable baseline on its
later side. Requiring at least three baseline points on both sides removes that
event, leaving five: quadratic S/N remains only 2.61 for PDC/2 d and 0.73 for
SAP/2 d. Individual events, counts and gap distances are preserved in JSON.

Broadening both known-planet masks to ±0.25 d or ±0.5 d leaves PDC/2 d global
S/N near 8, but quadratic event-scatter S/N only 2.07 or 2.39. Thus the present
lead does not become convincing when plausible known-transit leakage is
suppressed. This does not establish complete removal of every timing-varying
known transit.

## Narrower events and baseline-preservation safeguards

At the original frozen epoch, quadratic fits with durations 0.08, 0.11 and
0.16 d return PDC/2 d depths −161, −51 and +8 ppm, respectively. A broad BLS
epoch could still be misplaced relative to a genuine narrower event, so a
training-only refinement searched **945 combinations**: 15 periods across
4.53–4.6265 d, 21 epochs within ±0.15 d of the initial epoch, and durations
0.065, 0.09 and 0.13 d. The highest quadratic event-scatter S/N is only 2.36,
at P ≈ 4.55068 d and duration 0.09 d. Its S/N is 1.21 in SAP/2 d and 1.11
in sector 49 PDC/2 d. These many correlated trials are an exploratory ranking,
not an independently significant detection.

To check that curvature fitting itself does not automatically erase boxes, a
nonblind preservation control multiplies original native flux by a 1000 ppm
box at the original investigated ephemeris before normalization. The median
per-event recovered depth increments are:

| Processing | PDC | SAP |
| --- | --- | --- |
| Native | 999 ppm | 999 ppm |
| 1 d median | 1035 ppm | 1027 ppm |
| 2 d median | 976 ppm | 984 ppm |

The small over-recoveries reflect the interaction between median normalization
and the existing variability. This control supports the use of local curvature
fits to examine this particular lead. It is not blind injection recovery,
a sensitivity curve, or an exclusion limit for realistic limb-darkened planets.

## Sector 49 and period ambiguity

Sector 49 is separated from sector 22 by roughly 740 d. The sector-22/2 d
period profile has trials between 4.5300 and 4.6265 d exceeding 95% of its
maximum BLS ranking score. This is a descriptive plateau, **not a confidence
interval**. It spans more than three cycles of accumulated phase uncertainty
over the separation. Therefore the sector-49 null at a single extrapolated
epoch is not sufficient to falsify the lead by itself.

A 4.3–4.9 d diagnostic search allows phase and duration to vary in sector 49.
It finds moderate global box maxima, but with different periods, phases and
reduction dependence. A joint search across both sectors explicitly tests
possible cycle aliases rather than assuming exact knowledge of the period.
Parameters selected jointly are not independent holdout confirmation:

| Joint processing | Best period | Global diagonal S/N | Local quadratic S/N, sector 22 | Local quadratic S/N, sector 49 |
| --- | --- | --- | --- | --- |
| PDC, 1 d median | 4.768500 d | 5.58 | 3.33 | 2.10 |
| PDC, 2 d median | 4.574375 d | 9.46 | 2.07 | 0.47 |
| SAP, 1 d median | 4.756375 d | 6.55 | 1.70 | 1.07 |
| SAP, 2 d median | 4.660750 d | 9.93 | 0.71 | 1.53 |

The attractive joint scores disappear when baseline curvature and event scatter
are considered. Different reductions also choose inconsistent periods. No
coherent localized transit interpretation is supported by these diagnostics.

## Evidence, reproducibility, and limits

`measurements.csv` contains 16 sector/reduction/normalization rows.
`vetting.json` preserves individual-event measurements, native observation
identifiers and hashes, all mask variants, narrowed duration tests, period
profile spans, joint consistency checks, refinement trials, and injection
controls. `sector22_individual_events.png` shows all six event neighborhoods;
`frozen_phase_comparison.png` shows both sectors with alternate normalization
and SAP/PDC measurements. The source is `src/vet_toi1803.py`.

```bash
python src/vet_toi1803.py
```

The current interpreter is `/workspace/openrouter-check-venv/bin/python`.
The command makes no network or paid model requests. Original FITS observations
remain in the external data cache; prepared diagnostic arrays are retained as
`$TESS_DATA_DIR/intermediate/toi1803_vetting.npz`, defaulting to
`/workspace/tess-research-data/intermediate/`. Product URIs, retrieval URLs,
SPOC processing versions and hashes are preserved through the retrieval
manifest and the JSON report.

This investigation does not yield an upper limit on additional-planet radii,
does not reconstruct the known planets' timing variations, and does not establish
the exact origin of the broad variability. No pixel localization or physical
planetary model is warranted for this weak and reduction-dependent lead.
Future work on this system should use measured timing-varying masks, additional
sectors and detrending that protects candidate transit intervals before making
a planetary interpretation. The useful current research decision is to reject
this peak as a credible candidate and prioritize less ambiguous targets.
