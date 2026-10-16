# Brighter/cooler batch and final pooled-search follow-up

16 October 2026 UTC. The prospective Sectors 96/105 batch selected 963 observations of 850 distinct TICs using tighter pre-flux brightness, temperature and radius limits. All 963 observations were searched: 274 with producer errors and 689 in the separately labelled empirical-error branch. All 113 repeated stars were blindly pooled. All 71 TICs shared with the first batch also received an all-observation pooled search. These searches overlap the survey stars and do not increase the unique-star count.

**No credible new planet or previously unknown star is established.** [TIC 229146531](../../candidates/multisector/tic229146531/README.md) is a real recurrent target-area eclipse/transit signal whose short-period companion remains unclassified. Its already published 321-day astrometric companion does not identify the short-period source. Dark-companion and doubled-period grazing equal-star shapes fit essentially equally well, with severe radius/density/dilution degeneracy. Its single-sector odd/even warnings reverse sign between sectors and do not establish a binary classification either.

## Recognition and first-pass outcomes

| TIC | Investigated signal | Outcome |
|---|---|---|
| 229146531 | approximately 0.935024 d | Six-sector recurrence, stable native source-area deficit; unresolved planet/stellar/blend alternatives. Highest follow-up priority. |
| 71274048 | 1.359244 d | QLP raw local S/N 29.09, but small aperture −0.60 and large aperture −16.98. Native pixels show no convincing deficit: all maximum pixel scores <2.5, several bound-width fits. Reference+gradient fits are consistent with zero loss. Source assignment is unsupported. Older fixed-ephemeris nonmatches alone are not a timing-uncertainty-aware veto. |
| 147271811 | 0.57445 d | Twice an already public SPOC approximately 0.28746-day signal with opposing epochs; known binary-like harmonic, not a novel planet. |
| 139526205 | 2.922 d | Known Exoplanet Archive period. |
| 229112973 | approximately 1.048 d / 0.5245-d pooled harmonic | Published ROT period approximately 1.051 d; prior variable recognition. |
| 271640189 | 2.45996 d | Matches a known ASAS-SN EA variable at approximately 36 arcsec, P2.4603 d; contamination hypothesis rather than an unreported planet. |
| 79004991 | 1.1867 d | Twice a known nearby 0.59322-day eclipsing/SB1 source; the neighbor has independent VSX/Gaia recognition. |
| 80427865 | 2.1655 / 6.5012 d aliases | Already a Gaia astrometric acceleration binary; pooled 6.5-d local raw event scores only 0.77/0.77. Variability and unstable event depths explain the attractive pooled score; no additional companion established. |

Recognition records and archived responses are in `results/multisector_recognition/b02_096_105/`. A nearby known variable is a contamination hypothesis unless ephemerides/pixels assign the flux deficit. Variable classifications alone do not rule out planets.

## Additional final pooled queues

After all pools completed, a recorded selection chose the top eight distinct retained stars per final queue with at least two supported sectors and depth <5000 ppm. The already characterized TIC229146531/238182368 stayed in their own dossiers. **Fourteen additional stars /30 observations** received raw, alternate-aperture, local-baseline, secondary and variability checks. See [the searchable table](../../results/multisector_final_pool_followup/observations.csv), source JSON and immutable selection in that directory. Scores remain conditional on the pooled ephemeris; these are not withheld observations or calibrated false-alarm tests.

TIC92626695's pooled S/N24.91 drops to corrected local 3.25/4.81, with raw −1.42/15.63 and large orbital variability. TIC197951810's pooled P5.0947d is not consistently present across three sectors; raw local 7.13/8.77/−1.34 and strong baseline/aperture differences do not support a clean common transit. Other long-period peaks mostly depend on one event per sector and weaken in raw local fits. Sparse events themselves are not a physical veto; their periods are alias-ambiguous and require pixel/coverage checks.

Two cases received further investigation:

- **TIC307085108:** broad approximately 1.54872-day troughs recur in S96/S98. An independent Lomb–Scargle search outside event windows finds 0.38647/0.38682-day modulation, making the BLS period approximately four times that waveform. The earlier 1.16127-day trial is approximately three times the same period. Two-harmonic subtraction lowers raw residual event scores to 2.95/3.68, with secondaries and aperture dependence. S90 has a corresponding broad signal; S89 does not at the fixed ephemeris. Flexible subtraction can absorb real events and the period uncertainty was not propagated over that baseline, so these observations support a variability/alias explanation without proving absence of a smaller transit.
- **TIC183373155:** nominal P15.71043d connects only one deep event in each of S96/S105. S69 does not reproduce that fixed phase, while the S29 event is in a gap (nearest good point 0.306d); the latter is not a nondetection. Native S96 differences are sky dominated and collapse to maximum S/N2.1–4.4 after sky subtraction. S105 has a bipolar difference across the image, consistent with redistribution rather than pure stellar dimming. Reference-image+gradient+background fits give −101±315 ppm in S96 and +191±651 ppm in S105 (2.5-pixel region, quadratic sky-subtracted case), compared with much deeper QLP dips. The wider region likewise agrees with zero. A Gaussian-only `localization_valid` flag missed this dipole structure: its formal offset must not be called a secure contaminant position. In the same fits, the positive target-area control TIC229146531 retains 2992±158 ppm. No common astrophysical transit source or unique period is established for TIC183373155.

Dedicated scripts preserve all queries, traceback/zero-event fallbacks and native observations. New image-gradient tests are explanatory and can also absorb blended genuine variability; they omit full PRF/pixel covariance and are not formal model probabilities. The known zero-event helper exception was retained rather than silently converted to a physical null.

The strongest remaining science is classification of TIC229146531 and better shallow-signal sensitivity, followed by further preselected sector samples. The six-sector search is a finite sample, not an exhaustive examination of TESS.
