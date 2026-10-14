# Independent alternate-peak audit

Reviewed all 32 alternate-period rows and their existing main-period dossiers. No credible new planet is established. These rows belong to already counted stars; most represent aliases of existing dips. The next scientific work is independent-sector and source validation, with the rankings below.

Current TIC estimates retrieved on 2026-10-14 are used throughout; old QLP radius headers were not used. Frozen v3 code and numeric dossiers are unchanged. Single-sector S/N is conditional and lacks survey-wide false-alarm calibration.

## Ranked next checks

| Priority | TIC | Preferred period (d) | Status |
|---|---:|---:|---|
| 1 | 260820871 | 0.673095 | High priority within this subset |
| 2 | 290028406 | 6.772758 | High priority within this subset |
| 3 | 262608189 | 2.134621 | Revisit existing main-period signal |
| 4 | 232074299 | 0.966920 | Revisit existing main-period signal |
| 5 | 140803419 | 4.579859 | Revisit coverage-limited existing main-period signal |
| 6 | 89330579 | 4.073798 | Lower reserve |
| 7 | 260131114 | 1.356190 | Lower reserve existing main-period signal |

### TIC 260820871

28/28 positive raw events, SNR 13.25, depth 438 ppm, broadly consistent across apertures. Current R=2.40763 Rsun, M=2.85 Msun, Teff=10942 K. True-period half-phase secondary is -0.58 sigma and selected other phase 3.19 sigma.

Odd/even 3.21 sigma; grazing/long ingress. QLP x centroid is 12.14 sigma but only 0.000684 pixel. Stellar-companion, gravity-darkening/rotation and blend alternatives remain. The 1.347 d secondary alone does not reject a 0.673 d planet.

Retrieve independent sectors, common-ephemeris parity depths/shapes, Gaia/SIMBAD/VSX/literature and source pixels; compare planetary and binary models.

### TIC 290028406

The longer alias restores three consistent raw dips, SNR 15.44; the shorter alias mixes strong and shallow events and falls to 4.53. Current R=1.08849 Rsun, M=1.02 Msun, planetary depth scale possible.

Only three events, ingress fraction 0.453; aperture and centroid warnings. Half-phase and selected other-phase tests each sample only one event, so secondary absence remains weakly constrained.

Independent sectors, additional eclipses, difference-image source and aperture dilution model.

### TIC 262608189

Eight raw events, SNR 14.89, odd/even 0.50 sigma; current R=0.919909 Rsun and fitted undiluted size near 0.86 Rjup.

Small/raw/large depths 2273/7762/11708 ppm and x centroid 18.74 sigma; true-period secondary 2.93 sigma raw and 4.09 sigma large aperture. Double-period secondary is expected from any repeating short-period dip.

Independent-sector secondary tests and source pixels; centroid/aperture warnings alone do not falsify it.

### TIC 232074299

18 raw events, SNR 19.80, odd/even 0.79 sigma; true-P secondary 1.45 sigma half phase and 2.71 sigma selected phase. Current R=1.25167 Rsun gives fitted undiluted size around 1.50 Rjup.

Small/raw/large depths 4443/12017/17892 ppm; y centroid -0.01542 pixel at 32.3 sigma. A blended stellar source is an important unresolved alternative.

Independent sectors, source pixels, full-orbit and parity-shape comparison of short-period planet vs equal eclipses at twice the period.

### TIC 140803419

The 9.158 d alias is twice the 4.57986 d repeating dips. At least four deep dips are visible in raw flux, but frozen local fits retain three. Current R=0.83794 Rsun, M=0.92 Msun.

Zero recorded secondary at 9.158 d is untested coverage: the early sampled deep dip is excluded by a quality gap 0.062499881 d barely above three median cadences; two later secondary windows have no good points. Coherent raw 4.58 d modulation is 2371 +/-126 ppm, 18.8 sigma, favoring stellar variability/blend/ellipsoidal alternatives.

Independent sectors, documented gap-tolerance sensitivity kept separate from frozen v3, source attribution and full-orbit models.

### TIC 89330579

Five consistent primary dips, raw SNR 11.60, odd/even 0.29 sigma, current R=1.20506 Rsun.

Selected other-phase eclipse is 4.97 sigma raw, 6.40 sigma large aperture, from only two events; source warnings. It is not presently secondary-free.

Independent-sector all-phase secondary/source tests.

### TIC 260131114

Longer original period has 14 raw events, SNR 12.12 and no strong secondary; current R=1.20593 Rsun.

The 0.67684 d alternate has excessive duty cycle and odd/even differences. Original longer solution still has depth-aperture dependence and x centroid 21.2 sigma.

Independent sectors at original longer period and pixel source checks if higher priorities fail.

## Shared short-period signal

TIC 139199396 and TIC 152819299 are 1.971 degrees apart but share camera 2 / CCD 3, periods near 0.7012 d differing by 0.0000504 d, and epochs within 17 minutes. Across 1043 shared cadences their processed flux correlation is 0.8218, remaining 0.8194 after a 97-cadence median high-pass. Raw high-pass correlation is 0.6966. TIC 152819299 also has raw/background correlation 0.8487 and raw MAD 0.486 in normalized flux. Shared instrumental or optical contamination is favored; its physical origin is unidentified. These are not two independently supported planet candidates.

## Aliases and incomplete secondary tests

A secondary caused by folding repeated dips at twice their photometric period also appears for a true short-period planet. It does not independently reject that shorter period. Compare parity depths/shapes at a common ephemeris, off-phase secondary recurrences at the shorter period, full-orbit variability, physical geometry and source position. Centroid and aperture warnings warrant tests rather than automatic false-positive claims.

TIC 140803419 has a sampled deep dip at the secondary phase of its 9.158 d alias. Frozen local fitting excludes it because a 0.062499881 d gap barely exceeds three median cadences; two other secondary windows lack good data. Its stored zero secondary score therefore represents missing usable measurements. Strong raw 4.58 d modulation independently favors a stellar alternative.

TIC 183642107 also cannot be interpreted by the doubled-period secondary alone: the short-period single trapezoid has chi-square 1682 for 279 points, while the longer alias isolates deeper narrow events. A common physical eclipse/shape comparison is required before promoting it.

## Evidence and reproducibility

`independent_audit.json` retains all 32 rows, current stellar catalog fields, sibling dossier hashes, period-specific raw depths/event counts, exact common-mode correlations, input FITS hashes/paths and the secondary-coverage edge case. The cross-target check uses quality-zero finite native cadences, aligns CADENCENO, and computes Pearson correlations before/after `scipy.ndimage.median_filter(size=97)` subtraction. All raw inputs remain in their original cache. No paid calls, new downloads or commits were made by this audit.

## All alternate rows

| TIC | Extra period (d) | Raw S/N | Events | Manual status |
|---:|---:|---:|---:|---|
| 260820871 | 0.673095 | 13.25 | 28 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 260131114 | 0.676836 | 9.61 | 24 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 309791057 | 1.482779 | 11.50 | 9 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 141122202 | 2.699531 | 8.76 | 7 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 140688841 | 3.605312 | 9.50 | 5 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 89330579 | 4.073798 | 11.60 | 5 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 167808509 | 1.864116 | 12.58 | 8 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 232074299 | 1.934296 | 14.19 | 10 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 278730067 | 2.491329 | 10.01 | 7 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 355652280 | 2.768703 | 6.55 | 5 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 262608189 | 4.266461 | 12.33 | 5 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 141980042 | 4.641977 | 20.97 | 6 | Independent 23-sector followup finds a robust secondary and strong blend context; stellar eclipses favored. |
| 126874570 | 2.663596 | 1.30 | 3 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 152819299 | 0.701202 | 7.04 | 29 | Shared 0.7012 d signal with TIC 139199396, extreme raw scatter and raw/background correlation; do not promote as an independent planet. |
| 300563153 | 6.134979 | 9.93 | 1 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 300384329 | 9.449327 | 9.80 | 3 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 370247818 | 3.522500 | 14.35 | 5 | Independent followup finds displaced phase ~0.475 secondary and stellar/source alternatives favored. |
| 139199396 | 0.701252 | 11.45 | 29 | Shared nearly identical 0.7012 d signal with TIC 152819299; same CCD and strong high-pass correlation across 1.971 degrees favor instrumental/optical contamination. |
| 290028406 | 6.772758 | 15.44 | 3 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 238202207 | 4.643371 | 15.07 | 5 | Primary followup already active; double-period alias, not an additional planet. |
| 261107520 | 7.129890 | 11.94 | 4 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 140803419 | 9.158127 | 19.05 | 3 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 97408876 | 2.167238 | 5.04 | 6 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 39163780 | 9.026684 | 3.51 | 2 | Both raw solutions lack enough positive independent events. |
| 183642107 | 3.143099 | 22.22 | 5 | A double-period secondary alone cannot reject the shorter period, but the short-period single-shape fit is very poor (chi2 1682 for 279 points) and longer alias folds deeper narrow events. Common eclipse/shape model needed. |
| 269994642 | 3.492953 | 15.13 | 6 | Primary followup already active; double-period alias, not an additional planet. |
| 204497901 | 3.519317 | 8.00 | 5 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 181050705 | 2.308373 | 7.69 | 4 | Excessive duty cycle and large orbital variability; no compact robust transit lead. |
| 269772299 | 4.290918 | 5.40 | 3 | Alternate period is not a demonstrated additional planet. Compare alias geometry, parity depths and true-short-period secondary tests; a doubled-period secondary alone does not falsify a shorter-period planet. |
| 181043205 | 6.003692 | 2.86 | 2 | Weak raw signal and duty/variability concerns; low priority. |
| 29634861 | 2.861289 | 0.00 | 0 | Extra peak has zero local events; original raw SNR only 2.98. |
| 197936765 | 5.906035 | 2.62 | 3 | Weak raw events and variability; no robust periodic transit. |
