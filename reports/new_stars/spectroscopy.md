# Repeat-velocity search for unseen companions

This examines existing visible stars, seeking overlooked multiplicity rather
than new positions. GALAH DR4 is public spectroscopy independent of the earlier
TESS screening. The retained quality predicate selects dwarf/subgiant-like
temperatures 4000–7000 K, logg>3.7, high S/N, usable errors and unflagged visits.
3430 stars with at least three selected visits supply 13176 RV measurements.
An aggregate 356869 quality-visit count includes single-visit stars and is not
the analyzed-system count.

Forty amplitude-ranked cases exceed 20 km/s over at least two days.39 are already
binary-flagged in GALAH's all-star table. The other has reduction problems and
only one completely clean visit. A stricter all-star/visit-unflagged search
contains 1799 repeated stars, with seven ranges above 3 km/s. Three are Gaia NSS
solutions, two have strong published Gaia RV variability; two remain provisional:

| Source | Existing identifiers | Catalogue range | Main limitation |
|---|---|---:|---|
|4690783166942834816|UCAC 4 093-001240;2MASS J 00592586-7129339;TIC 51939231|7.75 km/s|Sparse 2014/2023 epochs; known rotational variable; wavelength-dependent line shifts |
|5437472031746145536|TYC 7162-1658-1;2MASS J 09272542-3507283;TIC 24731098|5.80 km/s|Three clean nights; broad F-star lines; stronger masking/region dependence |

The first is already a VSX ROT source with a 20.766-day period. SB9, actual
Traven 2020's 12760 rows, Gaia aliases/NSS, WDS, VSX, ASAS-SN and SIMBAD bibliography
checks have successful positive retrieval controls. No exact recognized binary
was found for these two in that finite scope. A known Traven SB2 is 111 arcsec
from the second target and is a different Gaia source. An initial erroneous
ASAS-SN table name produced a service error, not a valid absence; corrected
catv2021 cones and a true positive control replace that interpretation.

2388 same-field spectral rows on 971 sources provide 94 and 54 eligible old/new
paired controls. Median offsets are+0.141/−0.105 km/s, robust scatter 0.502/0.449 km/s.
A uniform epoch zero point is insufficient, but individual spectral failures
and true binaries remain in the tails; these controls do not calibrate global
false positives.

Four original published observed-flux products plus their covariance/result
files were retained. **Their wavelengths are already in the fitted rest frame.**
The published transformation was inverted with matching fit velocities before
comparing data in separate blue/green/red regions. Values are conditional on
that reconstruction, not independent instrument calibration. Regional shifts
span−7.48..−8.83 km/s for the first object; the second spans−7.20..−10.12 km/s with
masking changes. Formal errors of order 0.03–0.13 km/s substantially understate
these discrepancies. Synthetic 8 km/s shifts only validate sign and units.
No precise new binary RV solution is justified.

Three original QLP products (TIC 24731098 sectors 35/62;TIC 51939231 sector 68)
were checked for rapid photometric variability using raw/local and corrected
flux. The hotter object's dominant frequencies do not reproduce consistently;
no clear stable high-amplitude pulsation is established. One-day raw filtering
attenuates slow rotation, so it cannot rule out the known 20.8-day modulation.
These same-aperture data do not localize variability or exclude all activity.
Exact target-named SPOC queries returned no rows; that is a finite query result,
not proof no other reduction exists. A QLP sector 27 product 404 is preserved.

An unseen companion remains possible, along with calibration/profile/activity
alternatives. Sparse velocities provide no orbital period, mass function or
component mass. Therefore **neither is evidence for a newly discovered star**.
Read the complete data in `candidates/new_stars/rv/` and recognition snapshots
in `data/new_stars/rv_recognition/`.

```sh
python src/spectroscopic_companion_scout.py
python scripts/rv_lead_recognition.py
python src/rv_spectral_validation.py
python src/rv_photometric_veto.py
python scripts/summarize_new_stars.py
```

Retrieval flags restore missing static products. Query response snapshots and
their failed variants remain retained; no paid models were used.
