# Independent stellar context and attempts to falsify the two system claims

Public snapshots checked 2026-10-14 UTC. This analysis extends the previous
recognition audit with Gaia DR3 astrophysical parameters, actual RV-statistic
definitions, quality diagnostics, infrared photometry, modern spectroscopy
catalogs and a controlled wide-binary query. It does not replace pixel source
attribution or a physical eclipse fit.

**The context supports a dwarf-like background eclipsing source near TIC
273879360 and independently establishes that TIC 260352037 is a red giant.
The giant's 0.8808-day eclipses therefore require a different interpretation
from an ordinary detached binary on that giant. A faint close source is
physically plausible, but neither its identity nor binding to the giant is
established. A deeper check confirms that TIC 273879360 was previously
publicly flagged for eclipse-like features by machine learning, although no
vetted binary classification or period was found.**

## TIC 273879360: a catalogued dwarf-like source, without component measurements

The Gaia parallax implies 812.35 ± 11.05 pc by direct inversion before
zero-point correction. Its high parallax signal/error makes this a useful
distance scale, not a full Bayesian distance determination. Absolute G before
extinction is 5.155 mag. The previous TIC catalog estimated T_eff = 5,420 K,
radius 0.9595 solar radii and mass 0.94 solar masses. Gaia's single-star
GSP-Phot model instead gives T_eff = 5,832.9 K, log g = 4.4695, [M/H] = −0.5114,
A_G = 0.3765 mag, distance 747.8 pc and radius 0.8577 solar radii. FLAME gives
radius 1.0313 solar radii and luminosity 0.6758 solar luminosities; its flag
`20` means mass/age are unavailable and the parallax was used for the FLAME
distance scale. It does not itself flag the input atmospheric parameters.

These are different catalog models of the combined source. They consistently
suggest a broadly solar-size dwarf-like object, but their differences exceed
some internal model intervals. Once two luminous components are plausible,
single-star radius, metallicity and temperature estimates cannot be treated
as measured properties of either component. No mass ratio, age or precise
component spectral type follows from these entries.

The exact 2MASS counterpart J02090051−7119235 has good `AAA` photometry:
J = 13.432 ± 0.024, H = 12.976 ± 0.024, K_s = 12.902 ± 0.023 mag. AllWISE
J020900.57−711923.4 has W1 = 12.908 ± 0.024 and W2 = 12.926 ± 0.023 mag,
with no artifact flags in those bands. W1−W2 = −0.018 ± 0.033 mag supplies
no convincing warm infrared excess. W3 has 0.369-mag error and W4 is an upper
limit, so neither establishes a disk or unusual infrared source. These
unresolved, non-simultaneous measurements do not decompose the binary SED.

Gaia photometric contamination counters are clean for this source: no BP/RP
contaminated or blended transits, RUWE 0.957, IPD multiple-peak fraction 0.
The Riello et al. (2021) color-corrected BP/RP flux excess is only +0.00358.
The normal astrometric fit does not exclude a close binary whose orbital
photocenter is small or inadequately sampled.

**A tempting RV argument fails:** Gaia reports RV = 48.06 ± 9.68 km/s from
20 transits, but `rv_method_used = 2` means a faint-star combined
cross-correlation measurement. The error is not epoch RV scatter. Gaia
provides no RV-constancy statistic or orbital solution here. It is therefore
incorrect to call the large RV error an independent spectroscopic-binary
detection or to use the catalog RV as a well-measured systemic velocity.

## TIC 260352037: independent spectroscopic giant evidence

The bright source has inverse-parallax scale 1,674.25 ± 26.89 pc and absolute
G = 0.683 before extinction. Gaia GSP-Spec, based on RVS spectroscopy rather
than the TESS eclipse, reports:

| Parameter | Published value and internal 16th–84th percentile bounds |
|---|---|
| T_eff | 4,834 K, 4,812–4,895 K |
| log g, cgs | 2.51, 2.40–2.65 |
| [M/H] | −0.12, −0.16 to −0.06 |
| FLAME spectroscopic-input radius | 9.447 solar radii, 9.176–9.697 |
| FLAME spectroscopic-input luminosity | 44.523 solar luminosities, 42.694–46.073 |

The alternative Gaia ANN spectroscopic solution gives 4,626 K and log g
2.75. Its disagreement cautions against treating the quoted internal bounds
as complete errors, but both solutions support a giant. The FLAME_spec flag
`2` means no available mass/age; its evolution-stage value is not a measured
component age. This spectroscopic evidence supports the old TIC radius
10.1065 solar radii rather than rescuing a solar-size eclipsing host.

For illustration, at P = 0.880788 d Keplerian semimajor axes for total masses
1–3 solar masses are approximately 3.87–5.58 solar radii, smaller than the
giant. At twice that period they remain approximately 6.14–8.86 solar radii.
These calculations use an assumed mass range, not measured masses. An
ordinary detached orbit around the catalogued giant is physically implausible;
a faint blended system, misidentified source, or unusual configuration must
be considered. Roche-lobe constraints would be stricter than the simple
separation-versus-radius comparison.

The giant's Gaia RV is 54.869 ± 0.362 km/s from 21 transits. Here
`rv_method_used = 1`, RV-constancy p = 0.8066 and renormalized goodness of fit
= −1.134. It does not meet Gaia's suggested variable-RV selection
(p < 0.01, goodness of fit > 4, at least 10 transits). This is no reported
significant epoch-RV variation on the measured giant, not proof of singleness
or an upper limit on a faint neighbor's binary orbit. A 5.395-km/s robust
peak-to-peak range alone is not an orbital amplitude or a detection.

Its good 2MASS J/H/K_s = 10.221/9.707/9.558 and AllWISE W1/W2 = 9.502/9.581
belong to the bright giant's unresolved image. The faint sources at 4–8
arcseconds lack secure separate infrared photometry in these products.
Nearby low-quality 2MASS detections with upper limits/confusion flags must
not be substituted for the close companion's SED.

## The possible faint source: informative astrometry, uncertain binding and colors

TIC 737466840, Gaia DR3 5496090260793489024, is 4.2151 arcseconds from the
giant. At the nominal giant distance that is 7,057 AU projected separation.
Its parallax differs by 1.174 combined standard errors. Accounting for each
star's published RA/Dec proper-motion covariance gives a zero-relative-motion
chi-square of 2.333 for two dimensions (p = 0.311). This is compatible
astrometry under independent source errors; shared Gaia systematics and
orbital motion are not modeled.

Compatibility must not be upgraded to proof of gravitational binding. The
nominal transverse relative velocity is 1.086 km/s at the giant's distance.
Escape speeds evaluated at the projected separation are only 0.501, 0.709
and 0.868 km/s for assumed total masses 1, 2 and 3 solar masses. Actual
three-dimensional separation could be larger. The proper-motion errors are
large enough that this central-value comparison cannot show the pair is
unbound either; a radial velocity of the faint star and improved relative
astrometry would matter. No binding probability is estimated.

The faint source's colors warrant strong caution. Gaia BP/RP flux excess is
1.60385, with color-corrected excess **+0.37539**, versus approximately zero
for the giant and source 273879360. Its BP/RP products include contaminated
transits and one blended BP transit. Single-star GSP-Phot gives 4,793.7 K,
distance 996.7 pc and radius 0.5282 solar radii, whereas direct parallax
inversion gives approximately 1,901 pc and the alternative PHOENIX model
radius is 0.8524 solar radii. FLAME radius is 0.9332. Close bright-star
contamination and model assumptions can affect these estimates; they cannot
define component masses or establish coevality with the giant.

The TIC T-band magnitudes allow a conditional light-budget test. With equal
aperture throughput, the giant plus TIC 737466840 mixture gives the faint
source 0.006685 of total flux. A 1,710-ppm central eclipse needs a plausible
25.6% intrinsic loss of the faint source's light. The measured primary and
secondary window means would need 14.5% and 6.66%. This shows feasibility,
not localization. Catalog colors, differential aperture throughput and other
blend light are not measured here.

TIC 260352041 at approximately 7.7 arcseconds is fainter still, T = 18.1961.
Equal throughput would require 110.8% central loss, impossible under that
restricted two-source model. Only approximately 11% higher relative aperture
coupling could remove this contradiction, so it is a disfavored example,
not a definitive exclusion. Its different proper motion also argues against
a simple comoving association with the giant.

## Prior recognition: controls and a misleading assumed-binary product

**The strongest novelty challenge is a prior automatic eclipse flag.** TIC
273879360's exact row is present in Kostov et al. (2025) Table 2, the roughly
900,000 unvetted/unvalidated targets. Reading only the table's four columns
(TIC, position and T magnitude) can mislead: this is not a generic inventory
of observed stars. The [paper abstract](https://arxiv.org/abs/2506.05631)
explicitly says the neural network found eclipse-like features with score
greater than 0.9 for these targets. Section III.4 calls it a catalog of
light curves with eclipse-like features and explicitly distinguishes it from
an EB catalog; many entries can be contaminated, instrumental or planetary.

Our current evidence is a substantially stronger independent characterization
of a **previously machine-flagged source**, including an identified background
source, recurrence period, primary and secondary eclipses and independent
ground photometry. It should not be presented as the first recognition of
any eclipse-like variability or a completely unflagged system. No matching
entry was found in the vetted Tables 3/4, and this finite absence does not
establish first publication of a validated binary system.

Gaia MSC contains temperatures for two hypothetical components of source
273879360 and TIC 737466840. It also contains a poor fit to the giant, with
distance 244.7 pc inconsistent with its parallax. The [official MSC scope
documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Data_analysis/chap_cu8par/sec_cu8par_apsis/ssec_cu8par_apsis_msc.html)
states that it processed all G ≤ 18.25 sources with parallax and BP/RP spectra:

> “We do not provide any information whether a respective source actually
> is an unresolved binary; we just provide the best parameter fits assuming
> it is.”

**These entries are not prior binary recognitions or independent companion
measurements.** Gaia DSC's combined binary-class probabilities are only
0.000104 for source 273879360, 0.0000144 for the giant and 0.00000442 for TIC
737466840. Those low-resolution classifications cannot override resolved
eclipse evidence or prove any source is single.

Ten-arcsecond target cones returned no source spectra/catalog rows in GALAH
DR3, GALAH DR4, APOGEE DR17 or RAVE DR6, including the giant's close-source
neighborhood. The identical cone query endpoint and parameters recovered
positive control sources in every catalog (2, 2, 4 and 3 returned rows,
respectively, including catalog subtables). These nulls mean no data in
those specific checked releases, not no previous spectroscopy anywhere.

Exact Gaia IDs returned no pair in El-Badry, Rix & Heintz (2021)'s large
Gaia EDR3 wide-binary catalog. An exact-ID positive pair was successfully
recovered. Crucially, that catalog selects parallax > 1 mas: **the giant and
its close faint sources are outside its distance selection**, so their null
cannot support novelty. Source 273879360 is within the nominal distance
selection but its absence concerns resolved wide companions, not close
eclipsing multiplicity. The previous Tian/WDS/VSX and finite literature
audit remains the relevant bounded recognition check; no first-discovery
claim follows from the new nulls.

## Independent ground photometry and reproducibility

The companion [ground-photometry report](../../candidates/stellar_systems/deep_context/ground_photometry/README.md)
retrieved 552 archived ASAS-SN V-band observations of AP36443078, 0.158
arcseconds from TIC 273879360, from May 2014 to September 2018. Exported CSV
and independent embedded HTML photometry agree exactly. At the fixed TESS
ephemeris, all 11 primary-window observations in nine distinct events dim;
eight of those events predate TESS. A descriptive fixed-shape primary
amplitude is 35.2%, with night-block bootstrap interval 28.4–39.7% and
baseline alternatives approximately 32.9–35.2%. The secondary is weaker and
camera-dependent. This independent earlier optical recurrence substantially
weakens instrumental explanations, while different bandpass, dilution and
unresolved blends still prevent component measurements.

The same work retrieved resolved SkyMapper observations, but none samples
either predicted eclipse, so it is not a confirmation or nondetection.
The current SkyPatrol light-curve service announced maintenance beginning
13 October 2026; the older public photometry release remained accessible.
An accessible target inventory is not itself a variable-star classification.

Run the new retained-data derivation:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python scripts/stellar_deep_context.py
# Add --download only to retrieve missing public snapshots.
```

The [machine-readable context](../../candidates/stellar_systems/deep_context/context_summary.json)
and [source table](../../candidates/stellar_systems/deep_context/source_properties.csv)
retain formulas, numerical results, catalog-model values and SHA-checked
provenance. Raw response bytes and public retrieval manifests are in
`data/stellar_systems/deep_context/`. A failed Gaia metadata probe and an
APOGEE positive-control cone using unsupported scientific-notation RA are
retained; corrected retries succeeded and the failures are never interpreted
as zero catalog membership. Positive-control validation, failed-response
rejection and source hash verification run in the default derivation.

Useful external sources: Gaia DR3 [data-model definitions](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html),
Riello et al. 2021 [photometric calibration](https://doi.org/10.1051/0004-6361/202039587),
[El-Badry2021 selection](https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/506/2269/ReadMe).
No paid API calls, external submissions or observation requests were made.
