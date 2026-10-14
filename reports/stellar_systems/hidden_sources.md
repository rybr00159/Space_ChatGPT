# A short-period eclipsing signal hidden in the image of giant TIC 260352037

Checked 2026-10-14. This is a targeted stellar-system investigation of three
previous broad-search leads, not a new survey of all 159,505 stars.

**TIC 260352037 has a convincing repeated eclipse signal in its blended TESS
image. Its stellar components are not spatially identified or newly discovered
by this analysis.** Two other apparent giant-star box signals fail a fixed-phase
holdout test and are not promoted.

## Independent evidence

The raw QLP SAP observations for TIC 260352037 span 25 sectors: 1–13, 27–29,
and 31–39. All original observations were retained and their recorded byte
hashes and FITS integrity were checked before processing.

A narrow eclipse template, period and epoch were selected using only the first
two sectors. That fixed ephemeris predicts dips in every one of the 23 later
sectors, with individual sector raw local-baseline signal/error ratios between
about 7 and 14. This is stronger evidence than independently optimizing an
epoch in every sector. The noisy bright giant has many oscillation peaks; the
coherent, short dip need not be the largest Lomb–Scargle peak.

Seventy-six independently measurable eclipse times give a photometric
recurrence period **0.880788466 ± 0.000002074 days** and primary epoch
**BTJD 1326.593218 ± 0.001879**, where BTJD = BJD(TDB) − 2457000. These are
conditional fixed-template uncertainties, inflated by timing scatter. The
initial all-sector cadence-integrated trapezoid has central depth
**1,710 ± 64 ppm**, duration **1.85 ± 0.04 hours**, and ingress fraction 0.434.
It is a phenomenological shape model, not a measurement of companion radius.
The shape fit used the early-sector trial period; the independently timed
linear ephemeris is the recommended timing solution.

Raw independent polynomial baselines at primary and pre-specified half phase
produce these pooled window-average depths:

| Local baseline | Primary depth (ppm) | Half-phase depth (ppm) |
|---|---:|---:|
| Linear | 980 ± 24 | 483 ± 30 |
| Quadratic | 972 ± 26 | 445 ± 38 |
| Cubic | 977 ± 26 | 437 ± 39 |

The primary and secondary measurements contain 517 and 522 covered events,
respectively, in the quadratic analysis. The error calculation includes
per-event scatter and sector heterogeneity; the resulting secondary ratio is
11.65. This is a targeted conditional ratio, **not** a global false-alarm
probability or an independent discovery significance. The primary central
shape depth exceeds its window mean because ingress and egress dilute that
mean.

Twenty-four other orbital phases were inspected; none had a pooled positive
dip larger than 171 ppm. They share data and some control baselines contain
the actual eclipses, producing negative depths down to −574 ppm. They must
not be treated as an independent null sample or used to claim a p value.

The full-phase stacked curve shows two isolated dips, with a relatively flat
inter-eclipse baseline after a two-day median normalization. The raw local
measurements above do not depend on that median-normalized display.

## Where the signal is

Three public 11 × 11 native TessCut cubes, sectors 1, 28 and 39, add spatial
information. The target falls on verified science silicon in all three; no
prescan or collateral rows were used. Independent per-event quadratic,
background and pointing-nuisance fits produce deficit images. Cadence
integration is included in the eclipse template.

| Sector | Primary events | Quadratic source offset from giant catalog position (pixels) | Half-phase events | Half-phase offset (pixels) |
|---|---:|---:|---:|---:|
| 1 | 22 | 0.145 | — | — |
| 28 | 28 | 0.042 | 26 | 0.109 |
| 39 | 29 | 0.056 | 31 | 0.108 |

A TESS pixel is about 21 arcseconds. Thus both eclipse types occur in or near
the unresolved image of the bright giant, rather than on the bright resolved
field stars tens of arcseconds away. The images do **not** establish that the
giant itself eclipses. Whole-event plus synchronized five-cadence residual
bootstraps describe statistical centroid dispersion, but exclude PRF, WCS,
proper-motion and model systematic uncertainties. Pointing nuisance columns
shift the sector-28 primary fit by 0.29 pixel, demonstrating that systematic
uncertainty cannot be ignored. Small-region Gaussian fits to the reference
image also move with the region selected; subpixel source attribution is not
secure.

The currently catalogued giant has radius 10.1065 solar radii, temperature
4,809 K, Gaia G = 11.8018 and parallax 0.59728 ± 0.00959 mas. A 0.8808-day
Keplerian orbit has semimajor axis only about 3.9–5.6 solar radii for an
illustrative total mass of 1–3 solar masses. Even a doubled period remains
inside this giant radius at the upper mass of that illustrative range.
The mass is not measured here. These figures make an ordinary detached
binary on the catalogued giant physically implausible, and motivate an
unresolved faint stellar source. They do not prove a particular source or
rule out errors in stellar inference.

There is an interesting resolved near neighbor:

- Gaia DR3 **5496090260793489024**, G = 17.1632, lies 4.21 arcseconds from
  giant Gaia DR3 **5496090260794096128**. Its exact TIC crosswalk is
  **TIC 737466840**.
- Its proper motion (−5.001, +17.608) mas/year agrees with the giant's
  (−4.864, +17.614) within 1.42 standard errors in right ascension and
  0.06 in declination. Parallaxes differ by 1.17 combined standard errors.
- If associated at the nominal giant distance, projected separation is
  approximately 7,000 AU. This is compatible astrometry for a wide pair,
  not proof of gravitational binding.
- Its Gaia G flux is approximately 0.0072 times the giant's. A deep eclipse
  on this faint star could yield the observed diluted signal. Gaia G is
  not the TESS band and aperture coupling is unmeasured; this is only a
  plausibility calculation.

**The deficit centroids do not consistently favor that 4.21-arcsecond neighbor.**
A closer unresolved binary, the neighbor, or a more complex blend remains
possible. Gaia DR3 5496090256498964096 = **TIC 260352041** at 7.73 arcseconds
is another faint source, but its proper motion is markedly different from
the giant.

## Recognition and the other two targets

An independent recognition audit found no exact TIC match for the giant or
the two nearby faint sources in the Shi (2022), Tang (2025), Battley (2022)
or actual EB tables 3/4 of Kostov (2025). Their exact Gaia identifiers have
no checked NSS, eclipsing-binary or variability-classifier solutions. The
120-arcsecond cones are empty in VSX, live ASAS-SN, Rowan (2022) detached
binaries, WISE (2021), SB9, WDS and Tian (2020) wide pairs. Current ExoFOP
TOI/CTOI, companions and uploaded files are empty for all three TICs;
Villanova reports each is not in its catalog. The giant's SIMBAD alias is
TYC 8542-438-1.

The complete source manifests, crosswalk and finite scope are recorded in
`reports/stellar_systems/prior_audit_supplement.md`. These are meaningful
non-recognition results within the checked resources, **not novelty proof**.
Finite arXiv alias queries do not cover all full-text papers or supplementary
tables. Any positive recognition takes precedence over this exploratory
assessment. No confirmed new star, new bound system or first published EB
claim is made here.

For **TIC 271697430**, the first-sector approximately 3.059-day box peak does
not maintain a fixed phase across its 25 cached sectors. Independent later
optimizations mostly reproduce weak troughs at unrelated phases. This is
consistent with its previously published oscillating-giant context, rather
than evidence for an additional star.

For **TIC 406280158**, three cached sectors give inconsistent fixed-phase
support at approximately 10.12 days. One later sector is negative at the
predicted eclipse, while the other has incomplete sampling. Its Gaia RV
variability is interesting but does not identify an eclipse orbit or companion.
The early-only ephemeris uncertainty is substantial; these observations do not
exclude every binary interpretation, but offer no credible new system claim.

## Reproduction, preservation and follow-up

Code: `scripts/stellar_hidden_sources.py`. Figures and quantitative outputs:
`candidates/stellar_systems/hidden_sources/`. New astrometry and exact request:
`data/stellar_systems/hidden_sources/tic260352037_astrometry*`.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/workspace/openrouter-check-venv/bin/python scripts/stellar_restore_lightcurves.py --download
/workspace/openrouter-check-venv/bin/python scripts/stellar_hidden_sources.py --mode coherence
/workspace/openrouter-check-venv/bin/python scripts/stellar_hidden_sources.py --mode photometry
/workspace/openrouter-check-venv/bin/python scripts/stellar_hidden_sources.py --mode controls
/workspace/openrouter-check-venv/bin/python scripts/stellar_hidden_sources.py --mode pixels --sectors 1 28 39 --bootstrap 200
/workspace/openrouter-check-venv/bin/python scripts/stellar_hidden_sources.py --mode pixels --sectors 28 39 --bootstrap 150 --event-phase .5
```

The coherence audit reads the retained light curves listed with exact URLs and
SHA-256 hashes in the old giant-survivor QA dossiers. If the cache is lost,
the restoration command above retrieves only missing files and verifies their original hashes before running. TessCut
retrieval is built into the pixel mode; raw hashes, full headers, retrieval
URLs, checksums and FFI identifiers are preserved in each sector directory.
Raw cubes are retained outside Git under
`/workspace/tess-research-data/stellar_systems/hidden_sources/tic260352037/`.
Their total size is about 21.6 MB. The environment's independent storage
persistence is unverified; all retrieval information and compact diagnostics
are in Git.

The most useful next observation is high-resolution differential photometry
covering both predicted eclipse types, resolving the giant from the star at
4.21 arcseconds. Alternatively, calibrated TESS PRF scene fitting across
several camera orientations may help. Spectroscopy of the faint source,
after spatial identification, can test an inner spectroscopic binary. A
resolved known-wide-pair association and a new inner eclipsing binary are
separate claims requiring separate evidence.
