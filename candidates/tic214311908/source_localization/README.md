# TIC 214311908: source assignment remains unresolved

The native S66 and S93 difference images support light loss near the target
and its close neighbors at the independently selected 5.6915-day ephemeris.
They **do not identify TIC 214311908 as the host**. A nearby dwarf and faint
background sources remain plausible alternatives. This audit uses only the
S66/S93 maps; no held-out-sector observations were inspected.

## Proper motion and conditional source fits

Gaia DR3 positions were propagated from J2016.0 to each sector's median native
observing time, using the measured angular proper motion. The raw cutout WCS
then converts each propagated position to zero-based pixels. The original
`RA_OBJ` target cross is displaced from the propagated position by **0.161
pixel in S66 and 0.175 pixel in S93**, about 3.4–3.7″. Its coordinates are close
to the target's J2000 position. This offset matters when considering neighbors
only 3–4″ away.

| Source | S66 position (x,y) | S93 position (x,y) |
|---|---|---|
| Target TIC 214311908 | (4.0685, 4.9079) | (4.0809, 4.6714) |
| Near-comoving TIC 214311910 | (4.6324, 4.9073) | (4.6452, 4.6734) |
| TIC 1515709327, 3.15″ neighbor | (4.1646, 5.0240) | (4.1885, 4.7803) |
| TIC 1515709326, 4.27″ neighbor | (4.2980, 4.9706) | (4.3219, 4.7270) |
| Mira TIC 214311906 | (4.9634, 4.9556) | (4.9872, 4.7148) |

All 16 Gaia sources within 12″ and the named Mira/two neighboring EBs were
compared. Each fixed source receives the same positive pixel-integrated
elliptical Gaussian plus background plane, with its widths and amplitude
fitted. Two prescribed-width alternatives are also retained. There are **304
fixed-source comparisons**: 19 possible sources × two sectors × eight
overlapping reduction variants. These are conditional descriptive fits, not
304 independent tests or calibrated source probabilities.

The table gives fixed-source χ² minus target χ² across the reduction variants;
negative values favor the alternative within this model.

| Alternative source | S66 Δχ² range | S93 Δχ² range |
|---|---:|---:|
| Near-comoving TIC 214311910 | −5.37 to +2.39 | −0.33 to +1.17 |
| 3.15″ TIC 1515709327 | −2.22 to −0.57 | −0.80 to −0.05 |
| 4.27″ TIC 1515709326 | −3.77 to −0.03 | −1.99 to −1.56 |
| Mira TIC 214311906 | −3.51 to +5.49 | +10.68 to +12.95 |
| Known 0.316083-day EB at 26.63″ | +7.71 to +17.85 | +33.91 to +46.72 |
| Known 0.260926-day EB at 62.12″ | +28.70 to +52.42 | +50.40 to +75.30 |

The target, its near-comoving neighbor and the close faint sources fit almost
equally well. S66's preferred center and shape change with background and local
baseline treatment. In S93 the free center is 0.268–0.359 pixel from the
propagated target; one close faint source is even nearer, while the target and
near-comoving source remain statistically indistinguishable under the adopted
model. The known EBs are a worse spatial explanation of these particular dips;
that is consistent with the separate much shorter-period variability. The Mira
is also a worse spatial fit in S93, but Gaussian/diagonal-noise approximations
are insufficient for an absolute exclusion or a quoted sigma.

A pixel peak S/N around 5–6 is not a 5–6σ host identification. These maps reuse
the discovery photons and selected ephemeris. A calibrated TESS PRF, an actual
WCS residual measurement and better treatment of correlated pixel/event noise
would be needed before converting position fits into source probabilities.

## QLP normalization changes the depth interpretation

The official MAST QLP page links Huang et al. 2020, arXiv:2011.06459v1.
Its reduction step 7 states that QLP adds expected source flux based on the
TIC T magnitude and instrumental zeropoint to the measured difference flux.
The paper calls this **"equivalent to deblending the flux time series assuming
the variations observed are from the target star."** It also warns that
incorrect TIC magnitudes change amplitudes, that this uncertainty is absent
from the reported point errors, and that incomplete PSF containment in the
smallest aperture can underestimate variation amplitudes.

Therefore **the TIC contamination ratio 2.23 must not simply be applied as a
3.23× correction to a normalized QLP depth**. The QLP product already uses a
target-flux/deblending assumption. A different host requires its own flux,
aperture coupling, reference normalization and stellar properties. The
S66/S93 raw QLP headers have TESSMAG 13.0135, TICVER 8.2, and no measured
`CROWDSAP` or `FLFRCSAP` value. The cited paper describes the method; it is not
a new calibration of the current apertures.

This also makes a concrete false-positive scenario possible: the 3.15″ neighbor
TIC 1515709327 has T=16.8469, or approximately **2.93% of the target's TIC T-band
flux**. With similar aperture coupling, a roughly 50% eclipse of that star
could produce approximately a 1.5% target-normalized QLP dip. Its native-map
position fits at least as well as the target. The 4.27″ neighbor's nominal flux
is only about 1.19% of the target's; explaining a larger dip would challenge
that flux estimate, but a hard exclusion needs an actual response/calibration
model. Neither a centroid match nor a catalogue nonentry resolves these cases.

## Resolved public imaging and useful follow-up

A public Legacy Survey viewer **DECaPS2 grz cutout** was retrieved and retained.
The z-band image clearly resolves the target, the bright Mira, and the
near-comoving neighbor; the close faint sources appear in the surrounding
field. It confirms substantial real crowding. The coadd has no unique epoch,
so its Gaia labels deliberately use J2016 positions without a false exact
epoch claim. This static image identifies neither an eclipse nor a new star.

![Resolved public field](decaps2_resolved_context.png)

The most decisive practical next observation is resolved photometry during
the predicted event, monitoring **the target, TIC 214311910 and especially
TIC 1515709327 simultaneously**. About 1″ seeing and appropriate apertures/PSF
fitting could separate the 3.15″ neighbor; a deep eclipse on that faint source
would be much easier to distinguish than another blended TESS transit.
Spectroscopy of the bright target alone could miss a background EB.

## Reproduction

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  /workspace/openrouter-check-venv/bin/python scripts/tic214311908_source_localization.py
```

This uses the preserved S66/S93 raw cutouts and maps, checks their hashes, and
checks the frozen Gaussian helper's hash. Results are in
[`fixed_source_comparisons.csv`](fixed_source_comparisons.csv), the two sector
JSON files, and [`summary.json`](summary.json). Source overlays are
[`sector66_gaia_overlays.png`](sector66_gaia_overlays.png) and
[`sector93_gaia_overlays.png`](sector93_gaia_overlays.png). The archive image,
PDF, extracted methodology text and HTTP URL/checksum manifests are under
`data/tic214311908/source_localization/`; their manifests permit retrieval if
the environment is replaced. A failed optional DECaPS homepage request is
preserved; it did not prevent retrieval of the actual public image.
