# Independent follow-up of public sector-107 TCEs

These are already-public threshold-crossing signals, independently examined here for possible unrecognized planetary interpretations. **No credible new planet is established.** Read the per-target waveform plots and measurements; a large formal diagonal score is insufficient.

Native QUALITY=0 finite positive SAP/PDC flux is preserved and time references are verified. Separate event baselines use linear and quadratic models, central box depths and symmetric nearby out-of-transit coverage. Uncertainties have a native adjacent-sample scatter floor. Alternative one/two-day gap-aware medians, a targeted fine-period BLS, odd/even and half-phase measurements, Lomb–Scargle and out-of-dip one/two-harmonic models expose variability and instability. Raw data and processing provenance are in `data/tce107_followup_observations.json`.

| TIC | Sector | Flux | Local linear depth ppm | Formal S/N | Event-scatter S/N | Events | Local quadratic depth ppm | Strong waveform period d |
|---|---|---|---|---|---|---|---|---|
|159923666|105|PDCSAP_FLUX|-1652.5|-32.6|-0.4|2|1585.6|2.686345|
|159923666|105|SAP_FLUX|-1688.5|-29.6|-0.4|2|1568.1|unavailable|
|159923666|106|PDCSAP_FLUX|-7903.0|-204.9|-5.3|3|-1717.7|2.793732|
|159923666|106|SAP_FLUX|-7960.8|-209.7|-5.3|3|-1768.8|unavailable|
|159923666|107|PDCSAP_FLUX|2128.8|42.1|1.2|3|1434.3|2.734219|
|159923666|107|SAP_FLUX|2158.7|46.2|1.2|3|1428.3|unavailable|
|178155030|105|PDCSAP_FLUX|22783.4|203.0|80.8|44|14632.2|0.486693|
|178155030|105|SAP_FLUX|22309.9|200.0|78.6|44|14368.7|unavailable|
|178155030|106|PDCSAP_FLUX|21214.6|177.3|35.4|37|14342.7|0.486529|
|178155030|106|SAP_FLUX|20790.8|177.1|35.6|37|14059.4|unavailable|
|178155030|107|PDCSAP_FLUX|21691.0|158.3|72.8|40|14349.8|0.486638|
|178155030|107|SAP_FLUX|20590.4|158.4|73.0|40|13630.0|unavailable|
|100777592|33|PDCSAP_FLUX|15098.8|38.9|2.6|4|5913.2|0.430795|
|100777592|33|SAP_FLUX|6056.2|39.3|2.5|4|2356.8|unavailable|
|100777592|98|PDCSAP_FLUX|5193.2|17.0|0.5|7|3643.0|0.430323|
|100777592|98|SAP_FLUX|2057.6|16.8|0.5|7|1462.0|unavailable|
|100777592|107|PDCSAP_FLUX|14866.7|26.9|1.2|2|11652.2|0.430674|
|100777592|107|SAP_FLUX|5836.6|26.5|1.2|2|4599.9|unavailable|

TIC 159923666: a large approximately 2.734-day waveform explains most variance; the public approximately 8.07-day interpretation groups different waveform phases. Individual sector-107 depths vary strongly, and the predicted earlier-sector primary events are often brightenings. Stable planetary transits are unsupported.

TIC 100777592: a strong approximately 0.431-day waveform and highly inconsistent selected-event depths undermine the approximately 5.85-day periodic interpretation. CROWDSAP is approximately 0.397; PDC depths are scaled relative to raw SAP as expected for strong crowding correction. This does not identify the actual variable source. Pixel localization would be needed to attribute the variability; no planet is promoted.

TIC 178155030: the approximately 0.487-day modulation repeats in sectors 105–107, with a large asymmetric full-phase waveform and an additional narrow minimum. The dip remains real after one/two-harmonic removal, at approximately 12,163–12,679 ppm for two harmonics; stellar modulation alone is therefore insufficient to dismiss the localized feature. **The decisive prior check identifies it as a previously published complex periodic variable:** the Bouma et al. author table matches both TIC 178155030 and Gaia DR3 5050451154310013056, lists period 11.67 hours, Tucana–Horologium association, age 45 Myr, radius 0.49 solar, and four earlier observed sectors. The paper is [Transient Corotating Clumps Around Adolescent Low-Mass Stars From Four Years of TESS](https://arxiv.org/abs/2309.06471), initially posted 2023-09-12, published in The Astronomical Journal 167(1), 38 ([DOI 10.3847/1538-3881/ad0c4c](https://doi.org/10.3847/1538-3881/ad0c4c)) on 2023-12-28. This recognized waveform and its narrow dips must not be promoted as a newly recognized planetary signal. The public author table, manuscript, source metadata and byte hashes are preserved locally; the table URL is pinned to a Git commit in its manifest. This does not exclude an unrelated additional planet around the star.

The coordinate/catalog check remains narrower than a literature novelty investigation. The source TCE table was already public, so none of these signals can be called an original new detection here. Physical planet fitting, injection completeness and centroid/pixel diagnostics cannot repair a waveform interpretation that fails initial screening; any later revival needs independent reduction and targeted physical/literature tests.

Reproduce using `python src/triage_tce107.py` then `python src/vet_tce107.py`. The script does not download, overwrite or remove source FITS. All measurements, event-level values, hashes linked through the follow-up manifest and diagnostic plots are retained. These three outcomes do not exclude planets around their hosts or validate a population-level false-positive rate.
