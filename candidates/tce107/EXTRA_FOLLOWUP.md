# Further independent follow-up of sector-107 public TCEs

Two additional stars were selected from the frozen public-TCE metadata ranking before observing their flux: TIC 49617264 (MES 17.84, four reported events, no severe initial DV flag) and TIC 89187590 (MES 14.83, three events, geometric-model S/N 7.79). Eight native SPOC observations totaling **18,095,040 bytes** were retrieved and integrity checked. This extends time-series follow-up from three to five stars and from nine to seventeen curves. The 846 public TCE metadata rows are separate from independent light-curve analysis counts.

**TIC 49617264 recovers a previously reported neighboring planet candidate, not a new object.** The approximately 4.832-day dips repeat in sectors 97 and 105–107, persist in raw SAP and corrected PDC, and survive quadratic local baselines and alternative gap-aware detrending. The restricted sector-107 fine search gives P≈4.832023 days, epoch BTJD≈4270.47597 and box duration≈2.29 hours; these are selected follow-up estimates, not a physical orbital posterior. A wider prior search then found original bright neighbor TIC 49617263 = TOI 2417.01 (FP) only 43.18 arcseconds away, outside the default 30-arcsecond match radius, at nearly the same period. This neighboring period alone does not identify an eclipsing binary: the decisive public 2024 follow-up notes already call the actual faint source TIC 49617264.11 / TOI 2417.11 and explicitly describe it as an NPC to the original bright candidate. The source was already reported and followed as a planetary candidate before this analysis.

The [Brierfield 2024-10-05 report](https://exofop.ipac.caltech.edu/tess/get_file.php?id=1004492) describes an on-time, on-target approximately 13.1 ppt event in a claimed uncontaminated 6.6-arcsecond aperture, with a reported 2h46m transit duration, and cites TRAPPIST 2024-09-15 zp follow-up at approximately 12.8 ppt and revised period 4.83261 days. It quotes '[NPC to original candidate TIC 49617263 (TOI-2417.01)] SAFFIs constrain event to this star relative to Gaia DR2 stars bright enough to cause the TESS detection.' These are prior observers' claims, not new confirmation by this project; the notes also contain transcription/date/telescope inconsistencies. The raw note, URL/hash manifest and context are preserved in `data/broad_catalogs/tic49617264_followup20241005_notes.txt` and `reports/tic49617264_prior_check.md`. An empty TOI field on the actual target and the original bright target's FP disposition must not be converted into either a novelty claim or an EB classification of the relocated candidate.

Crowding is substantial (CROWDSAP approximately 0.46–0.61), corrected depth changes between sectors, and independent pixel localization is useful as a check on the prior relocation. PDC linear-local depths are approximately 0.64%, 0.86%, 1.03%, 1.15% for sectors 97,105,106,107, while SAP depths are approximately 0.30%, 0.54%, 0.52%, 0.70%. Formal odd/even comparisons are weak in the three recent sectors; an older-sector formal difference is sensitive to one questionable event and must not be treated as a calibrated binary classification. A physical radius inferred by assuming the catalog target and its dilution remains conditional. The known candidate is a useful positive recovery and catalog-recognition lesson, without a new planetary designation or discovery claim.

**TIC 89187590 fails the proposed approximately 8.726-day planetary interpretation.** Every sector contains a strong approximately 1.09-day broad waveform. Quadratic native local baselines reduce sector-107's selected depth from approximately 1518 ppm to approximately 5 ppm, and sector-106's from approximately 278 ppm to approximately −616 ppm. The restricted fine BLS selects a different phase about 3.28 days from the public epoch at only diagonal S/N≈8.6. This is evidence of baseline curvature/variability rather than stable local transits. Earlier-sector fixed-phase negative depths alone are not used as decisive falsification because the public uncertainty columns are unreliable over seven years; recent native baselines already undermine the hypothesis. This does not exclude other planets around the star.

| TIC | Sector | Flux | Local linear depth ppm | Formal S/N | Event-scatter S/N | Events | Local quadratic depth ppm |
|---|---|---|---|---|---|---|---|
|49617264|97|PDCSAP_FLUX|6368.15|15.35|6.63|9|6833.42|
|49617264|97|SAP_FLUX|2975.92|14.95|14.40|9|2844.97|
|49617264|105|PDCSAP_FLUX|8622.38|18.80|11.45|5|8942.37|
|49617264|105|SAP_FLUX|5378.69|16.34|5.37|5|5743.67|
|49617264|106|PDCSAP_FLUX|10347.14|25.32|25.32|5|10522.95|
|49617264|106|SAP_FLUX|5212.59|23.98|23.98|5|5263.09|
|49617264|107|PDCSAP_FLUX|11499.98|21.90|19.01|4|11512.55|
|49617264|107|SAP_FLUX|7031.18|21.12|21.12|4|7059.36|
|89187590|4|PDCSAP_FLUX|-3622.21|-30.59|-6.82|3|-467.63|
|89187590|4|SAP_FLUX|-3565.87|-29.27|-5.98|3|-414.23|
|89187590|31|PDCSAP_FLUX|3621.23|29.05|1.93|3|87.88|
|89187590|31|SAP_FLUX|3621.17|27.94|1.91|3|53.18|
|89187590|106|PDCSAP_FLUX|277.77|1.94|1.94|2|-616.14|
|89187590|106|SAP_FLUX|-148.11|-1.04|-0.32|2|175.94|
|89187590|107|PDCSAP_FLUX|1518.23|12.50|2.39|3|5.21|
|89187590|107|SAP_FLUX|1559.86|12.95|2.49|3|58.86|

Reproduce with `python src/tce107_followup_extra.py`. If the cache is lost, use `--retrieve`; exact URLs, MAST identifiers, sectors, processing versions, SHA-256 and byte counts are frozen in `data/tce107_extra_observations.json`. The original nine-product manifest remains unchanged. Raw inputs are retained outside Git; diagnostic plots and native event-level measurements are in this directory. The wrapper checks hashes and raw FITS checksums before each analysis. These follow-up ranking statistics and simple covariance errors are not global false-alarm probabilities.

Next: preserve the known-candidate recovery and independent pixel comparison; focus discovery effort on other unrecognized signals. Future shortlisted targets need a match radius appropriate to their pixel aperture, nearby-star period checks, aliases and ground-follow-up uploads, even when the target's own catalog row is empty. This investigation found three variability/baseline aliases, one previously published complex periodic variable, and one previously reported neighboring planet candidate among the five independently followed public-TCE stars. No credible new planet emerged from this track.
