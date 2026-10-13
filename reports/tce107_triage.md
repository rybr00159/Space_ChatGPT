# Public sector-107 TCE triage

Checked 2026-10-13: **846 previously published SPOC TCEs on 623 unique stars**, not 623 new detections and not an independent 160k-star search. 71 rows meet the initial metadata screen. This is prioritization for independent time-series retrieval, not a credible candidate list.

Source: [https://archive.stsci.edu/missions/tess/catalogs/tce/tess2026221140108-s0107-s0107_dvr-tcestats.csv](https://archive.stsci.edu/missions/tess/catalogs/tce/tess2026221140108-s0107-s0107_dvr-tcestats.csv); created 2026-10-02; SHA-256 `94a19264c671a0d3b6c1857ba4870d596bd882e05517bd599977645189947382`. The pinned raw CSV is `data/broad_tce107.csv`. Reproduce with `python src/triage_tce107.py`.

The archive MES (`tce_max_mult_ev`) and geometric-model S/N are distinct. No official bootstrap FAP is treated as a survey-level probability. Crossmatching uses confirmed planets including RV detections, TOI, CTOI, DIAmante, explicit TESS EB membership and cached published candidates via `CatalogIndex`. TIC/period/harmonic matches are conservative ranking flags. This CSV contains no sky coordinates or CROWDSAP; position-based contamination screening and quantitative dilution remain unavailable until follow-up. A nonmatch is not proof of novelty. Literature and alias checks remain pending.

Screen: pipeline MES≥9; at least three reported transits; pipeline radius 0.3–20 Earth radii; depth 30–20,000 ppm; duration 0.5–8 hours; stellar radius 0.1–3 solar; Teff<7500 K; converged geometric fit with positive model S/N (excluding missing/failed −1 sentinels); no cached known-period/harmonic or explicit EB host; no likely secondary of another same-star TCE. The ranking uses log(1+MES), with explicit penalties for odd/even, weak-secondary and centroid statistics, very large radius, long duration, Solar placeholders and duplicate-period TCEs. These choices are transparent heuristics, not probabilities, completeness estimates or physical validation.

**Within-star TCE comparison matters:** TIC 146435300's appealing third row has the same period as a deeper first row and is displaced by approximately half an orbit; it is a likely secondary eclipse. TIC 178288581 similarly has a larger-than-planet first TCE plus a weaker half-phase row. Both would be misleading if only one individually planet-sized row were selected. These flags and all 846 rows remain searchable in the output CSV.

**Export caution:** many geometric-fit uncertainty columns have the same value in all 846 rows, while trapezoidal uncertainties are uniformly −1. Do not use these as trustworthy target-specific errors. Geometric transit epochs are rounded (often to 0.01 day). The table supplies full-precision `tcet_period` / `tcet_time0bt` as follow-up starting points, with the geometric values retained; both require native-FITS refitting.

Constant uncertainty columns: `tce_period_err`, `tce_time0bt_err`, `tce_time0_err`, `tce_ror_err`, `tce_dor_err`, `tce_incl_err`, `tce_impact_err`, `tce_duration_err`, `tce_ingress_err`, `tce_depth_err`, `tce_eccen_err`, `tce_longp_err`, `tce_prad_err`, `tce_sma_err`, `tce_eqt_err`, `tce_insol_err`, `tce_ditco_jra_err`, `tce_ditco_jdec_err`, `tce_ditco_jsky_err`.

| tic | tceid | period_days_followup | epoch_btjd_followup | pipeline_mes | tce_num_transits | tce_prad | tce_depth | tce_duration | diagnostic_flags |
|---|---|---|---|---|---|---|---|---|---|
| 72729777 | 00072729777-01 | 7.98610734939575 | 4266.87977897325 | 282.9372 | 3 | 4.51866 | 1439.98 | 5.63856 | centroid_unavailable;geometric_fit_SNR_below_9 |
| 159923666 | 00159923666-01 | 8.07499599456787 | 4267.47530361319 | 26.16343 | 3 | 6.66143 | 3597.18 | 5.12775 |  |
| 178155030 | 00178155030-01 | 0.4865727424621579 | 4262.139708944929 | 24.73565 | 42 | 5.94182 | 15296.2 | 1.09404 |  |
| 100777592 | 00100777592-01 | 5.84583044052124 | 4262.32385474538 | 22.06946 | 3 | 7.51538 | 13814.4 | 0.8269149999999998 |  |
| 77437461 | 00077437461-01 | 0.607908666133881 | 4262.67014129924 | 72.40169 | 34 | 5.92589 | 14566.4 | 0.506768 | DV_centroid_offset |
| 44790652 | 00044790652-01 | 7.48888540267944 | 4267.54598430491 | 21.51999 | 4 | 6.9976 | 11631.2 | 2.60906 | weak_secondary_DV_MES;centroid_unavailable |
| 139827347 | 00139827347-01 | 5.96944141387939 | 4266.03702977374 | 21.26053 | 4 | 4.11103 | 958.463 | 7.56228 | long_duration_for_circular_central_transit;geometric_fit_SNR_below_9 |
| 89187590 | 00089187590-01 | 8.72638511657715 | 4266.58384602016 | 14.82587 | 3 | 4.25286 | 2108.42 | 1.18092 | geometric_fit_SNR_below_9 |
| 49617264 | 00049617264-01 | 4.83194208145142 | 4265.64815338287 | 17.84495 | 4 | 13.1446 | 14601.6 | 2.77069 |  |
| 47596872 | 00047596872-01 | 0.208089023828506 | 4262.19566585068 | 12.71762 | 99 | 5.42966 | 12985.1 | 0.687025 |  |
| 114851071 | 00114851071-01 | 8.74860668182373 | 4266.08348257439 | 12.10056 | 3 | 3.41801 | 4582.3 | 1.87357 |  |
| 175491080 | 00175491080-01 | 8.13749599456787 | 4268.11070964651 | 12.04215 | 3 | 6.40384 | 3627.35 | 1.5655 |  |
| 92846244 | 00092846244-01 | 0.256125777959824 | 4262.26901925819 | 11.65198 | 87 | 6.86188 | 14964.1 | 1.04681 |  |
| 168849065 | 00168849065-01 | 8.93610668182373 | 4266.68888069036 | 16.77543 | 3 | 2.388 | 6668.17 | 3.07313 | centroid_unavailable;long_duration_for_circular_central_transit |
| 229080275 | 00229080275-01 | 6.6347188949585 | 4266.83539520171 | 11.17813 | 4 | 3.85188 | 657.887 | 1.95831 | geometric_fit_SNR_below_9 |
| 594723 | 00000594723-01 | 2.0513880252838095 | 4263.80378030715 | 23.80437 | 12 | 15.5692 | 16534.4 | 2.85378 |  |
| 123727832 | 00123727832-02 | 8.72499582905766 | 4264.48424535481 | 10.51326 | 3 | 6.54055 | 5972.96 | 2.38471 | centroid_unavailable;geometric_fit_SNR_below_9 |
| 159973423 | 00159973423-01 | 5.73055267333984 | 4266.4302134559 | 12.49808 | 3 | 5.86318 | 1841.19 | 1.15846 | weak_secondary_DV_MES |
| 179106508 | 00179106508-01 | 0.370439827442169 | 4262.21681886624 | 10.24004 | 61 | 8.11739 | 8518.59 | 1.58368 |  |
| 257644579 | 00257644579-01 | 8.76388454437256 | 4263.716622076209 | 17.63603 | 3 | 7.67305 | 3218.07 | 1.50257 | weak_secondary_DV_MES |

Next: independently retrieve native SPOC sector-107 light curves for the highest-ranked clean targets; measure events in SAP and PDCSAP with local baselines, assess individual recurrence, refit epoch/period/depth/duration and uncertainties, query prior sectors, then inspect pixels/centroids and catalogs/literature for any survivor. Public TCE recognition must remain explicit even if a planetary interpretation later survives. No new planet is claimed from this table.
