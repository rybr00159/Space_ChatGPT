J/other/PASA/42.51               GALAH survey DR4                 (Buder+, 2025)
================================================================================
The GALAH survey: Data release 4.
    Buder S., Kos J., Wang X.E., McKenzie M., Howell M., Martell S.,
    Hayden M.R., Zucker D.B., Nordlander T., Montet B., Traven G.,
    Bland-Hawthorn J., De Silva G.M., Freeman K., Lewis G., Lind K., Sharma S.,
    Simpson J.D., Stello D., Zwitter T., Amarsi A.M., Armstrong J.J., Banks K.,
    Beavis M., Beeson K.-L., Chen B., Ciuca I., Da Costa G.S., de Grijs R.,
    Martin B., Nataf D.M., Ness M., Rains A.D., Scarr T., Vogrincic R.,
    Purmortal Wang Z., Wittenmyer R.A., Xie Y.A., GALAH Coll.
    <Publ. Astron. Soc. Australia, 42, 51 (2025)>
    =2025PASA...42...51B        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Surveys ; Milky Way ; Photometry, infrared ; Abundances
Keywords: surveys - the Galaxy - methods: observational -
          methods: data analysis - stars: fundamental parameters -
          stars: abundances

Abstract:
    The stars of the Milky Way carry the chemical history of our Galaxy in
    their atmospheres as they journey through its vast expanse. Like
    barcodes, we can extract the chemical fingerprints of stars from
    high-resolution spectroscopy. The fourth data release (DR4) of the
    Galactic Archaeology with HERMES (GALAH) Survey, based on a decade of
    observations, provides the chemical abundances of up to 32 elements
    for 917 588 stars that also have exquisite astrometric data from the
    Gaia satellite. For the first time, these elements include
    life-essential nitrogen to complement carbon, and oxygen as well as
    more measurements of rare-earth elements critical to modern-life
    electronics, offering unparalleled insights into the chemical
    composition of the Milky Way. For this release, we use neural networks
    to simultaneously fit stellar parameters and abundances across the
    whole wavelength range, leveraging synthetic grids computed with
    Spectroscopy Made Easy. These grids account for atomic line formation
    in non-local thermodynamic equilibrium for 14 elements. In a
    two-iteration process, we first fit stellar labels to all 1 085 520
    spectra, then co-add repeated observations and refine these labels
    using astrometric data from Gaia and 2MASS photometry, improving the
    accuracy and precision of stellar parameters and abundances. Our
    validation thoroughly assesses the reliability of spectroscopic
    measurements and highlights key caveats. GALAH DR4 represents yet
    another milestone in Galactic archaeology, combining detailed chemical
    compositions from multiple nucleosynthetic channels with kinematic
    information and age estimates. The resulting dataset, covering nearly
    a million stars, opens new avenues for understanding not only the
    chemical and dynamical history of the Milky Way but also the broader
    questions of the origin of elements and the evolution of planets,
    stars, and galaxies.

Description:
    The GALAH Survey uses the 3.9-m Anglo-Australian Telescope at Siding
    Spring Observatory on Gamilaraay Country and its Two-Degree Field
    positioning system (2dF) top end (Lewis et al., 2002MNRAS.333..279L).
    2dF magnetically places up to 400 fibre buttons on one of two metal
    field plates, which can be tumbled to allow observing with one set of
    fibres while configuring the other. Light is delivered through the
    fibres to the High Efficiency and Resolution Multi-Element
    Spectrograph (HERMES) spectrograph (Barden et al., 2010, SPIE, 7735,
    09; Brzeski Case & Gers 2011, SPIE, 8125, 04; Heijmans et al. 2012,
    SPIE, 8446, 0W; Farrell et al. 2014, SPIE, 9152, 23; Sheinis et al.
    2015, JATIS, 1, 035002 and dispersed into four non-contiguous
    wavelength bands in the optical that cover ~1000{AA} in the range of
    4713-4903 (blue CCD or CCD1), 5648-5873 (green/CCD2), 6478-6737
    (red/CCD3), and 7585-7887{AA} (infrared IR/CCD4). The data used in
    this data release is primarily based on observations of stars with
    this setup, but also makes use of auxiliary photometric and
    astrometric information for the stars where available.

    allstar.dat is our main results catalogue. It contains results for
    917,588 stars observed as part of the GALAH, K2-HERMES, TESS-HERMES,
    and other related surveys that used the HERMES spectrograph on the
    Anglo-Australian Telescope between December 2013 and August 2023. For
    all targets we provide stellar parameters, radial velocities, and
    elemental abundances.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
allstar.dat     1756   917588   GALAH DR4 all stars catalog (24-07-05 version)
allspec.dat     1721  1085520   GALAH DR4 all spectra catalog (24-07-05 version)
--------------------------------------------------------------------------------

See also:
   https://galah-survey.org/ : GALAH Home Page


Byte-by-byte Description file: allstar.dat
--------------------------------------------------------------------------------
   Bytes   Format Units   Label           Explanations
--------------------------------------------------------------------------------
    1-  15  I15   ---       GALAH         GALAH identifier (sobject_id)
   17-  33  A17   ---       2MASS         2MASS identifier (tmass_id)
   35-  53  I19   ---       GaiaDR3       Gaia DR3 source_id (gaiadr3_source_id)
   55-  71  A17   ---       SurveyName    HERMES-2dF Survey/Program/Pointing
                                           (other if not available)
                                           (survey_name)
   73-  83  I11   ---       FieldId       ?=-1 GALAH Field ID
                                           (-1 if not available) (field_id)
   85-  91  A7    ---       setup         Analysis setup: single/binary/coadds
                                           (setup)
   93- 104  F12.6 d         MJD           Modified Julian Date (mjd)
  106- 120 F15.11 deg       RAdeg         Right ascension (J2000) propagated
                                           from Gaia DR3 (ra)
  122- 136 F15.11 deg       DEdeg         Declination (J2000) propagated
                                           from Gaia DR3 (dec)
  138- 142  I5    ---       flagSp        Major spectroscopic quality
                                           bitmask flag (flag_sp) (G1)
  144- 145  I2    ---       flagSpFit     Major spectroscopic fitting
                                           quality bitmask flag (flag_sp_fit)
  147- 148  I2    ---       Optloop       Nr of optimisation loops used for
                                           fitting (opt_loop)
  150- 156  I7    ---       flagRed       Quality bitmask flag of reduction
                                           pipeline (flag_red)
  158- 168  F11.6 ---       SNRpxCCD1     ?=- Average signal-to-noise ratio
                                           (per pixel) of CCD1 (snr_px_ccd1)
  170- 180  F11.6 ---       SNRpxCCD2     Average signal-to-noise ratio
                                           (per pixel) of CCD2 (snr_px_ccd2)
  182- 193  F12.6 ---       SNRpxCCD3     ?=- Average signal-to-noise ratio
                                           (per pixel) of CCD3 (snr_px_ccd3)
  195- 206  F12.6 ---       SNRpxCCD4     ?=- Average signal-to-noise ratio
                                           (per pixel) of CCD4 (snr_px_ccd4)
  208- 219  E12.6 ---       chi2sp        ?=- Chi2 value of spectroscopic
                                           fitting (chi2_sp)
  221- 223  I3    %         pxUsedPerc    Percentage of spectrum used for
                                           spectroscopic fit (px_used_perc)
  225- 239  A15   ---       ModelName     Neural network model used for creating
                                           synthetic spectra (model_name)
  241- 255  A15   ---       ClosestModel  Neural network model closest for
                                           Teff/logg/[Fe/H] combination
                                           (closest_model)
  257- 268  F12.6 s         compTime      ?=- Computation time spent on spectrum
                                           (comp_time)
       270  I1    ---       fitGlobalRV   Are spectra coadded (without RV
                                           shifts) and then a global RV fitted?
                                           (fit_global_rv)
  272- 283  F12.6 km/s      RVcomp1       ?=- Radial velocity of primary source
                                           (rv_comp_1)
  285- 296  E12.6 km/s    e_RVcomp1       ?=- Uncertainty on RVcomp1
                                           (e_rv_comp_1)
  298- 309  F12.6 km/s      RVcomp2       ?=- Radial velocity of potential
                                           secondary source (rv_comp_2)
  311- 320  F10.6 km/s    e_RVcomp2       ?=- Uncertainty on RVcomp1
                                           (e_rv_comp_2)
  322- 333  F12.6 km/s      RVGDR3        ?=- Radial velocity in Gaia DR3
                                           (rv_gaia_dr3)
  335- 343  F9.6  km/s    e_RVGDR3        ?=- Uncertainty on rv_gaia_dr3
                                           (e_rv_gaia_dr3)
  345- 354  F10.6 km/s      VbaryEff      ?=- Barycentric velocity applied to
                                           reduced spectra (v_bary_eff)
  356- 362  F7.2  K         Teff          ?=- Spectroscopic effective
                                           temperature (used for fitting) (teff)
  364- 375  E12.6 K       e_Teff          ?=- Uncertainty on teff (e_teff)
  377- 385  F9.6  [cm/s2]   logg          ?=- Surface gravity adjusted via
                                           parallax information (logg)
  387- 394  F8.6  [cm/s2] e_logg          ?=- Uncertainty on logg (e_logg)
  396- 404  F9.6  ---       [Fe/H]        ?=- Abundance of Fe and all other
                                           elements not fitted in GALAH
                                           (Fe: 1D-NLTE) (fe_h)
  406- 417  E12.6 ---     e_[Fe/H]        ?=- Uncertainty on [Fe/H] (e_fe_h)
       419  I1    ---     f_[Fe/H]        Quality flag on [Fe/H]
                                           (flag_fe_h) (G2)
  421- 428  F8.6  km/s      vmic          ?=- Microturbulence velocity (fitted)
                                           (vmic)
  430- 442  F13.6 km/s    e_vmic          ?=- Uncertainty on vmic (e_vmic)
  444- 451  F8.5  km/s      vsini         ?=- Broadening velocity
                                           (fitted sme.vsini with sme.vmac=0)
                                           (vsini)
  453- 465  F13.6 km/s    e_vsini         ?=- Uncertainty on vsini (e_vsini)
  467- 475  F9.6  ---       [Li/Fe]       ?=- Elemental abundance for [Li/Fe]
                                          (li_fe)
  477- 488  E12.6 ---     e_[Li/Fe]       ?=- Uncertainty on [Li/Fe] (e_li_fe)
  490- 491  I2    ---     f_[Li/Fe]       Quality bitmask flag of [Li/Fe]
                                           (flag_li_fe) (G2)
  493- 501  F9.6  ---       [C/Fe]        ?=- Elemental abundance for [C/Fe]
                                           (c_fe)
  503- 514  E12.6 ---     e_[C/Fe]        ?=- Uncertainty on [C/Fe] (e_c_fe)
  516- 517  I2    ---     f_[C/Fe]        Quality bitmask flag of [C/Fe]
                                           (flag_c_fe) (G2)
  519- 527  F9.6  ---       [N/Fe]        ?=- Elemental abundance for [N/Fe]
                                           (n_fe)
  529- 537  F9.6  ---     e_[N/Fe]        ?=- Uncertainty on [N/Fe] (e_n_fe)
  539- 541  I3    ---     f_[N/Fe]        Quality bitmask flag of [N/Fe]
                                           (flag_n_fe) (G2)
  543- 551  F9.6  ---       [O/Fe]        ?=- Elemental abundance for [O/Fe]
                                           (o_fe)
  553- 564  E12.6 ---     e_[O/Fe]        ?=- Uncertainty on [O/Fe] (e_o_fe)
  566- 567  I2    ---     f_[O/Fe]        Quality bitmask flag of [O/Fe]
                                           (flag_o_fe) (G2)
  569- 577  F9.6  ---       [Na/Fe]       ?=- Elemental abundance for [Na/Fe]
                                           (na_fe)
  579- 588  F10.6 ---     e_[Na/Fe]       ?=- Uncertainty on [Na/Fe] (e_na_fe)
  590- 591  I2    ---     f_[Na/Fe]       Quality bitmask flag of [Na/Fe]
                                           (flag_na_fe) (G2)
  593- 601  F9.6  ---       [Mg/Fe]       ?=- Elemental abundance for [Mg/Fe]
                                           (mg_fe)
  603- 614  E12.6 ---     e_[Mg/Fe]       ?=- Uncertainty on [Mg/Fe] (e_mg_fe)
  616- 617  I2    ---     f_[Mg/Fe]       Quality bitmask flag of [Mg/Fe]
                                           (flag_mg_fe) (G2)
  619- 627  F9.6  ---       [Al/Fe]       ?=- Elemental abundance for [Al/Fe]
                                           (al_fe)
  629- 637  F9.6  ---     e_[Al/Fe]       ?=- Uncertainty on [Al/Fe] (e_al_fe)
  639- 640  I2    ---     f_[Al/Fe]       Quality bitmask flag of [Al/Fe]
                                           (flag_al_fe) (G2)
  642- 650  F9.6  ---       [Si/Fe]       ?=- Elemental abundance for [Si/Fe]
                                           (si_fe)
  652- 660  F9.6  ---     e_[Si/Fe]       ?=- Uncertainty on [Si/Fe] (e_si_fe)
  662- 663  I2    ---     f_[Si/Fe]       Quality bitmask flag of [Si/Fe]
                                           (flag_si_fe) (G2)
  665- 673  F9.6  ---       [K/Fe]        ?=- Elemental abundance for [K/Fe]
                                           (k_fe)
  675- 683  F9.6  ---     e_[K/Fe]        ?=- Uncertainty on [K/Fe] (e_k_fe)
  685- 686  I2    ---     f_[K/Fe]        Quality bitmask flag of [K/Fe]
                                           (flag_k_fe) (G2)
  688- 696  F9.6  ---       [Ca/Fe]       ?=- Elemental abundance for [Ca/Fe]
                                           (ca_fe)
  698- 709  E12.6 ---     e_[Ca/Fe]       ?=- Uncertainty on [Ca/Fe] (e_ca_fe)
  711- 712  I2    ---     f_[Ca/Fe]       Quality bitmask flag of [Ca/Fe]
                                           (flag_ca_fe) (G2)
  714- 722  F9.6  ---       [Sc/Fe]       ?=- Elemental abundance for [Sc/Fe]
                                           (sc_fe)
  724- 735  E12.6 ---     e_[Sc/Fe]       ?=- Uncertainty on [Sc/Fe] (e_sc_fe)
  737- 738  I2    ---     f_[Sc/Fe]       Quality bitmask flag of [Sc/Fe]
                                           (flag_sc_fe) (G2)
  740- 748  F9.6  ---       [Ti/Fe]       ?=- Elemental abundance for [Ti/Fe]
                                           (ti_fe)
  750- 761  E12.6 ---     e_[Ti/Fe]       ?=- Uncertainty on [Ti/Fe] (e_ti_fe)
  763- 764  I2    ---     f_[Ti/Fe]       Quality bitmask flag of [Ti/Fe]
                                           (flag_ti_fe) (G2)
  766- 774  F9.6  ---       [V/Fe]        ?=- Elemental abundance for [V/Fe]
                                           (v_fe)
  776- 785  F10.6 ---     e_[V/Fe]        ?=- Uncertainty on [V/Fe] (e_v_fe)
  787- 788  I2    ---     f_[V/Fe]        Quality bitmask flag of [V/Fe]
                                           (flag_v_fe) (G2)
  790- 798  F9.6  ---       [Cr/Fe]       ?=- Elemental abundance for [Cr/Fe]
                                           (cr_fe)
  800- 811  E12.6 ---     e_[Cr/Fe]       ?=- Uncertainty on [Cr/Fe] (e_cr_fe)
  813- 814  I2    ---     f_[Cr/Fe]       Quality bitmask flag of [Cr/Fe]
                                           (flag_cr_fe) (G2)
  816- 824  F9.6  ---       [Mn/Fe]       ?=- Elemental abundance for [Mn/Fe]
                                           (mn_fe)
  826- 835  F10.6 ---     e_[Mn/Fe]       ?=- Uncertainty on [Mn/Fe] (e_mn_fe)
  837- 838  I2    ---     f_[Mn/Fe]       Quality bitmask flag of [Mn/Fe]
                                           (flag_mn_fe) (G2)
  840- 848  F9.6  ---       [Co/Fe]       ?=- Elemental abundance for [Co/Fe]
                                           (co_fe)
  850- 859  F10.6 ---     e_[Co/Fe]       ?=- Uncertainty on [Co/Fe] (e_co_fe)
  861- 862  I2    ---     f_[Co/Fe]       Quality bitmask flag of [Co/Fe]
                                           (flag_co_fe) (G2)
  864- 872  F9.6  ---       [Ni/Fe]       ?=- Elemental abundance for [Ni/Fe]
                                           (ni_fe)
  874- 885  E12.6 ---     e_[Ni/Fe]       ?=- Uncertainty on [Ni/Fe] (e_ni_fe)
  887- 888  I2    ---     f_[Ni/Fe]       Quality bitmask flag of [Ni/Fe]
                                           (flag_ni_fe) (G2)
  890- 898  F9.6  ---       [Cu/Fe]       ?=- Elemental abundance for [Cu/Fe]
                                           (cu_fe)
  900- 909  F10.6 ---     e_[Cu/Fe]       ?=- Uncertainty on [Cu/Fe] (e_cu_fe)
  911- 912  I2    ---     f_[Cu/Fe]       Quality bitmask flag of [Cu/Fe]
                                           (flag_cu_fe) (G2)
  914- 922  F9.6  ---       [Zn/Fe]       ?=- Elemental abundance for [Zn/Fe]
                                           (zn_fe)
  924- 933  F10.6 ---     e_[Zn/Fe]       ?=- Uncertainty on [Zn/Fe] (e_zn_fe)
  935- 936  I2    ---     f_[Zn/Fe]       Quality bitmask flag of [Zn/Fe]
                                           (flag_zn_fe) (G2)
  938- 946  F9.6  ---       [Rb/Fe]       ?=- Elemental abundance for [Rb/Fe]
                                           (rb_fe)
  948- 956  F9.6  ---     e_[Rb/Fe]       ?=- Uncertainty on [Rb/Fe] (e_rb_fe)
  958- 959  I2    ---     f_[Rb/Fe]       Quality bitmask flag of [Rb/Fe]
                                           (flag_rb_fe) (G2)
  961- 969  F9.6  ---       [Sr/Fe]       ?=- Elemental abundance for [Sr/Fe]
                                           (sr_fe)
  971- 979  F9.6  ---     e_[Sr/Fe]       ?=- Uncertainty on [Sr/Fe] (e_sr_fe)
  981- 982  I2    ---     f_[Sr/Fe]       Quality bitmask flag of [Sr/Fe]
                                           (flag_sr_fe) (G2)
  984- 992  F9.6  ---       [Y/Fe]        ?=- Elemental abundance for [Y/Fe]
                                           (y_fe)
  994-1005  E12.6 ---     e_[Y/Fe]        ?=- Uncertainty on [Y/Fe] (e_y_fe)
 1007-1008  I2    ---     f_[Y/Fe]        Quality bitmask flag of [Y/Fe]
                                           (flag_y_fe) (G2)
 1010-1018  F9.6  ---       [Zr/Fe]       ?=- Elemental abundance for [Zr/Fe]
                                           (zr_fe)
 1020-1029  F10.6 ---     e_[Zr/Fe]       ?=- Uncertainty on [Zr/Fe] (e_zr_fe)
 1031-1032  I2    ---     f_[Zr/Fe]       Quality bitmask flag of [Zr/Fe]
                                           (flag_zr_fe) (G2)
 1034-1042  F9.6  ---       [Mo/Fe]       ?=- Elemental abundance for [Mo/Fe]
                                           (mo_fe)
 1044-1052  F9.6  ---     e_[Mo/Fe]       ?=- Uncertainty on [Mo/Fe] (e_mo_fe)
      1054  I1    ---     f_[Mo/Fe]       Quality bitmask flag of [Mo/Fe]
                                           (flag_mo_fe) (G2)
 1056-1064  F9.6  ---       [Ru/Fe]       ?=- Elemental abundance for [Ru/Fe]
                                           (ru_fe)
 1066-1074  F9.6  ---     e_[Ru/Fe]       ?=- Uncertainty on [Ru/Fe] (e_ru_fe)
      1076  I1    ---     f_[Ru/Fe]       Quality bitmask flag of [Ru/Fe]
                                           (flag_ru_fe) (G2)
 1078-1086  F9.6  ---       [Ba/Fe]       ?=- Elemental abundance for [Ba/Fe]
                                           (ba_fe)
 1088-1099  E12.6 ---     e_[Ba/Fe]       ?=- Uncertainty on [Ba/Fe] (e_ba_fe)
 1101-1102  I2    ---     f_[Ba/Fe]       Quality bitmask flag of [Ba/Fe]
                                           (flag_ba_fe) (G2)
 1104-1112  F9.6  ---       [La/Fe]       ?=- Elemental abundance for [La/Fe]
                                           (la_fe)
 1114-1123  F10.6 ---     e_[La/Fe]       ?=- Uncertainty on [La/Fe] (e_la_fe)
 1125-1126  I2    ---     f_[La/Fe]       Quality bitmask flag of [La/Fe]
                                           (flag_la_fe) (G2)
 1128-1136  F9.6  ---       [Ce/Fe]       ?=- Elemental abundance for [Ce/Fe]
                                           (ce_fe)
 1138-1147  F10.6 ---     e_[Ce/Fe]       ?=- Uncertainty on [Ce/Fe] (e_ce_fe)
 1149-1150  I2    ---     f_[Ce/Fe]       Quality bitmask flag of [Ce/Fe]
                                           (flag_ce_fe) (G2)
 1152-1160  F9.6  ---       [Nd/Fe]       ?=- Elemental abundance for [Nd/Fe]
                                           (nd_fe)
 1162-1171  F10.6 ---     e_[Nd/Fe]       ?=- Uncertainty on [Nd/Fe] (e_nd_fe)
 1173-1174  I2    ---     f_[Nd/Fe]       Quality bitmask flag of [Nd/Fe]
                                           (flag_nd_fe) (G2)
 1176-1184  F9.6  ---       [Sm/Fe]       ?=- Elemental abundance for [Sm/Fe]
                                           (sm_fe)
 1186-1194  F9.6  ---     e_[Sm/Fe]       ?=- Uncertainty on [Sm/Fe] (e_sm_fe)
 1196-1197  I2    ---     f_[Sm/Fe]       Quality bitmask flag of [Sm/Fe]
                                           (flag_sm_fe) (G2)
 1199-1207  F9.6  ---       [Eu/Fe]       ?=- Elemental abundance for [Eu/Fe]
                                           (eu_fe)
 1209-1216  F8.6  ---     e_[Eu/Fe]       ?=- Uncertainty on [Eu/Fe] (e_eu_fe)
 1218-1219  I2    ---     f_[Eu/Fe]       Quality bitmask flag of [Eu/Fe]
                                           (flag_eu_fe) (G2)
 1221-1228  F8.6  Msun      Mass          ?=- Mass used for calculating
                                           logg (plx) (mass)
 1230-1238  F9.6  Gyr       Age           ?=- Age estimated when calculating
                                           mass (age)
 1240-1247  F8.6  mag       BCKs          ?=- Bolometric Correction
                                           Ks, BC(Ks), used for calculating
                                           logg (plx) (bc_ks)
 1249-1256  F8.6  mag       AKs           ?=- Attenuation in Ks-band A(Ks) used
                                           for calculating logg (plx) (a_ks)
 1258-1269  E12.6 Lsun      Lbol          ?=- Bolometric Luminosity used for
                                           calculating logg (plx) (lbol)
 1271-1282  E12.6 pc        rmed          ?=- Median Distance used for
                                           calculating logg (plx) (r_med)
 1284-1295  E12.6 pc      b_rmed          ?=- Lower Limit Distance used for
                                           calculating logg (plx) (r_lo)
 1297-1308  E12.6 pc      B_rmed          ?=- Higher Limit Distance used for
                                           calculating logg (plx) (r_hi)
 1310-1321  E12.6 km/s      sb2RV16       ?=- 16th perc. radial velocity of fit
                                           to syn-obs residuals (sb2_rv_16)
 1323-1334  E12.6 km/s      sb2RV50       ?=- 50th perc. radial velocity of fit
                                           to syn-obs residuals (sb2_rv_50)
 1336-1347  E12.6 km/s      sb2RV84       ?=- 84th perc. radial velocity of fit
                                           to syn-obs residuals (sb2_rv_84)
 1349-1359  F11.6 0.1nm     EWHbeta       ?=- Equivalent Width of observed Hbeta
                                           core (ew_h_beta)
 1361-1372  F12.6 0.1nm     EWHalpha      ?=- Equivalent Width of observed
                                           Halpha core (ew_h_alpha)
 1374-1384  F11.6 0.1nm     ResHbeta      ?=- Residual Equivalent Width of
                                           syn-obs residuals at Hbeta core
                                           (res_h_beta)
 1386-1397  F12.6 0.1nm     ResHalpha     ?=- Residual Equivalent Width of
                                           syn-obs residuals at Halpha core
                                           (res_h_alpha)
 1399-1406  F8.6  0.1nm     EWKis         ?=- Equivalent Width of fit for K7699
                                           Interstellar Line (ew_k_is)
 1408-1415  F8.6  ---       sigmaKis      ?=- Sigma of Gaussian fit for K7699
                                           Interstellar Line (sigma_k_is)
 1417-1427  F11.6 km/s      RVKis         ?=- Radial velocity of fit to syn-obs
                                           residuals around K7699 line (rv_k_is)
 1429-1437  F9.6  0.1nm     EWdib5780     ?=- Equivalent Width of fit for 5780
                                           Diffiuse Interstellar Band
                                           (ew_dib5780)
 1439-1446  F8.6  ---       sigdib5780    ?=- Sigma of Gaussian fit for 5780
                                           DIB (sigma_dib5780)
 1448-1459  F12.6 km/s      RVdib5780     ?=- Radial velocity of fit to syn-obs
                                           residuals around 5780 DIB
                                           (rv_dib5780)
 1461-1469  F9.6  0.1nm     EWdib5797     ?=- Equivalent Width of fit for 5797
                                           Diffiuse Interstellar Band
                                           (ew_dib5797)
 1471-1478  F8.6  ---       sigdib5797    ?=- Sigma of Gaussian fit for 5797 DIB
                                           (sigma_dib5797)
 1480-1491  F12.6 km/s      RVdib5797     ?=- Radial velocity of fit to syn-obs
                                           residuals around 5797 DIB
                                           (rv_dib5797)
 1493-1501  F9.6  0.1nm     EWdib6613     ?=- Equivalent Width of fit for
                                           6613 Diffiuse Interstellar Band
                                           (ew_dib6613)
 1503-1510  F8.6  ---       sigdib6613    ?=- Sigma of Gaussian fit for
                                           6613 DIB (sigma_dib6613)
 1512-1523  F12.6 km/s      RVdib6613     ?=- Radial velocity of fit to syn-obs
                                           residuals around 6613 DIB
                                           (rv_dib6613)
 1525-1534  F10.6 mag       E(B-V)        ?=- Extinction E(B-V) (ebv)
 1536-1543  F8.5  mag       Gmag          ?=- Mean Gaia DR3 G-band apparent
                                           magnitude (phot_g_mean_mag)
 1545-1553  F9.6  mag       BP-RP         ?=- Color of BP-RP bands (bp_rp)
 1555-1561  F7.3  mag       Jmag          ?=- 2MASS J-band magnitude (j_m)
 1563-1574  E12.6 mag     e_Jmag          ?=- Uncertainty on Jmag (j_msigcom)
 1576-1582  F7.3  mag       Hmag          ?=- 2MASS H-band magnitude (h_m)
 1584-1595  E12.6 mag     e_Hmag          ?=- Uncertainty on Hmag (h_msigcom)
 1597-1603  F7.3  mag       Ksmag         ?=- 2MASS Ks-band magnitude (ks_m)
 1605-1616  E12.6 mag     e_Ksmag         ?=- Uncertainty on Ksmag (ks_msigcom)
 1618-1624  F7.3  mag       W2mag         ?=- AllWISE W2-band magnitude (W2mag)
 1626-1637  E12.6 mag     e_W2mag         ?=- uncertainty of W2mag (e_W2mag)
 1639-1647  F9.6  ---       RUWE          ?=- RUWE reported by Gaia DR3 (ruwe)
 1649-1658  F10.6 mas       plx           ?=- Astrometric parallax used for
                                           GALAH DR4 (parallax)
 1660-1667  F8.6  mas     e_plx           ?=- Uncertainty on astrometric
                                           parallax used for GALAH DR4
                                           (parallax_error)
 1669-1680  F12.6 0.1nm     EWLi          ?=- Equivalent width of LiI 6708 line
                                           region (ew_li)
 1682-1693  F12.6 0.1nm   e_EWLi          ?=- Lower uncertainty limit of
                                           ew_li (e_ew_li_low)
 1695-1706  F12.6 0.1nm   E_EWLi          ?=- Upper uncertainty limit of
                                           ew_li (e_ew_li_upp)
 1708-1716  F9.6  ---       ALi           ?=- Absolute Lithium abundance (a_li)
 1718-1726  F9.6  ---       ALiul         ?=- Upper limit of a_li (a_li_upp_lim)
 1728-1735  F8.6  ---     e_ALi           ?=- Lower uncertainty of a_li
                                           (e_a_li_low)
 1737-1744  F8.6  ---     E_ALi           ?=- Upper uncertainty of a_li
                                           (e_a_li_upp)
 1746-1753  F8.6  ---     e_ALiTeff       ?=- Teff-correlated uncertainty of
                                           a_li (e_a_li_teff)
 1755-1756  I2    ---     f_ALi           ?=- Quality bitmask flag of a_li
                                           (flag_a_li) (G2)
--------------------------------------------------------------------------------

Byte-by-byte Description of file: allspec.dat
--------------------------------------------------------------------------------
   Bytes   Format Units     Label        Explanations
--------------------------------------------------------------------------------
    1-  15  I15   ---       GALAH        GALAH identifier (sobject_id)
   17-  33  A17   ---       2MASS        2MASS identifier (tmass_id)
   35-  53  I19   ---       GaiaDR3      Gaia DR3 source_id (gaiadr3_source_id)
   55-  71  A17   ---       SurveyName   HERMES-2dF Survey/Program/Pointing
                                          (other if not available) (survey_name)
   73-  83  I11   ---       FieldId      ?=-1 GALAH Field ID
                                          (-1 if not available) (field_id)
   85-  90  A6    ---       setup        Analysis setup: single/binary/coadds
                                          (setup)
   92- 103  F12.6 d         MJD          Modified Julian Date (mjd)
  105- 119 F15.11 deg       RAdeg        Right ascension (J2000) propagated
                                          from Gaia DR3 (ra)
  121- 135 F15.11 deg       DEdeg        Declination (J2000) propagated
                                          from Gaia DR3 (dec)
  137- 141  I5    ---       flagSp       Major spectroscopic quality
                                          bitmask flag (flag_sp) (G1)
  143- 144  I2    ---       flagSpFit    Major spectroscopic fitting
                                          quality bitmask flag (flag_sp_fit)
  146- 147  I2    ---       Optloop      Nr of optimisation loops used for
                                          fitting (opt_loop)
  149- 155  I7    ---       flagRed      Quality bitmask flag of reduction
                                          pipeline (flag_red)
  157- 167  F11.6 ---       SNRpxCCD1    ?=- Average signal-to-noise ratio
                                          (per pixel) of CCD1 (snr_px_ccd1)
  169- 179  F11.6 ---       SNRpxCCD2    ?=- Average signal-to-noise ratio
                                          (per pixel) of CCD2 (snr_px_ccd2)
  181- 192  F12.6 ---       SNRpxCCD3    ?=- Average signal-to-noise ratio
                                          (per pixel) of CCD3 (snr_px_ccd3)
  194- 205  F12.6 ---       SNRpxCCD4    ?=- Average signal-to-noise ratio
                                          (per pixel) of CCD4 (snr_px_ccd4)
  207- 218  E12.6 ---       chi2sp       ?=- Chi2 value of spectroscopic
                                          fitting (chi2_sp)
  220- 221  I2    %         pxUsedPerc   Percentage of spectrum used for
                                          spectroscopic fit (px_used_perc)
  223- 237  A15   ---       ModelName    Neural network model used for creating
                                          synthetic spectra (model_name)
  239- 253  A15   ---       ClosestModel Neural network model closest for
                                          Teff/logg/[Fe/H] combination
                                          (closest_model)
  255- 266  F12.6 s         compTime     ?=- Computation time spent on spectrum
                                          (comp_time)
  268- 269  I2    ---       RVcompNr     Nr RV cross-correlation fubxtion peaks
                                          (rv_comp_nr)
  271- 283  F13.6 km/s      RVcomp1      ?=- Radial velocity of primary source
                                          (rv_comp_1)
  285- 296  E12.6 km/s    e_RVcomp1      ?=- Uncertainty of RVcomp1
                                          (e_rv_comp_1)
  298- 305  F8.6  ---       RVcomp1p     ?=- Prominence of RVcomp1 in CCF
                                          (rv_comp_1_p)
  307- 318  F12.6 km/s      RVcomp2      ?=- Radial velocity of potential
                                          secondary source (rv_comp_2)
  320- 327  F8.6  km/s      RVcomp2h     ?=- Height of recomp2 in CCF
                                          (rv_comp_2_h)
  329- 336  F8.6  km/s      RVcomp2p     ?=- Prominence of RVcomp2 in CCF
                                          (rv_comp_2_p)
  338- 349  F12.6 km/s      RVGDR3       ?=- Radial velocity in Gaia DR3
                                          (rv_gaia_dr3)
  351- 359  F9.6  km/s    e_RVGDR3       ?=- Uncertainty of rv_gaia_dr3
                                          (e_rv_gaia_dr3)
  361- 370  F10.6 km/s      VbaryEff     ?=- Barycentric velocity applied to
                                          reduced spectra (v_bary_eff)
  372- 378  F7.2  K         Teff         ?=- Spectroscopic effective
                                          temperature (used for fitting) (teff)
  380- 391  E12.6 K       e_Teff         ?=- Uncertainty on teff (e_teff)
  393- 401  F9.6  [cm/s2]   logg         ?=- Surface gravity adjusted via
                                          parallax information (logg)
  403- 410  F8.6  [cm/s2] e_logg         ?=- Uncertainty on logg (e_logg)
  412- 421  F10.6 ---       [Fe/H]       ?=- Abundance of Fe and all other
                                          elements not fitted in GALAH
                                          (Fe: 1D-NLTE) (fe_h)
  423- 434  E12.6 ---     e_[Fe/H]       ?=- Uncertainty on [Fe/H] (e_fe_h)
       436  I1    ---     f_[Fe/H]       Quality flag on [Fe/H] (flag_fe_h) (G2)
  438- 448  F11.6 km/s      vmic         ?=- Microturbulence velocity (fitted)
                                          (vmic)
  450- 461  E12.6 km/s    e_vmic         ?=- Uncertainty on vmic (e_vmic)
  463- 473  F11.5 km/s      vsini        ?=- Broadening velocity
                                          (fitted sme.vsini with sme.vmac=0)
                                          (vsini)
  475- 486  E12.6 km/s    e_vsini        ?=- Uncertainty of vsini (e_vsini)
  488- 497  F10.6 ---       [Li/Fe]      ?=- Elemental abundance for [Li/Fe]
                                          (li_fe)
  499- 508  F10.6 ---     e_[Li/Fe]      ?=- Uncertainty on [Li/Fe] (e_li_fe)
  510- 511  I2    ---     f_[Li/Fe]      Quality bitmask flag of [Li/Fe]
                                          (flag_li_fe) (G2)
  513- 521  F9.6  ---       [C/Fe]       ?=- Elemental abundance for [C/Fe]
                                          (c_fe)
  523- 534  E12.6 ---     e_[C/Fe]       ?=- Uncertainty on [C/Fe] (e_c_fe)
  536- 537  I2    ---     f_[C/Fe]       Quality bitmask flag of [C/Fe]
                                          (flag_c_fe) (G2)
  539- 547  F9.6  ---       [N/Fe]       ?=- Elemental abundance for [N/Fe]
                                          (n_fe)
  549- 557  F9.6  ---     e_[N/Fe]       ?=- Uncertainty on [N/Fe] (e_n_fe)
  559- 560  I2    ---     f_[N/Fe]       Quality bitmask flag of [N/Fe]
                                          (flag_n_fe) (G2)
  562- 570  F9.6  ---       [O/Fe]       ?=- Elemental abundance for [O/Fe]
                                          (o_fe)
  572- 583  E12.6 ---     e_[O/Fe]       ?=- Uncertainty on [O/Fe] (e_o_fe)
  585- 586  I2    ---     f_[O/Fe]       Quality bitmask flag of [O/Fe]
                                          (flag_o_fe) (G2)
  588- 596  F9.6  ---       [Na/Fe]      ?=- Elemental abundance for [Na/Fe]
                                          (na_fe)
  598- 609  E12.6 ---     e_[Na/Fe]      ?=- Uncertainty on [Na/Fe] (e_na_fe)
  611- 612  I2    ---     f_[Na/Fe]      Quality bitmask flag of [Na/Fe]
                                          (flag_na_fe) (G2)
  614- 622  F9.6  ---       [Mg/Fe]      ?=- Elemental abundance for [Mg/Fe]
                                          (mg_fe)
  624- 635  E12.6 ---     e_[Mg/Fe]      ?=- Uncertainty on [Mg/Fe] (e_mg_fe)
  637- 638  I2    --      f_[Mg/Fe]      Quality bitmask flag of [Mg/Fe]
                                          (flag_mg_fe) (G2)
  640- 648  F9.6  ---       [Al/Fe]      ?=- Elemental abundance for [Al/Fe]
                                          (al_fe)
  650- 658  F9.6  ---     e_[Al/Fe]      ?=- Uncertainty on [Al/Fe] (e_al_fe)
  660- 661  I2    ---     f_[Al/Fe]      Quality bitmask flag of [Al/Fe]
                                          (flag_al_fe) (G2)
  663- 671  F9.6  ---       [Si/Fe]      ?=- Elemental abundance for [Si/Fe]
                                          (si_fe)
  673- 684  E12.6 ---     e_[Si/Fe]      ?=- Uncertainty on [Si/Fe] (e_si_fe)
  686- 687  I2    ---     f_[Si/Fe]      Quality bitmask flag of [Si/Fe]
                                          (flag_si_fe) (G2)
  689- 697  F9.6  ---       [K/Fe]       ?=- Elemental abundance for [K/Fe]
                                          (k_fe)
  699- 710  E12.6 ---     e_[K/Fe]       ?=- Uncertainty on [K/Fe] (e_k_fe)
  712- 713  I2    ---     f_[K/Fe]       Quality bitmask flag of [K/Fe]
                                          (flag_k_fe) (G2)
  715- 723  F9.6  ---       [Ca/Fe]      ?=- Elemental abundance for [Ca/Fe]
                                          (ca_fe)
  725- 736  E12.6 ---     e_[Ca/Fe]      ?=- Uncertainty on [Ca/Fe] (e_ca_fe)
  738- 739  I2    ---     f_[Ca/Fe]      Quality bitmask flag of [Ca/Fe]
                                          (flag_ca_fe) (G2)
  741- 749  F9.6  ---       [Sc/Fe]      ?=- Elemental abundance for [Sc/Fe]
                                          (sc_fe)
  751- 762  E12.6 ---     e_[Sc/Fe]      ?=- Uncertainty on [Sc/Fe] (e_sc_fe)
  764- 765  I2    ---     f_[Sc/Fe]      Quality bitmask flag of [Sc/Fe]
                                          (flag_sc_fe) (G2)
  767- 775  F9.6  ---       [Ti/Fe]      ?=- Elemental abundance for [Ti/Fe]
                                          (ti_fe)
  777- 788  E12.6 ---     e_[Ti/Fe]      ?=- Uncertainty on [Ti/Fe] (e_ti_fe)
  790- 791  I2    ---     f_[Ti/Fe]      Quality bitmask flag of [Ti/Fe]
                                          (flag_ti_fe) (G2)
  793- 801  F9.6  ---       [V/Fe]       ?=- Elemental abundance for [V/Fe]
                                          (v_fe)
  803- 811  F9.6  ---     e_[V/Fe]       ?=- Uncertainty on [V/Fe] (e_v_fe)
  813- 814  I2    ---     f_[V/Fe]       Quality bitmask flag of [V/Fe]
                                          (flag_v_fe) (G2)
  816- 824  F9.6  ---       [Cr/Fe]      ?=- Elemental abundance for [Cr/Fe]
                                          (cr_fe)
  826- 837  E12.6 ---     e_[Cr/Fe]      ?=- Uncertainty on [Cr/Fe] (e_cr_fe)
  839- 840  I2    ---     f_[Cr/Fe]      Quality bitmask flag of [Cr/Fe]
                                          (flag_cr_fe) (G2)
  842- 850  F9.6  ---       [Mn/Fe]      ?=- Elemental abundance for [Mn/Fe]
                                          (mn_fe)
  852- 863  E12.6 ---     e_[Mn/Fe]      ?=- Uncertainty on [Mn/Fe] (e_mn_fe)
  865- 866  I2    ---     f_[Mn/Fe]      Quality bitmask flag of [Mn/Fe]
                                          (flag_mn_fe) (G2)
  868- 876  F9.6  ---       [Co/Fe]      ?=- Elemental abundance for [Co/Fe]
                                          (co_fe)
  878- 886  F9.6  ---     e_[Co/Fe]      ?=- Uncertainty on [Co/Fe] (e_co_fe)
  888- 889  I2    ---     f_[Co/Fe]      Quality bitmask flag of [Co/Fe]
                                          (flag_co_fe) (G2)
  891- 899  F9.6  ---       [Ni/Fe]      ?=- Elemental abundance for [Ni/Fe]
                                          (ni_fe)
  901- 912  E12.6 ---     e_[Ni/Fe]      ?=- Uncertainty on [Ni/Fe] (e_ni_fe)
  914- 915  I2    ---     f_[Ni/Fe]      Quality bitmask flag of [Ni/Fe]
                                          (flag_ni_fe) (G2)
  917- 925  F9.6  ---       [Cu/Fe]      ?=- Elemental abundance for [Cu/Fe]
                                          (cu_fe)
  927- 938  E12.6 ---     e_[Cu/Fe]      ?=- Uncertainty on [Cu/Fe] (e_cu_fe)
  940- 941  I2    ---     f_[Cu/Fe]      Quality bitmask flag of [Cu/Fe]
                                          (flag_cu_fe) (G2)
  943- 951  F9.6  ---       [Zn/Fe]      ?=- Elemental abundance for [Zn/Fe]
                                          (zn_fe)
  953- 961  F9.6  ---     e_[Zn/Fe]      ?=- Uncertainty on [Zn/Fe] (e_zn_fe)
  963- 964  I2    ---     f_[Zn/Fe]      Quality bitmask flag of [Zn/Fe]
                                          (flag_zn_fe) (G2)
  966- 974  F9.6  ---       [Rb/Fe]      ?=- Elemental abundance for [Rb/Fe]
                                          (rb_fe)
  976- 984  F9.6  ---     e_[Rb/Fe]      ?=- Uncertainty on [Rb/Fe] (e_rb_fe)
  986- 987  I2    ---     f_[Rb/Fe]      Quality bitmask flag of [Rb/Fe]
                                          (flag_rb_fe) (G2)
  989- 997  F9.6  ---       [Sr/Fe]      ?=- Elemental abundance for [Sr/Fe]
                                          (sr_fe)
  999-1007  F9.6  ---     e_[Sr/Fe]      ?=- Uncertainty on [Sr/Fe] (e_sr_fe)
 1009-1010  I2    ---     f_[Sr/Fe]      Quality bitmask flag of [Sr/Fe]
                                          (flag_sr_fe) (G2)
 1012-1020  F9.6  ---       [Y/Fe]       ?=- Elemental abundance for [Y/Fe]
                                          (y_fe)
 1022-1033  E12.6 ---     e_[Y/Fe]       ?=- Uncertainty on [Y/Fe] (e_y_fe)
 1035-1036  I2    ---     f_[Y/Fe]       Quality bitmask flag of [Y/Fe]
                                          (flag_y_fe) (G2)
 1038-1046  F9.6  ---       [Zr/Fe]      ?=- Elemental abundance for [Zr/Fe]
                                          (zr_fe)
 1048-1057  F10.6 ---     e_[Zr/Fe]      ?=- Uncertainty on [Zr/Fe] (e_zr_fe)
 1059-1060  I2    ---     f_[Zr/Fe]      Quality bitmask flag of [Zr/Fe]
                                          (flag_zr_fe) (G2)
 1062-1070  F9.6  ---       [Mo/Fe]      ?=- Elemental abundance for [Mo/Fe]
                                          (mo_fe)
 1072-1080  F9.6  ---     e_[Mo/Fe]      ?=- Uncertainty on [Mo/Fe] (e_mo_fe)
 1082-1083  I2    ---     f_[Mo/Fe]      Quality bitmask flag of [Mo/Fe]
                                          (flag_mo_fe) (G2)
 1085-1093  F9.6  ---       [Ru/Fe]      ?=- Elemental abundance for [Ru/Fe]
                                          (ru_fe)
 1095-1103  F9.6  ---     e_[Ru/Fe]      ?=- Uncertainty on [Ru/Fe] (e_ru_fe)
 1105-1106  I2    ---     f_[Ru/Fe]      Quality bitmask flag of [Ru/Fe]
                                          (flag_ru_fe) (G2)
 1108-1116  F9.6  ---       [Ba/Fe]      ?=- Elemental abundance for [Ba/Fe]
                                          (ba_fe)
 1118-1129  E12.6 ---     e_[Ba/Fe]      ?=- Uncertainty on [Ba/Fe] (e_ba_fe)
 1131-1132  I2    ---     f_[Ba/Fe]      Quality bitmask flag of [Ba/Fe]
                                          (flag_ba_fe) (G2)
 1134-1142  F9.6  ---       [La/Fe]      ?=- Elemental abundance for [La/Fe]
                                          (la_fe)
 1144-1153  F10.6 ---     e_[La/Fe]      ?=- Uncertainty on [La/Fe] (e_la_fe)
 1155-1156  I2    ---     f_[La/Fe]      Quality bitmask flag of [La/Fe]
                                          (flag_la_fe) (G2)
 1158-1166  F9.6  ---       [Ce/Fe]      ?=- Elemental abundance for [Ce/Fe]
                                          (ce_fe)
 1168-1177  F10.6 ---     e_[Ce/Fe]      ?=- Uncertainty on [Ce/Fe] (e_ce_fe)
 1179-1180  I2    ---     f_[Ce/Fe]      Quality bitmask flag of [Ce/Fe]
                                          (flag_ce_fe) (G2)
 1182-1190  F9.6  ---       [Nd/Fe]      ?=- Elemental abundance for [Nd/Fe]
                                          (nd_fe)
 1192-1203  E12.6 ---     e_[Nd/Fe]      ?=- Uncertainty on [Nd/Fe] (e_nd_fe)
 1205-1206  I2    ---     f_[Nd/Fe]      Quality bitmask flag of [Nd/Fe]
                                          (flag_nd_fe) (G2)
 1208-1216  F9.6  ---       [Sm/Fe]      ?=- Elemental abundance for [Sm/Fe]
                                          (sm_fe)
 1218-1227  F10.6 ---     e_[Sm/Fe]      ?=- Uncertainty on [Sm/Fe] (e_sm_fe)
 1229-1230  I2    ---     f_[Sm/Fe]      Quality bitmask flag of [Sm/Fe]
                                          (flag_sm_fe) (G2)
 1232-1240  F9.6  ---       [Eu/Fe]      ?=- Elemental abundance for [Eu/Fe]
                                          (eu_fe)
 1242-1250  F9.6  ---     e_[Eu/Fe]      ?=- Uncertainty on [Eu/Fe] (e_eu_fe)
 1252-1253  I2    ---     f_[Eu/Fe]      Quality bitmask flag of [Eu/Fe]
                                          (flag_eu_fe) (G2)
 1255-1266  E12.6 km/s      sb2RV16      ?=- 16th perc. radial velocity of fit
                                          to syn-obs residuals (sb2_rv_16)
 1268-1279  E12.6 km/s      sb2RV50      ?=- 50th perc. radial velocity of fit
                                          to syn-obs residuals (sb2_rv_50)
 1281-1292  E12.6 km/s      sb2RV84      ?=- 84th perc. radial velocity of fit
                                          to syn-obs residuals (sb2_rv_84)
 1294-1304  F11.6 0.1nm     EWHbeta      ?=- Equivalent Width of observed Hbeta
                                          core (ew_h_beta)
 1306-1317  F12.6 0.1nm     EWHalpha     ?=- Equivalent Width of observed
                                          Halpha core (ew_h_alpha)
 1319-1329  F11.6 0.1nm     ResHbeta     ?=- Residual Equivalent Width of
                                          syn-obs residuals at Hbeta core
                                          (res_h_beta)
 1331-1342  F12.6 0.1nm     ResHalpha    ?=- Residual Equivalent Width of
                                          syn-obs residuals at Halpha core
                                          (res_h_alpha)
 1344-1353  F10.6 0.1nm     EWKis        ?=- Equivalent Width of fit for K7699
                                          Interstellar Line (ew_k_is)
 1355-1362  F8.6  ---       sigmaKis     ?=- Sigma of Gaussian fit for K7699
                                          Interstellar Line (sigma_k_is)
 1364-1374  F11.6 km/s      RVKis        ?=- Radial velocity of fit to syn-obs
                                          residuals around K7699 line (rv_k_is)
 1376-1384  F9.6  0.1nm     EWdib5780    ?=- Equivalent Width of fit for 5780
                                          Diffiuse Interstellar Band
                                          (ew_dib5780)
 1386-1393  F8.6  ---       sigdib5780   ?=- Sigma of Gaussian fit for 5780
                                          DIB (sigma_dib5780)
 1395-1406  F12.6 km/s      RVdib5780    ?=- Radial velocity of fit to syn-obs
                                          residuals around 5780 DIB (rv_dib5780)
 1408-1417  F10.6 0.1nm     EWdib5797    ?=- Equivalent Width of fit for 5797
                                          Diffiuse Interstellar Band
                                          (ew_dib5797)
 1419-1426  F8.6  ---       sigdib5797   ?=- Sigma of Gaussian fit for 5797 DIB
                                          (sigma_dib5797)
 1428-1439  F12.6 km/s      RVdib5797    ?=- Radial velocity of fit to syn-obs
                                          residuals around 5797 DIB (rv_dib5797)
 1441-1450  F10.6 0.1nm     EWdib6613    ?=- Equivalent Width of fit for
                                          6613 Diffiuse Interstellar Band
                                          (ew_dib6613)
 1452-1459  F8.6  ---       sigdib6613   ?=- Sigma of Gaussian fit for
                                          6613 DIB  (sigma_dib6613)
 1461-1472  F12.6 km/s      RVdib6613    ?=- Radial velocity of fit to syn-obs
                                          residuals around 6613 DIB
                                          (rv_dib6613)
 1474-1482  F9.6  [cm/s2]   loggspec     ?=- Spectrsocopic surface gravity
                                          estimate (logg_spec)
 1484-1495  E12.6 [cm/s2] e_loggspec     ?=- Uncertainty on loggspec
                                          (e_logg_spec)
 1497-1504  F8.5  mag       Gmag         ?=- Mean Gaia DR3 G-band apparent
                                          magnitude (phot_g_mean_mag)
 1506-1513  F8.5  mag       BPmag        ?=- Mean Gaia DR3 GBP-band apparent
                                          magnitude (phot_bp_mean_mag)
 1515-1523  F9.6  mag       BP-RP        ?=- Color of BP-RP bands (bp_rp)
 1525-1531  F7.3  mag       Hmag         ?=- 2MASS H-band magnitude (h_m)
 1533-1537  F5.3  mag     e_Hmag         ?=- Uncertainty on Hmag (h_msigcom)
 1539-1545  F7.3  mag       Ksmag        ?=- 2MASS Ks-band magnitude (ks_m)
 1547-1551  F5.3  mag     e_Ksmag        ?=- Uncertainty on Ksmag (ks_msigcom)
 1553-1560  F8.6  mag       E(B-V)       ?=- Extinction E(B-V) (ebv)
 1562-1571  F10.6 mag       W2mag        ?=- AllWISE W2-band magnitude (W2mag)
 1573-1580  F8.6  mag     e_W2mag        ?=- uncertainty of W2mag (e_W2mag)
 1582-1592  F11.6 mas       plx          ?=- Astrometric parallax used for
                                          GALAH DR4 (parallax)
 1594-1601  F8.6  mas     e_plx          ?=- Uncertainty on astrometric
                                          parallax used for GALAH DR4
                                          (e_parallax)
 1603-1613  F11.6 mas       plxGDR3      ?=- Gaia DR3 parallax
                                          (parallax_gaia_edr3)
 1615-1622  F8.6  mas     e_plxGDR3      ?=- uncertainty of plxGDR3
                                          (e_parallax_gaia_edr3)
 1624-1632  F9.6  ---       RUWE         ?=- RUWE reported by Gaia DR3 (ruwe)
 1634-1645  E12.6 pc        rmed         ?=- Median Distance used for
                                          calculating logg (plx) (r_med)
 1647-1658  E12.6 pc      b_rmed         ?=- Lower Limit Distance used for
                                          calculating logg (plx) (r_lo)
 1660-1671  E12.6 pc      B_rmed         ?=- Higher Limit Distance used for
                                          calculating logg (plx) (r_hi)
 1673-1680  F8.6  mag       AKs          ?=- Attenuation in Ks-band A(Ks) used
                                          for calculating logg (plx) (a_ks)
 1682-1689  F8.6  Msun      Mass         ?=- Mass used for calculating
                                          logg(plx) (mass)
 1691-1699  F9.6  Gyr       Age          ?=- Age estimated when calculating
                                          mass (age)
 1701-1708  F8.6  mag       BCKs         ?=- Bolometric Correction
                                          Ks, BC(Ks), used for calculating
                                          logg(plx) (bc_ks)
 1710-1721  E12.6 Lsun      Lbol         ?=- Bolometric Luminosity used for
                                          calculating logg (plx) (lbol)
--------------------------------------------------------------------------------

Global notes:
Note (G1): Major spectroscopic quality bitmask flag as follows:
      0 = No flag rasied
      1 = Emission in Halpha/Hbeta detected
      2 = Not all 4 CCDs available
      4 = Spectroscopic Binary SB1 (RV changes beyond 2sigma)
      8 = Spectroscopic Binary SB2 (Double line splitting)
     16 = Chi2 unusually low/high by 3 sigma
     32 = Broadening (vsini) warning
     64 = Microturbulence (vmic) warning
    128 = Double line splitting detected (SB2)
    256 = Temperature (teff) warning
    512 = Gravity (logg) warning
   1024 = [Fe/H] (fe_h) warning
   2048 = Signal-to-noise indicates spectra not reliable
   4096 = Not converged within 4 iterations
   8192 = Extrapolating spectrum model
  16384 = No spectroscopic analysis results available
Note (G2):
      0 = No flag rasied
      1 = Upper limit
      2 = No measurement available
      4 = No convergence
      8 = Measurement above limit: A(Li) >
           4; [X/Fe] > 1.5 for C,N,O,Y,Ba,La,Ce,Nd; [X/Fe] > 1.0 for the rest
     16 = Measurement below limit: A(Li) < 0; [X/Fe] < -0.5 for Mg, Si, Ti;
           [X/Fe] < -1.0 for the rest
     32 = Measurement issue of CNO
     64 = Fitting likely stuck at initial value


Acknowledgements:
    Copied at https://www.galah-survey.org/dr4/the_catalogues/

License: CC-BY-4.0

================================================================================
(End)                                      Patricia Vannier [CDS]    29-Jul-2026
