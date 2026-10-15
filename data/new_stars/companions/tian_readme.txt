J/ApJS/246/4    Catalog of ultrawide binary stars from Gaia DR2    (Tian+, 2020)
================================================================================
The separation distribution of ultrawide binaries across Galactic populations.
    Tian H.-J., El-Badry K., Rix H.-W., Gould A.
   <Astrophys. J. Suppl. Ser., 246, 4 (2020)>
   =2020ApJS..246....4T
================================================================================
ADC_Keywords: Stars, double and multiple; Positional data; Photometry; Optical;
              Radial velocities
Keywords: Binary stars ; Visual binary stars ; Astrometric binary stars ;
          Galaxy stellar content ; Galaxy stellar disks ; Galaxy stellar halos ;
          Star formation ; Astrostatistics

Abstract:
    We present an extensive and pure sample of ultrawide binary stars with
    separations of 0.01<~s/pc<~1 in the solar neighborhood. Using data
    from Gaia DR2, we define kinematic subpopulations via the systems'
    tangential velocities, i.e., disk-like (v_{perp},tot_<=40km/s),
    intermediate (v_{perp},tot_=40-85km/s), and halo-like
    (v_{perp},tot_>=85km/s) binaries, presuming that these velocity cuts
    represent a rough ordering in the binaries' age and metallicity.
    Through stringent cuts on astrometric precision, we can obtain pure
    binary samples at such wide separations with thousands of binaries in
    each sample. Fitting a smoothly broken power law for the separation
    distribution, we find that its slope at s=10^2.5-4^au is the same for
    all subpopulations, p(s){propto}s^{gamma}^ with {gamma}~-1.54.
    However, the logarithmic slope of p(s) steepens at s>~10^4^au. We find
    some evidences that the degree of steepening increases with the
    binaries' age, with a slope change of only {Delta}{gamma}~0.5 for
    disk-like stars, but {Delta}{gamma}>~1 for halo-like stars. This trend
    is contrary to what might be expected if steepening at wide
    separations were due to gravitational perturbations by molecular
    clouds or stars, which would preferentially disrupt disk binaries. If
    we were to interpret steepening at s>~10^4^au as a consequence of
    disruption by MAssive Compact Halo Objects (MACHOs), we would have to
    invoke a MACHO population inconsistent with other constraints. As a
    more plausible alternative, we propose a simple model to predict the
    separation distribution of wide binaries formed in dissolving star
    clusters. This model generically predicts {gamma}~-1.5 as observed,
    with steepening at larger separations due to the finite size of
    binaries' birth clusters.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records  Explanations
--------------------------------------------------------------------------------
ReadMe         80        .  This file
catalog.dat   638   807611  Catalog of all candidate wide binaries (a<1pc)
                             selected from Gaia DR2 in the solar neighborhood
                             with distances d<4kpc
--------------------------------------------------------------------------------

Description of files:
    Also available on the FTP, the file "All_candidates.fits" is the
    original catalog in FITS format.
    "WBs_Disklike.csv", "WBs_Intermediate.csv" and "WBs_Halolike.csv" are
    the threee subsets in original CSV format for disk-like
    (v_{perp},tot_<=40km/s), intermediate (v_{perp},tot_=40-85km/s), and
    halo-like (v_{perp},tot_>=85km/s) binaries, respectively.

See also:
 I/345           : Gaia DR2 (Gaia Collaboration, 2018)
 J/AJ/133/889    : Faint companions of Hipparcos stars (Lepine+, 2007)
 J/AJ/139/2566   : SLoWPoKES catalog (Dhital+, 2010)
 J/ApJ/790/158   : Candidate halo wide binaries (Allen+, 2014)
 J/MNRAS/472/675 : Wide binaries in Tycho-Gaia: search method (Andrews+, 2017)
 J/AJ/153/259    : The GAMBLES extension of the SLoWPoKES cat. (Oelkers+, 2017)
 J/AJ/155/149    : Co-moving stars observed by Gaia (Bochanski+, 2018)
 J/AJ/156/259    : Robo-AO detected close binaries in Gaia DR2 (Ziegler+, 2018)
 J/AJ/157/78     : Multiple star systems from GaiaDR2 (Jimenez-Esteban+, 2019)
 J/other/JDSO/15.434 : New Binary Systems from Data of Gaia DR2 (Romanov, 2019)

Byte-by-byte Description of file: catalog.dat
--------------------------------------------------------------------------------
   Bytes Format Units    Label     Explanations
--------------------------------------------------------------------------------
   1- 19  I19   ---      Gaia2     Gaia DR2 identifier of Star 2
                                    (785,036 uniques) (source_id2)
  21- 35 F15.11 deg      RA2deg    Gaia DR2 barycentric right ascension (ICRS)
                                    at Ep=2015.5 of Star 2 (ra2)
  37- 51 F15.11 deg      DE2deg    Gaia DR2 barycentric declination (ICRS)
                                    at Ep=2015.5 of Star 2 (dec2)
  53- 65  F13.6 ---      chi2AL2   [22/135246] Star 2 astrometric goodness of
                                    fit ({chi}^2^) in the along-scan direction
                                    (astrometric_chi2_al_2)
  67- 70  I4    ---      NgAL2     [44/1100] Star 2 number of good CCD transits
                                    (astrometric_n_good_obs_al2)
  72- 80  F9.3  ---      RFG2      [50/12061] Star 2 G-band mean flux divided
                                    by its error (phot_g_mean_flux_over_error2)
  82- 91  F10.3 ---      RFRP2     [10/192345] Star 2 integrated RP mean flux
                                    divided by its error
                                    (phot_rp_mean_flux_over_error2)
  93-101  F9.3  ---      RFBP2     [10/73309] Star 2 integrated BP mean flux
                                    divided by its error
                                    (phot_bp_mean_flux_over_error2)
 103-107  F5.3  ---      E(BR/RP)2 [1/2.5] Star 2 BP/RP excess factor
                                    (phot_bp_rp_excess_factor2)
 109-117  F9.3  mas/yr   pmRA2     [-2251/4106] Star 2 proper motion in right
                                    ascension direction (pmRA*cosDE) (pmra2)
 119-123  F5.3  mas/yr e_pmRA2     [0.01/4.5] pmRA2 standard error (pmra_error2)
 125-133  F9.3  mas/yr   pmDE2     [-3541/3156] Star 2 proper motion in
                                    declination (pmdec2)
 135-139  F5.3  mas/yr e_pmDE2     [0.01/4] pmDE2 standard error (pmdec_error2)
 141-147  F7.4  mag      Gmag2     [3.2/21] Star 2 G-band mean magnitude (Vega)
                                    (phot_g_mean_mag2)
 149-155  F7.4  mag      BPmag2    [3.6/24.8] Star 2 integrated BP mean
                                    magnitude (Vega) (phot_bp_mean_mag2)
 157-163  F7.4  mag      RPmag2    [2.7/21] Star 2 integrated RP mean
                                    magnitude (Vega) (phot_rp_mean_mag2)
 165-172  F8.4  mas      Plx2      [0.2/286.2] Star 2 absolute parallax
                                    (parallax2)
 174-182  F9.4  ---      RPlx2     [2/5416] Parallax divided by its error
                                    (parallax_over_error2)
 184-202  I19   ---      Gaia      Gaia DR2 identifier of Star 1 
                                    (782,120 uniques) (source_id)
 204-218 F15.11 deg      RAdeg     Gaia DR2 barycentric right ascension (ICRS)
                                    at Ep=2015.5 of Star 1 (ra)
 220-234 F15.11 deg      DEdeg     Gaia DR2 barycentric declination (ICRS)
                                    at Ep=2015.5 of Star 1 (dec)
 236-244  F9.3  mas/yr   pmRA      [-2241/4165] Star 1 proper motion in right
                                    ascension direction (pmRA*cosDE) (pmra)
 246-250  F5.3  mas/yr e_pmRA      [0.01/1.2] pmRA1 standard error (pmra_error)
 252-260  F9.3  mas/yr   pmDE      [-3543/3250] Star 1 proper motion in
                                    declination (pmdec)
 262-266  F5.3  mas/yr e_pmDE      [0.01/1.3] pmDE1 standard error (pmdec_error)
 268-274  F7.4  mag      Gmag      [3.2/20.1] Star 1 G-band mean magnitude
                                    (Vega) (phot_g_mean_mag)
 276-282  F7.4  mag      BPmag     [3.6/21.2] Star 1 integrated BP mean
                                    magnitude (Vega) (phot_bp_mean_mag)
 284-290  F7.4  mag      RPmag     [2.6/19.3] Star 1 integrated RP mean
                                    magnitude (Vega) (phot_rp_mean_mag)
 292-299  F8.4  mas      Plx       [0.2/286] Star 1 absolute parallax (parallax)
 301-309  F9.4  ---      RPlx      [20/6543] Parallax divided by its error
                                    (parallax_over_error)
 311-319  F9.2  ---      chi2AL    [15/189370] Star 1 astrometric goodness of
                                    fit ({chi}^2^) in the along-scan direction
                                    (astrometric_chi2_al)
 321-324  I4    ---      NgAL      [45/1089] Star 1 number of good CCD transits
                                    (astrometric_n_good_obs_al)
 326-334  F9.3  ---      RFG       [60/13735] Star 1 G-band mean flux divided
                                    by its error (phot_g_mean_flux_over_error)
 336-345  F10.3 ---      RFRP      [20/842507] Star 1 integrated RP mean flux
                                    divided by its error
                                    (phot_rp_mean_flux_over_error)
 347-357  F11.3 ---      RFBP      [20/5104061] Star 1 integrated BP mean flux
                                    divided by its error
                                    (phot_bp_mean_flux_over_error)
 359-363  F5.3  ---      E(BR/RP)  [1/2.1] Star 1 BP/RP excess factor
                                    (phot_bp_rp_excess_factor)
 365-373  F9.7  deg      Sepdeg    [0.0003/1.1] Angular separation between
                                    star 1 and star 2 (pairdistance)
 375-384  F10.3 au       Sepau     [10.9/201600] Projected physical separation
                                    between two stars (s_AU)
     386  I1    ---      Np        [1/9] Number of pairs in neighborhood
                                    for each binary (num)
     388  I1    ---      Type      [1/3] Type of wide binary (binary_type) (1)
 390-404 F15.11 deg      GLON      Star 1 galactic longitude (gl)
 406-420 F15.11 deg      GLAT      Star 1 galactic latitude (gb)
 422-436 F15.11 deg      GLON2     Star 2 galactic longitude (gl2)
 438-452 F15.11 deg      GLAT2     Star 2 galactic latitude (gb2)
 454-463  F10.5 arcsec   Sep       [1.2/3858] Angular separation between
                                    the two stars (theta)
 465-472  F8.6  mas/yr e_dmu       [0.02/1.5] dmu uncertainty (sig_dmu)
 474-484  F11.6 mas/yr   dmu       [0.0003/184] Total scalar difference in
                                    proper motions between the two stars
                                    ({Delta}{mu}) (see Section 2.1)
 486-496  F11.6 mas/yr   dmuOrb    [0.008/821] Maximum proper-motion difference
                                    permissible for a circular orbit of total
                                    mass 5M_{sun}_ (dmu_orbit)
 498-504  F7.4  mag      GMag2     [-5.2/16] Star 2 absolute G-band magnitude
                                    (MG2)
 506-512  F7.4  mag      GMag      [-3/16] Star 1 absolute G-band magnitude (MG)
 514-520  F7.2  km/s     RV2       [-325/725]? Star 2 spectroscopic barycentric
                                    radial velocity (radial_velocity2)
 522-526  F5.2  km/s   e_RV2       [0.1/20]? RV2 standard error
                                    (radial_velocity_error2)
 528-534  F7.2  km/s     RV        [-716/565]? Star 1 spectroscopic barycentric
                                    radial velocity (radial_velocity)
 536-540  F5.2  km/s   e_RV        [0.1/20]? RV1 standard error
                                    (radial_velocity_error)
 542-549  F8.3  km/s     UGaia     [-290/547]? Star 1 radial velocity in the
                                    Cartesian coordinates (positive to the
                                    Galactic center) (U_gaia)
 551-558  F8.3  km/s     VGaia     [-498/652]? Star 1 azimuthal velocity in the
                                    Cartesian coordinates (V_gaia)
 560-567  F8.3  km/s     WGaia     [-255/338]? Star 1 vertical velocity in the
                                    Cartesian coordinates (W_gaia)
 569-574  F6.3  km/s   e_UGaia     [0.03/25]? UGaia error (eU_gaia)
 576-581  F6.3  km/s   e_VGaia     [0.03/20]? VGaia error (eV_gaia)
 583-588  F6.3  km/s   e_WGaia     [0.02/20]? WGaia error (eW_gaia)
 590-597  F8.3  km/s     UGaia2    [-534/335]? Star 2 radial velocity in the
                                    Cartesian coordinates (positive to the
                                    Galactic center) (U_gaia2)
 599-606  F8.3  km/s     VGaia2    [-500/532]? Star 2 azimuthal velocity in the
                                    Cartesian coordinates (V_gaia2)
 608-615  F8.3  km/s     WGaia2    [-381/175]? Star 2 vertical velocity in the
                                    (W_gaia2)
 617-622  F6.3  km/s   e_UGaia2    [0.05/20]? UGaia2 error (eU_gaia2)
 624-629  F6.3  km/s   e_VGaia2    [0.03/20]? VGaia2 error (eV_gaia2)
 631-636  F6.3  km/s   e_WGaia2    [0.04/19]? WGaia2 error (eW_gaia2)
     638  A1    ---      Pop       [dhi ] Binary kinematic subpopulation;
                                    column added by CDS (2)
--------------------------------------------------------------------------------
Note (1): Type of wide binary as follows:
   1 = MS-WD ("main-sequence"-white dwarf; 42,839 occurrences),
   2 = WD-WD (440 occurrences),
   3 = MS-MS (764,332 occurrences).
Note (2): Kinematic subpopulations define via the systems' tangential
           velocities as follows:
   d = disk-like (v_{perp},tot_<=40km/s; 4362 occurrences),
   i = intermediate (v_{perp},tot_=40-85km/s; 10090 occurrences), and
   h = halo-like (v_{perp},tot_>=85km/s; 4351 occurrences).

          See section 3 for further explanations.
--------------------------------------------------------------------------------

History:
    Catalog (Table 2 description in the paper) in FITS and subsets in CSV
    sent by the author.
    27-May-2020:Insert into VizieR
    12-Nov-2020: Column "Pop" corrected

Acknowledgements:
    Haijun Tian [hjtian at lamost.org]

================================================================================
(End)                                      Emmanuelle Perret [CDS]   17-Apr-2020
