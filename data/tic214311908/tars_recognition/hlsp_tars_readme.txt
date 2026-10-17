README for the TESS All-Sky Rotation Survey (TARS)
MAST webpage: https://archive.stsci.edu/hlsp/tars 
Refer to this HLSP with DOI: https://doi.org/10.17909/2sff-fn29 

Contributors: Andrew Boyle (UNC Chapel Hill), Luke Bouma (Carnegie
              Observatories), Andrew Mann (UNC Chapel Hill)
Date: 2026

================================================================================
General Information
================================================================================

The TESS All-Sky Rotation Survey (TARS) provides systematics-corrected
light curves and Lomb-Scargle periodograms for approximately 39 million
star-sector combinations drawn from 7.5 million unique target stars
observed by TESS between Sectors 1 and 96 (2018 July - 2025 September).
From these measurements, TARS produces a catalog of stellar rotation
periods for 1,046,317 stars with TESS magnitude T < 16 and distances
within 500 pc. The survey, data products, and catalog are described in
Boyle, Bouma, & Mann (2026).

================================================================================
Target Selection
================================================================================

Targets were selected from the TESS Input Catalog (TIC v8.2; Paegert
et al. 2021) satisfying T < 16 and Gaia DR2 parallax > 2 mas
(corresponding to d < 500 pc). After removing TIC duplicates, stars
lacking Gaia DR2 identifiers, and stars without Gaia DR3 G, BP, or RP
photometry, the team cross-matched with Gaia DR3 using the
dr2_dr3_neighbourhood table. Each remaining star was checked against
TESS pointing using tess-point (Burke et al. 2020) to determine sector
coverage. This procedure yielded 7,481,412 unique targets observed in at
least one TESS sector, corresponding to 39,061,674 star-sector
combinations.

================================================================================
Data Processing
================================================================================

For each star-sector combination, a light curve was extracted from the
TESS Calibrated Full Frame Images using the unpopular package (Hattori
et al. 2021). unpopular uses a causal pixel model that treats signals
common to many pixels as instrumental systematics and signals unique to
individual pixels as astrophysical. The team used a 50x50 pixel cutout, a
background estimated from the 300 faintest pixels, a systematics model
built from 100 pixels of similar brightness (excluding pixels within 5
of the target), a regularization parameter of 0.1, and 100 temporal
segments for cross-validation.

After systematics correction, a Lomb-Scargle periodogram was computed for each 
light curve using astropy. The periodogram was evaluated on a linearly spaced 
frequency grid from 0.1 to 20 days. The team recorded the periods, powers, and 
signal-to-noise ratios of the five highest peaks, and computed additional 
derived parameters including the best-fit sinusoid amplitude, point-to-point 
RMS, reduced chi-squared, and Bayesian Information Criterion values for one-
term and two-term sinusoidal models. Full definitions of all recorded 
parameters are given in Appendix A of Boyle, Bouma, & Mann (2026).

================================================================================
What is Being Released
================================================================================

This HLSP release comprises two types of data products for each of the
~39 million star-sector combinations:

  1. Light curve files (FITS): Systematics-corrected TESS light curves extracted with the unpopular package.

  2. Vetting plots (JPEG): Single-page visual summaries of each
     periodicity measurement.

The data are organized into directories by TESS sector.

================================================================================
File Naming Convention
================================================================================

Files follow the MAST HLSP naming convention:

  hlsp_tars_tess_ffi_s<SSSS>-<TTTTTTTTTTTTTTTT>_tess_<VV>_<product>.<ext>

where:

  tars        = HLSP project name
  tess        = mission
  ffi         = instrument configuration (full-frame images)
  <SSSS>      = TESS sector number, zero-padded to 4 digits
  <TTTTTTTTTTTTTTTT> = TIC ID, zero-padded to 16 digits
  tess        = filter (TESS bandpass)
  <VV>        = data release version (v01 for this release)
  <product>   = data product type (see below)
  <ext>       = file extension

Data product types:

  lc          = light curve (FITS)
  preview     = vetting plot (JPEG)

Example filenames:

  hlsp_tars_tess_ffi_s0083-0000000375546214_tess_v01_lc.fits
  hlsp_tars_tess_ffi_s0083-0000000375546214_tess_v01_preview.jpeg

================================================================================
FITS File Structure
================================================================================

Each light curve FITS file contains two HDUs:

  HDU 0 - PRIMARY (Header only)
  HDU 1 - LIGHTCURVE (Binary table)

--- HDU 0: PRIMARY ---

Contains target metadata in the header. No data array.

Keyword     Format    Description
--------    ------    -----------
TICID       int       TESS Input Catalog identifier
RA_TARG     float     Right ascension of target (deg, J2016.0, ICRS)
DEC_TARG    float     Declination of target (deg, J2016.0, ICRS)
TESSMAG     float     TESS magnitude (mag)
SECTOR      int       TESS sector number
CAMERA      int       TESS camera number (1-4)
CCD         int       TESS CCD number (1-4)
GAIADR2     int       Gaia DR2 source identifier
GAIADR3     int       Gaia DR3 source identifier
TSTART      float     Start time of observations (BTJD)
TSTOP       float     End time of observations (BTJD)
TELAPSE     float     Elapsed time of observations (days)
XPOSURE     float     Median cadence (days)
TIMESYS     string    Time system (TDB)
DATE-BEG    string    ISO start date of observations (UTC)
DATE-END    string    ISO end date of observations (UTC)
DATE-AVG    string    ISO midpoint date of observations (UTC)
HLSPNAME    string    HLSP name (TESS All-Sky Rotation Survey)
HLSPID      string    HLSP identifier (TARS)
HLSPLEAD    string    HLSP lead (Andrew Boyle)
HLSPVER     string    HLSP version (1.0)
HLSPTARG    string    HLSP target name (TIC <TICID>)
DOI         string    DOI for the HLSP
LICENSE     string    Data license (CC BY 4.0)
LICENURL    string    License URL

--- HDU 1: LIGHTCURVE ---

Systematics-corrected light curve. Contains two columns:

Column    Format    Unit    Description
------    ------    ----    -----------
TIME      float64   d       Barycentric TESS Julian Date (BTJD = BJD - 2457000)
FLUX      float32   --      Normalized flux (median ~ 1.0)

The LIGHTCURVE extension header also contains the following metadata
keywords:

  APERTURE  = Aperture size (TESS pixels)
  RA_TARG   = Right ascension of target (deg, J2016.0)
  DEC_TARG  = Declination of target (deg, J2016.0)
  DETECTOR  = Detector name (TESS Photometer)
  FILTER    = Filter name (TESS)
  RADESYS   = Coordinate reference frame (ICRS)

Periodogram statistics are stored as additional header keywords in this
extension:

Keyword     Format    Description
--------    ------    -----------
PER1        float     Period of highest periodogram peak (days)
POW1        float     Power of highest periodogram peak
PER2        float     Period of 2nd-highest peak (days)
POW2        float     Power of 2nd-highest peak
PER3        float     Period of 3rd-highest peak (days)
POW3        float     Power of 3rd-highest peak
PER4        float     Period of 4th-highest peak (days)
POW4        float     Power of 4th-highest peak
PER5        float     Period of 5th-highest peak (days)
POW5        float     Power of 5th-highest peak
SNR1        float     Peak 1 power / median power
SNR2        float     Peak 2 power / median power
SNR3        float     Peak 3 power / median power
SNR4        float     Peak 4 power / median power
SNR5        float     Peak 5 power / median power
AMPFIT      float     Amplitude of sine fit at PER1 (relative flux)
TMIN1       float     Time of first minimum in sine fit (BTJD)
PKWID1      float     Width of highest periodogram peak (days)
P2PRMS      float     Point-to-point noise estimate (relative flux)
A9010       float     90th - 10th percentile amplitude (relative flux)
REDCHI2     float     Reduced chi-squared of sine fit
SNR         float     SNR of light curve at PER1
BIC1T       float     BIC for 1-term sine model at PER1
BIC2T       float     BIC for 2-term sine model at 2x PER1
BIC2TS      float     BIC for 2-term sine model at PER1

Full definitions of these parameters are given in Appendix A of Boyle,
Bouma, & Mann (2026).

================================================================================
Vetting Plots
================================================================================

Each vetting plot is a single-page JPEG with four panels and a text header:


  Top:           Summary of stellar and observational properties,
                 including TIC and Gaia identifiers, sector/camera/CCD,
                 Gaia photometry, coordinates, proper motions, parallax,
                 effective temperature, RUWE, non-single-star flag,
                 classifier scores, and the top five periodogram
                 peaks.

  Top left:      Systematics-corrected light curve (flux vs. time in
                 BTJD). A horizontal bar marks the length of the detected period.

  Top right:     Lomb-Scargle periodogram with vertical dashed lines
                 at the detected period and twice the detected period.

  Bottom left:   Light curve phase-folded at the detected rotation
                 period.

  Bottom right:  Light curve phase-folded at twice the detected rotation
                 period.

In all light curve panels, gray points are individual measurements and black points
are phase-binned medians.

================================================================================
Related Data Products
================================================================================

The TARS rotation period catalog, a full table of sector-level
measurements, and a Python script for regenerating the catalog with user-defined
probability thresholds and quality cuts are available on Zenodo at:
  https://doi.org/10.5281/zenodo.19917941

================================================================================
References
================================================================================

The TARS paper is available from ADS at:
  https://ui.adsabs.harvard.edu/abs/2026arXiv260305586B/abstract

================================================================================
Contact
================================================================================

For questions about these data products, contact:
  Andrew Boyle (awboyle@unc.edu)

