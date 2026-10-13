Column,Units,Description
ticID,,TIC catalog ID number 
gaiaDR2id,,Gaia DR2 source ID number
gaiaDR3id,,Gaia DR3 source ID number
ra,deg,Right ascension from Gaia DR2
dec,deg,Declination from Gaia DR2
p,,Radius ratio 
pErr,,Error on radius ratio 
t0Fit,BTJD,Central transit time from transit fit 
t0FitErr,BTJD,Error on central transit time from transit fit
periodBLS,days,Orbital period from BLS analysis 
duration,days,Transit duration 
durationErr,days,Error on transit duration 
denFit,g cm^-3,Stellar density from transit fit 
denFitErr,g cm^-3,Error on stellar density from transit fit 
u1,,Linear limb darkening coefficient 
u1Err,,Error on linear limb darkening coefficient 
u2,,Quadratic limb darkening coefficient
u2Err,,Error on quadratic limb darkening coefficient 
const,,Constant term in transit fit 
constErr,,Error on constant term in transit fit 
chir,,Reduced chi squared of the transit fit 
deltarho,,Difference of theoretical and measured stellar density normalized by its error
trdepth,ppm,Transit depth 
rmsoot,,Root mean square of out of transit measurements 
ar,,Ratio of semi-major axis to stellar radius 
arErr,,Error on ratio of semi-major axis to stellar radius 
b,,Impact parameter 
bErr,,Error on impact parameter 
i,deg, Inclination 
iErr,deg,Error on inclination 
t12t14,,Ratio of ingress time to total transit time
t12t14Err,,Error on ratio of ingress time to total transit time
rpj,R_Jupiter,Candidate radius 
rpjErr,R_Jupiter,Error on candidate radius  
rpe,R_Earth,Candidate radius 
rpeErr,R_Earth,Error on candidate radius 
rank1,,First Mahalanobis rank 
rank2,,Second Mahalanobis rank 
rank3,,Third Mahalanobis rank 
rank4,,Fourth Mahalanobis rank 
prf,per cent,Random Forest probability (as a percentage)
pd,per cent,Centroid motion probability (as a percentage) 
peta,per cent,Local average stellar number density probability (as a percentage)
t0Lin,BTJD,Central transit time from linear fit
t0LinErr,BTJD,Error on central transit time from linear fit
periodLin,days,Orbital period from linear fit
periodLinErr,days,Error on orbital period from linear fit
alpha,days,alpha parameter (see M20 paper)
beta,days,beta parameter (see M20 paper)
teff,K,Stellar effective temperature 
teffErr,K,Error on stellar effective temperature 
radius,R_solar,Stellar radius 
radiusErr,R_solar,Error on stellar radius
mass,M_solar,Stellar mass 
massErr,M_solar,Error on stellar mass 
rho,g cm^-3,Stellar density 
rhoErr,g cm^-3, Error on stellar density 
ContRatioTIC,,Contamination ratio from TIC 
tmag,,TESS magnitude from TIC
vmag,,V magnitude 
gmag,,G magnitude from Gaia DR2
rv,km s-1,Radial velocity from Gaia DR3
rvErr,km s^-1,Error on radial velocity from Gaia DR3
nrv,,Number of eligible transits used to derive the median radial velocity from Gaia DR3
rvChisqPvalue,,P-value for constanct based on a chi-squared criterion from Gaia DR3
rvRenormalisedGof,,Radial velocity renormalised goodness of fit from Gaia DR3
gofAl,,Goodness of fit in the along scan direction from Gaia DR3
astroExcess,mas,Excess of astrometric noise from Gaia DR3	   
astroExcessSig,,Significance of Excess of astrometric noise from Gaia DR3
ruwe,,Renormalised unit weight error from Gaia DR3
snI,,S/N ratio of primary transit
snII,,S/N ratio of secondary eclipse
snIII,,S/N ratio of tertiary eclipse
snOE,,Odd/Even S/N ratio
snOEFit,,Odd/Even S/N ratio from transit fit 
r2oot,,R^2_OOT parameter
q,,Fractional transit duration
sde,,Signal Detection Efficiency
sdeAL,,Signal Detection Efficiency of Alias peaks
p2pio,,In/Out of transit point-to-point noise
p2ps,,P2P_s parameter
r,R_Jupiter,Estimated candidate radius
ntr,,Number of transits
binflag,,Binarity flag (binflag=1, binary; binflag=0, not binary)
aperture,pix,Aperture photometry used to analyze the target
apstat,,Difference of transit depth on two adjacent apertures normalized by the square root of the quadratic sum of their errors (see M22)
btier,,A level flag based on the impact parameter value: btier=1; b<=0.2 | btier=2; 0.2<b<=0.4 | btier=3; 0.4<b<=0.6 | btier=4; 0.6<b<=0.8 | btier=5; b>0.8
refmask,,Bitmask flag identifying objects already present in other sources: bitmask=0 - New candidate | bitmask=1 - NASA exoplanet archive confirmed planet | bitmask=2 - TOI | bitmask=4 - cTOI | bitmask=8 - Olmschenk et al. 2021, AJ, 161, 273 | bitmask=16 - Eisner et al. 2021, MNRAS, 501, 4669 | bitmask=32 - Yu et al. 2019, AJ, 158, 25 
