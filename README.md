# MRC

#Abstract:

Several studies use the emission lines fluxes of the star forming galaxies to derive the gas-phase metallicity by applying the so-called strong-line methods. These metallicity indicators are empirically calibrated using direct abundance measurements. Our main goal in this study is to present a new approach to derive oxygen abundances using all the possible combinations of strong line ratios sensitive to the metallicity, covering a wide range interval of this parameter,  through a large compilation of HII regions with Te-based estimations (Marino et al. 2013). Our method called MRC (Multiple Ratio Calibrator) provides with a calibration with a systematic error of 0.13 dex. We also present a comparison between this new method and the O3N2, N2 calibrations (Marino et al. 2013), the ONS method (Pilyugin et al. 2010) and the C method (Pilyugin et al. 2012). Furthermore, we perform an analysis of the abundance gradient of the galaxy NGC628 confirming previous works. 


# SUMMARY

pyMRC is a code that use multi emission line ratios to derive the oxygen abundance through a large sample of HII regions anchored to Te-based measurements. 

# MRC DIAGRAM 

Next, the elements and operation of the pyMRC code is presented. 


O/H [MRC] = F(Pi,Pj) 

1. Te-based HII regions sample
2. O/H (Te-based)  
3. Multiple emission line ratios

4. Natural Neighbour Interpolation 
5. 2D polynomial fitting  	
6. 2D polynomial coefficients
7. O/H (MRC) 

# HOW TO USE THE CODES


A. polyfitter.py : Is the script that calculate the interpolation and the 2D polynomail fitting and evaluation.
B. coeff_#.fits : Are the 2D polynomial coefficients necessary to do the polynomial evaluation


1. poly_coeff_MC.py : to calculate the polynomial coeffcients  


poly_coeff_MC.py  catalog_M13.csv 0.15 5 500


Command line entries: 

a. Te-based sample : catalog_M13.csv 
b. estimated error: 0.15
c. polynomial degree: 1,2,3,4,5,...,12
d. Monte Carlo iterations: 500

#
2. plots_ratios_interp2D.py : to calculate the interpolated 2D maps for the multiple emission line ratios comparisons


plots_ratios_interp2D.py catalog_M13.csv 5


Command line entries:

a. Te-based sample: catalog_M13.csv
b. polynomial coefficient of degree 5 



3. pyMRC.py : to derive the oxygen abundances values with a 2D polynomial of degree 5


pyMRC.py catalog_M13.csv coeff_5.fits


Entries:

a. Te-based sample M13
b. polynomial coefficient of degree n




4. pyMRCeOH.py : to derive the oxygen abundances plus the errors in the OH 

pyMRCeOH.py catalog_M13+err.csv 1 500 coeff_5.fits

Entries:

a. Te-based sample plus errors flux lines 
b. Gaussian Width = 1  
c. MC iterations
d. polynomial coefficiente of degree 5



5. grad_rad_NGC628.py : to derive the oxygen abundances for the NGC628 galaxy and its radial gradients.



grad_rad_NGC628.py catalog_flux_ngc628.csv catalog_coords_ngc628.csv OH_628_cor.csv 5.fits


Entries:

a. flux_catalog_NGC628
b. coordinates_catalog_NGC628
c. OH_abundaces corrected by extintion
d. polynomial coefficient


# REFERENCES

