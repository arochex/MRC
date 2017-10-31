# MRC

#Abstract:

Several studies use the emission lines fluxes of the star forming galaxies to derive the gas-phase metallicity by applying the so-called strong-line methods. These metallicity indicators are empirically calibrated using direct abundance measurements. Our main goal in this study is to present a new approach to derive oxygen abundances using all the possible combinations of strong line ratios sensitive to the metallicity, covering a wide range interval of this parameter,  through a large compilation of HII regions with Te-based estimations (Marino et al. 2013). Our method called MRC (Multiple Ratio Calibrator) provides with a calibration with a systematic error of 0.13 dex. We also present a comparison between this new method and the O3N2, N2 calibrations (Marino et al. 2013), the ONS method (Pilyugin et al. 2010) and the C method (Pilyugin et al. 2012). Furthermore, we perform an analysis of the abundance gradient of the galaxy NGC628 confirming previous works. 


# SUMMARY

pyMRC is a code that use multi emission line ratios to derive the oxygen abundance through a large sample of HII regions anchored to Te-based measurements. 

# MRC STEPS

Next, the elements and operation of the pyMRC code is presented. 


A. Derivation of the oxygen abundance from a flux table


A.1.  pyMRC.py : to derive the oxygen abundances values with a 2D polynomial of degree 5

#To run pyMRC.py yu need

< polyfitter.py : Is the script that calculate the interpolation and the 2D polynomial fitting and evaluation >

#How to run the pyMRC.py

pyMRC.py catalog_M13.csv coeff_5.fits

where: 

a) catalog_M13.csv:  is a catalog that comprises the fluxes and ratios of the stronger emission lines detected in each HII region (Marino et al. 2013). It uses the following format:

COLUMN1: ID, string, the ID of the galaxy.
COLUMN2: reference literature papers, string, the reference strings are: CAL --->CALIFA; BER12 -->Berg et al. (2012); BRE12 --->Bresolin et al. (2012); CrB09 --->Crowther & Bibby (2009); Cro09 --->Croxall et al. (2009); Est13 --->Esteban et al. (2013); Gbe10 --->Garcia-Benito et al. (2010); Gus12 -->Guseva et al. (2012); Had11 --->Hadfield & Crowther (2007); Keh11 --->Kehrig et al. (2011); Mon12 ---> Monreal-Ibero et al. (2012a); PMC09r# --->Perez-Montero & Contini (2009); P12--->Pilyugin et al. (2012) ; San12 --->Sanders et al. (2012); Sta13 --->Stasinska et al. (2013); Wes13 --->Westmoquette et al. (2013); ZaB11--->Zahid & Bresolin (2011); ZB12 ---> Zurita & Bresolin (2012). 
COLUMN3: [OII]3727, float, dust-extinction-corrected line flux. I adopt a Cardelli 1989 extiction law (1989ApJ...345..245C) with Rv=3.1 and AV/AHβ= 1.164. 
COLUMN4: [OII]3729, float, dust-extinction-corrected line flux. 
COLUMN5: [OIII]4363, float, dust-extinction-corrected line flux.
COLUMN6: [OIII]4959, float, dust-extinction-corrected line flux.
COLUMN7: [OIII]5007, float, dust-extinction-corrected line flux.
COLUMN8: [NII]5755, float, dust-extinction-corrected line flux.
COLUMN9: [NII]6548, float, dust-extinction-corrected line flux.
COLUMN10: [NII]6584, float, dust-extinction-corrected line flux.
COLUMN11: [SIII]6312, float, dust-extinction-corrected line flux.
COLUMN12: [SII]6717, float, dust-extinction-corrected line flux. 
COLUMN13: [SII]6731, float, dust-extinction-corrected line flux. 
COLUMN14: jD, integer, the index jD indicates what doublets ([OII], [SII]) were resolved in the original source. (0=[OII], [SII] are resolved; 1= only [SII] is resolved and we assume [OII]3729/[OII]3726=1.5; 2= none is resolved, we assume [OII]3729/[OII]3726=1.5 and [SII]6717/[SII]6731=1.5)
COLUMN15: jT, integer, the index jT is 1 when the electron temperature is derived from the auroral and nebular lines of [O III], is 2 when the [N II] lines are used, and 3 for the case of [S III].
COLUMN16: T3 (in units of 10^4 K).
COLUMN17: 12+log(O/H), float, Abundance computed from the electron temperature.
 

b) coeff_5.fits: is the file.fits that contains the coefficients of the 2D polynomial to evaluate.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

A2. pyMRCeOH.py : to derive the oxygen abundances plus the errors in the OH abundances


#How to run the pyMRCeOH.py

pyMRCeOH.py catalog_M13+err.csv 1 500 coeff_5.fits




where: 

a) catalog_M13.csv:  is a catalog that comprises the fluxes  plus flux errors  in each HII region (Marino et al. 2013). It uses the following format:

COLUMNA1:  [OII]3727+[OII]3729
COLUMNA2:  e_([OII]3727+[OII]3729), error del flujo de linea
COLUMNA3:  Hbeta_4861, float, 10e-16 erg s^-1 cm^-2, flujo observado de la linea Hb
COLUMNA4:  e_Hbeta_4861, error del flujo observado de la linea Hb
COLUMNA5: [OIII]5007
COLUMNA6: e_[OIII]5007
COLUMNA7: Ha_6563
COLUMNA8: e_Ha_6563
COLUMNA9: [NII]6584
COLUMNA10: e_[NII]6584
COLUMNA11: [SII]6717
COLUMNA12: e_[SII]6717
COLUMNA13: [SII]6731
COLUMNA14: e_[SII]6731
COLUMNA15: OH_Te, Marino et al. (2013)


b) Gaussian Width : 1  
c) MC iterations :  100, 500, 1000, ...
d) polynomial coefficiente of degree 5: coeff_5.fits


The output is:

#print ' OH_out '+str(OH_out_final) +'  : e_OH_out '+str(e_OH_out_final)

or 

#ascii.write([OH_out_final,e_OH_out_final],'output.csv', names=['OH_mean_out','e_OH_mean_out'],format= 'csv',fast_writer=False,formats={'OH_mean_out':'%1.2f','e_OH_mean_out':'%0.2f'})


------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

B) Derivation of the coefficients of the polynomial  (for updated HII catalogs)


B1. poly_coeff_MC.py : to calculate the polynomial coeffcients


#How to run the code:

poly_coeff_MC.py  catalog_M13.csv 0.15 5 500


Command line entries: 

a. Te-based sample : catalog_M13.csv 
b. estimated error: 0.15
c. polynomial degree: 1,2,3,4,5,...,12
d. Monte Carlo iterations: 500


#Output:  

		outputfile='coeffs_'+str(n)+'.fits'
		hdulist.writeto(outputfile,clobber=True)



------------------------------------------------------------------------------------------------------------------------------

B2. ratios_nn_interp_comparisons.py : to calculate the interpolated 2D maps for the multiple emission line ratios comparisons



plots_ratios_interp2D.py catalog_M13.csv 5


Command line entries:

a. Te-based sample: catalog_M13.csv
b. polynomial coefficient of degree 5 

-----------------------------------------------------------------------------------------------------------------------------------------------


C) Example of its use on a particular galaxy (in this work: NGC628)



C1.  grad_rad_NGC628.py : to derive the oxygen abundances for the NGC628 galaxy and its radial gradients.


#How tu run the code:


grad_rad_NGC628.py catalog_flux_ngc628.csv catalog_coords_ngc628.csv OH_628_cor.csv 5.fits


Entries:
a. flux_catalog_NGC628
b. coordinates_catalog_NGC628
c. OH_abundaces corrected by extintion
d. polynomial coefficient

----------------------------------------------------------------------------------------------------------------------


# REFERENCES

