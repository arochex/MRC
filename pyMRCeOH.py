
#!/usr/bin/env python
#-*- coding-utf-8 -*-


__autores__ = "A.A.Aroche & S.F. Sanchez"
__date__ = "2017"

#Abundances 

import sys
import csv
import argparse
import matplotlib.pyplot as plt
from matplotlib import mpl
from scipy.stats import gaussian_kde
import pyfits
import itertools
from numpy import *
import numpy as np
import math
import pandas as pd
from pylab import polyfit,poly1d
from scipy import ndimage
from matplotlib.mlab import griddata
from astropy.table import Table
import random as rnd
import numpy.ma as ma
import matplotlib.cm as cm
#import seaborn as sns
#sns.set(style="white", color_codes=True)
#sns.set()
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
from astropy.io import ascii
import polyfitter as plf

#======================DEFINITIONS======================#
#SLEs
'''
R2 = OII3727
R3 = OIII5007 
N2= NII6584 
S2a = SII6717  
S2b = SII6731  
'''
#=======================CONSTANTS=======================#

R_25 = 14.1  #Optical radius in (kpc) at the B_25 mag arcsec(-2). Rosales-Ortega et al. (2010)
rho = 314   #arcsec #Size of the optical radius $\rho_{25}$ in arsec
A = 58   #arcsec
B = 5    #adimensionless
C = 44   #1/rads


#--------------------------------------------------------------------------#

#-------Command line input files--------------#

usage = 'Usage: %s infile' % sys.argv[0]

try:
    input_file = sys.argv[1]      #catalogue_Marino
    input_file = sys.argv[2]      #guess_error
    input_file = sys.argv[3]      #MC_iterations 
    input_file = sys.argv[4]      #coeffs_out
except:
    print usage; sys.exit(1)

guess_error = float(sys.argv[2]) 
n_mc = int(sys.argv[3]) 

#--------------------Reading the fits file--------------------#

input_fits = sys.argv[4]
hdu_coeffs = pyfits.open(input_fits)

coeffs_out = hdu_coeffs[0].data
e_coeffs_out = hdu_coeffs[1].data/np.sqrt(500)   #standard error

#e_coeffs_out = e_coeffs_out / np.sqrt(500)

#f, ax = plt.subplots()
#ax.imshow(coeffs_out)
#plt.show()

e_c_R3N2_R23  =e_coeffs_out[0,:]
e_c_R3N2_N2   =e_coeffs_out[1,:]
e_c_R23_N2    =e_coeffs_out[2,:]
e_c_R23_N2R2  =e_coeffs_out[3,:]
e_c_N2R2_N2   =e_coeffs_out[4,:]
e_c_N2_R3     =e_coeffs_out[5,:]
e_c_R3N2_R3R2 =e_coeffs_out[6,:]
e_c_R3R2_N2R2 =e_coeffs_out[7,:]
e_c_R3R2_N2   =e_coeffs_out[8,:]
e_c_N2S2_R3N2 =e_coeffs_out[9,:]
e_c_N2S2_R23  =e_coeffs_out[10,:]
e_c_N2S2_N2R2 =e_coeffs_out[11,:]
e_c_N2S2_N2   =e_coeffs_out[12,:]


c_R3N2_R23  =coeffs_out[0,:]
c_R3N2_N2   =coeffs_out[1,:]
c_R23_N2    =coeffs_out[2,:]
c_R23_N2R2  =coeffs_out[3,:]
c_N2R2_N2   =coeffs_out[4,:]
c_N2_R3     =coeffs_out[5,:]
c_R3N2_R3R2 =coeffs_out[6,:]
c_R3R2_N2R2 =coeffs_out[7,:]
c_R3R2_N2   =coeffs_out[8,:]
c_N2S2_R3N2 =coeffs_out[9,:]
c_N2S2_R23  =coeffs_out[10,:]
c_N2S2_N2R2 =coeffs_out[11,:]
c_N2S2_N2   =coeffs_out[12,:]


#---------------------HII regions catalogues----------------------#

#-----------Marino et al. (2013)---------------#
#Clipped catalogue
input_file = sys.argv[1]

g = pd.read_csv(input_file, comment = '#', header=None)


e_R2_in  = ma.masked_where(g[1] == 0,g[1])
e_R3_in  = ma.masked_where(g[5] == 0,g[5])
e_N2_in  = ma.masked_where(g[9] == 0,g[9])
e_S2a_in = ma.masked_where(g[11] == 0,g[11]) 
e_S2b_in = ma.masked_where(g[13] == 0,g[13])
e_S2_in  = ma.masked_where((e_S2a_in + e_S2b_in) == 0,(e_S2a_in + e_S2b_in))


R2_in  = ma.masked_where(g[0] == 0,g[0])
R3_in  = ma.masked_where(g[4] == 0,g[4])
N2_in  = ma.masked_where(g[8] == 0,g[8])
S2a_in = ma.masked_where(g[10] == 0,g[10])
S2b_in = ma.masked_where(g[12] == 0,g[12]) 
S2_in  = ma.masked_where((S2a_in+S2b_in) == 0,(S2a_in+S2b_in))



#print 'nHII ='+str(nHII)
#sys.exit(0)

nHII=len(R2_in)
matrix_OH=np.zeros((nHII,n_mc))


for repeticiones in range(0,11):
    for nc in range(0,n_mc):
    

        R2 = R2_in  + rnd.gauss(0,guess_error)*e_R2_in
        R3  = R3_in  + rnd.gauss(0,guess_error)*e_R3_in
        N2   = N2_in  + rnd.gauss(0,guess_error)*e_N2_in
        #S2a   = S2a_in + rnd.gauss(0,guess_error)*e_S2a_in
        #S2b   = S2b_in + rnd.gauss(0,guess_error)*e_S2b_in
        #S2 = (S2a_in+S2b_in) +  rnd.gauss(0,guess_error) * (e_S2a_in+e_S2b_in)
        #S2 = S2a_in +  rnd.gauss(0,guess_error) * e_S2a_in
        S2a = S2a_in
        S2b = S2b_in
        S2 = (S2a_in+S2b_in) +  rnd.gauss(0,guess_error) * (e_S2a_in+e_S2b_in)

        R3N2 = R3/N2
        R23  = R2+R3
        R3R2 = R3/R2
        N2R2 = N2/R2
        N2S2 = N2/S2

        c_R3N2_R23  =coeffs_out[0,:]  + rnd.gauss(0,guess_error) * e_c_R3N2_R23
        c_R3N2_N2   =coeffs_out[1,:]  + rnd.gauss(0,guess_error) * e_c_R3N2_N2
        c_R23_N2    =coeffs_out[2,:]  + rnd.gauss(0,guess_error) * e_c_R23_N2 
        c_R23_N2R2  =coeffs_out[3,:]  + rnd.gauss(0,guess_error) * e_c_R23_N2R2
        c_N2R2_N2   =coeffs_out[4,:]  + rnd.gauss(0,guess_error) * e_c_N2R2_N2 
        c_N2_R3     =coeffs_out[5,:]  + rnd.gauss(0,guess_error) * e_c_N2_R3
        c_R3N2_R3R2 =coeffs_out[6,:]  + rnd.gauss(0,guess_error) * e_c_R3N2_R3R2
        c_R3R2_N2R2 =coeffs_out[7,:]  + rnd.gauss(0,guess_error) * e_c_R3R2_N2R2
        c_R3R2_N2   =coeffs_out[8,:]  + rnd.gauss(0,guess_error) * e_c_R3R2_N2
        c_N2S2_R3N2 =coeffs_out[9,:]  + rnd.gauss(0,guess_error) * e_c_N2S2_R3N2
        c_N2S2_R23  =coeffs_out[10,:] + rnd.gauss(0,guess_error) * e_c_N2S2_R23
        c_N2S2_N2R2 =coeffs_out[11,:] + rnd.gauss(0,guess_error) * e_c_N2S2_N2R2
        c_N2S2_N2   =coeffs_out[12,:] + rnd.gauss(0,guess_error) * e_c_N2S2_N2

        #    print 'c_R3N2_R2 = '+str(c_R3N2_R23)
        #    print 'e_c_R3N2_R2 = '+str(e_c_R3N2_R23)
    
        OH_R3N2_R23  = plf.polyval2d(log10(R3N2),log10(R23),c_R3N2_R23)
        OH_R3N2_N2   = plf.polyval2d(log10(R3N2),log10(N2),c_R3N2_N2)
        OH_R23_N2    = plf.polyval2d(log10(R23),log10(N2),c_R23_N2)
        OH_R23_N2R2  = plf.polyval2d(log10(R23),log10(N2R2),c_R23_N2R2)
        OH_N2R2_N2   = plf.polyval2d(log10(N2R2),log10(N2),c_N2R2_N2)
        OH_N2_R3     = plf.polyval2d(log10(N2),log10(R3),c_N2_R3)
        OH_R3N2_R3R2 = plf.polyval2d(log10(R3N2),log10(R3R2),c_R3N2_R3R2)
        OH_R3R2_N2R2 = plf.polyval2d(log10(R3R2),log10(N2R2),c_R3R2_N2R2)
        OH_R3R2_N2   = plf.polyval2d(log10(R3R2),log10(N2),c_R3R2_N2)
        OH_N2S2_R3N2 = plf.polyval2d(log10(N2S2),log10(R3N2),c_N2S2_R3N2)
        OH_N2S2_R23  = plf.polyval2d(log10(N2S2),log10(R23),c_N2S2_R23)
        OH_N2S2_N2R2 = plf.polyval2d(log10(N2S2),log10(N2R2),c_N2S2_N2R2)
        OH_N2S2_N2   = plf.polyval2d(log10(N2S2),log10(N2),c_N2S2_N2)
        
        OH_R3N2_R23  =  ma.masked_invalid(OH_R3N2_R23)
        OH_R3N2_N2   =  ma.masked_invalid(OH_R3N2_N2)
        OH_R23_N2    =  ma.masked_invalid(OH_R23_N2)
        OH_R23_N2R2  =  ma.masked_invalid(OH_R23_N2R2)
        OH_N2R2_N2   =  ma.masked_invalid(OH_N2R2_N2)
        OH_N2_R3     =  ma.masked_invalid(OH_N2_R3)
        OH_R3N2_R3R2 =  ma.masked_invalid(OH_R3N2_R3R2)
        OH_R3R2_N2R2 =  ma.masked_invalid(OH_R3R2_N2R2)
        OH_R3R2_N2   =  ma.masked_invalid(OH_R3R2_N2)
        OH_N2S2_R3N2 =  ma.masked_invalid(OH_N2S2_R3N2)
        OH_N2S2_R23  =  ma.masked_invalid(OH_N2S2_R23)
        OH_N2S2_N2R2 =  ma.masked_invalid(OH_N2S2_N2R2)
        OH_N2S2_N2   =  ma.masked_invalid(OH_N2S2_N2)


        #print OH_N2S2_N2R2


        OH_array = np.array([OH_R3N2_R23,OH_R3N2_N2,OH_R23_N2,OH_R23_N2R2,OH_N2R2_N2,OH_N2_R3,OH_R3N2_R3R2,OH_R3R2_N2R2,OH_R3R2_N2,OH_N2S2_R3N2,OH_N2S2_R23,OH_N2S2_N2R2,OH_N2S2_N2])

    
        OH_mean_out = ma.masked_where(np.median(OH_array,axis=0) == 0,np.median(OH_array,axis=0))
        matrix_OH[:,nc]=OH_mean_out

        matrix_OH_masked = ma.masked_invalid(matrix_OH)
        OH_out_final=np.median(matrix_OH_masked,axis=1)
        e_OH_out_final=np.std(matrix_OH_masked,axis=1)

        #print ' OH_out '+str(OH_out_final) +'  : e_OH_out '+str(e_OH_out_final)

#ascii.write([OH_out_final,e_OH_out_final],'output.csv', names=['OH_mean_out','e_OH_mean_out'],format= 'csv',fast_writer=False,formats={'OH_mean_out':'%1.5f','e_OH_mean_out':'%0.5f'})


