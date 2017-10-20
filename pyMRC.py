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
from mpl_toolkits.mplot3d import Axes3D
import numpy.ma as ma
import matplotlib.cm as cm
#import seaborn as sns
#sns.set(style="white", color_codes=True)
#sns.set()
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
import polyfitter as plf


#============Command line input files============#

usage = 'Usage: %s infile' % sys.argv[0]

try:
    input_file = sys.argv[1]      #catalogo de regiones HII. M13
    input_file = sys.argv[2]      #coeficientes de orden n
except:
    print usage; sys.exit(1)


#--------------------Reading the fits file--------------------#

input_fits = sys.argv[2]
hdu_coeffs = pyfits.open(input_fits)

coeffs_out = hdu_coeffs[0].data
e_coeffs_out = hdu_coeffs[1].data
#print np.array(e_coeffs_out)

#f, ax = plt.subplots()
#ax.imshow(coeffs_out)
#plt.show()


#---------------------HII regions catalogue----------------------#

#-----------P12.Marino et al. (2013)---------------#


input_file = sys.argv[1]

df = pd.read_csv(input_file, comment = '#', header=None)


OII3727  = df[2]
OII3729  = df[3]
OIII4959 = df[5]
OIII5007 = df[6]
NII6548  = df[8]
NII6584  = df[9]
SII6717  = df[11]
SII6731  = df[12]
OH_Te = df[16]

'''
R2 = OII3727 + OII3729
R3 = OIII4959 + OIII5007
N2 = NII6548 + NII6584
S2 = SII6717 + SII6731
'''
R2 = OII3727 + OII3729
R3 = 1.33*OIII5007
N2 = 1.33*NII6584
S2 = SII6717 + SII6731

#=====Pares de cocientes propuestos=====#

N2R3 = N2/R3 
R3N2 = R3/N2
R23  = R2+R3
R3R2 = R3/R2
N2R2 = N2/R2
N2S2 = N2/S2

#========Excitation parameter=========#

P = R3/(R23)

#=========================CALIBRATORS======================#

#--------O3N2 & N2 methods. Marino et al. (2013)-----------#


M_O3N2 = 8.533 - 0.214*log10((2.87)*(R3_M/N2_M))

M_N2 = 8.743 + 0.462*log10((1/2.86)*N2_M)



#=========ONS method. Pilyugin et al. (2010)============#

OH_ONS = lambda a, R2, R3, N2, S2, P: (a[0] + a[1]*P + a[2]*np.log10(R3) + a[3]*np.log10(N2/R2)+ a[4]*np.log10(S2/R2))

a_cool = (8.277, 0.657, -0.399, -0.061, 0.005)
a_warm = (8.816, -0.733, 0.454, 0.710, -0.337)
a_hot = (8.774, -1.855, 1.517, 0.304, 0.328)

Pil = []

for j in range(0,len(df[6])):
    
    if log10(N2[j]) > -0.1:
        ONS = OH_ONS(a_cool, R2[j], R3[j], N2[j], S2[j], P[j])

    if log10(N2[j]) < -0.1 and log10(N2[j]/S2[j]) > -0.25:
        ONS = OH_ONS(a_warm, R2[j], R3[j], N2[j], S2[j], P[j])
    
    if log10(N2[j]) < -0.1 and log10(N2[j]/S2[j]) < -0.25:
        ONS = OH_ONS(a_hot, R2[j], R3[j], N2[j], S2[j], P[j])

    Pil.append(ONS)


#=========================MRC_Calibrator========================#

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

                
            
OH_vec = np.array([OH_R3N2_R23,OH_R3N2_N2,OH_R23_N2,OH_R23_N2R2,OH_N2R2_N2,OH_N2_R3,OH_R3N2_R3R2,OH_R3R2_N2R2,OH_R3R2_N2,OH_N2S2_R3N2,OH_N2S2_R23,OH_N2S2_N2R2,OH_N2S2_N2])

OH_MRC = np.median(OH_vec,axis=0)

Delta_OH =  OH_Te - OH_MRC
mean_Delta_OH = np.mean(Delta_OH)
sigma_Delta_OH = np.std(Delta_OH)

#print ' OH_out '+str(OH_MRC)
print 'offset O/H =', mean_Delta_OH, sigma_Delta_OH

#ascii.write([OH_MRC],'OH_MRC_P12.csv', names=['OH_MRC'],format= 'csv',fast_writer=False,formats={'OH_MRC':'%1.2f'})


def plot_OH(param1, param2, tag1, tag2):

    size=60
    alpha = 0.6
    cm = plt.cm.get_cmap('nipy_spectral')
    
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor('white')

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="7%", pad=0.01)
    #ax.axvline(x=7.8, c='r',linestyle='-',linewidth=0.5)
    #ax.axhline(y=8.7684, c='b',linestyle='-',linewidth=0.5)
    #ax.axhline(y=8.1692, c='b',linestyle='-',linewidth=0.5)
    #ax.axvline(x=8.1692, c='r',linestyle='--',linewidth=0.2)
    #ax.axvline(x=8.7684, c='r',linestyle='--',linewidth=0.2)
    #ax.axvline(x=8.0038, c='r',linestyle='--',linewidth=0.2)
    #ax.axvline(x=8.6506, c='r',linestyle='--',linewidth=0.2)
    #ax.axvline(x=2.0, c='g',linestyle='--',linewidth=1.5)
    #ax.axhline(y=7.8, c='k',linestyle='--',linewidth=1.5)
    #ax.axvline(x=-2.5, c='k',linestyle='--',linewidth=1.5)
    #ax.axvline(x=-0.3, c='k',linestyle='--',linewidth=1.5)
    #ax.axvline(x=-1.6, c='r',linestyle='--',linewidth=1.5)
    #ax.axvline(x=-0.2, c='r',linestyle='--',linewidth=1.5)
    ax.set_ylim((6.5,9.5))
    ax.set_xlim((6.5,9.5))
    #ax.set_ylim((6.5,9.5))
    ax.plot((6.5,9.5), (6.5,9.5), 'k',linewidth=1)
    ax.plot((6.72,9.72), (6.5,9.5), 'k--',linewidth=1.5)
    ax.plot((6.5,9.5), (6.72,9.72), 'k--',linewidth=1.5)
    #ax.plot((6.64,9.64), (6.5,9.5), 'r--',linewidth=1.5)
    #ax.plot((6.5,9.5), (6.64,9.64), 'r--',linewidth=1.5)

    scar = ax.scatter(param1, param2, c=OH_MRC,s=90,marker='8', alpha=alpha,cmap=cm)
    cb = plt.colorbar(scar,cax=cax)
    cb.ax.tick_params(labelsize=16)
    cb.set_label('12+log(O/H)', fontweight='bold',fontsize=12)
    ax.set_xlabel(r'{}'.format(tag1),size=16,fontweight='bold')
    ax.set_ylabel(r'{}'.format(tag2),size=16,fontweight='bold')
    #ax.set_title("OH_out_NGC628_12",fontweight='bold',fontsize=18)
    ax.minorticks_on()
    ax.tick_params(axis='both',which='minor',length=5,width=2,labelsize=18)
    ax.tick_params(axis='both',which='major',length=11,width=2,labelsize=18)
    #ax.plot((),(),'k--',label='$\sigma$ = 0.09 dex')
    ax.legend(loc = 'best', numpoints=1, fancybox=True,fontsize=15)
    
    plt.savefig('/home/manga/Google Drive/negra/Maestria/sem4/seminario_Titula/Resultados_plots/CAbundan_plots/OH_M13-P12_5.pdf',format='pdf',dpi=100)
    
    plt.tight_layout()
    plt.show()
    plt.ion()
    return scar
  

sc1 = plot_OH(OH_Te,OH_MRC,'12+log(O/H)[Te]','12+log(O/H)[MRC]')

#sc1 = plot_OH(M_O3N2,OH_mean_out,'12+log(O/H) [O3N2]','12+log(O/H) [MRC]')
#sc1 = OH_plotterarino(M_O3N2,OH_Te,'12+log(O/H) [O3N2]','12+log(O/H) [Te]')
#sc1 = OH_plotterarino(np.log10(2.8*R3N2),OH_Te,'log(O3N2)','12+log(O/H) [Te]')
#sc1 = OH_plotterarino(np.log10(2.87*R3N2),OH_mean_out,'log(O3N2)','12+log(O/H) [out]')

#sc1 = OH_plotterarino(M_Pil,OH_mean_out,'12+log(O/H) [ONS]','12+log(O/H) [MRC]')
#sc1 = OH_plotterarino(M_Pil,OH_Te,'12+log(O/H) [ONS]','12+log(O/H) [Te]')
#sc1 = OH_plotterarino(np.log10(R3N2),M_Pil,'log(O3N2)','12+log(O/H) [ONS]')

#sc1 = OH_plotterarino(M_N2,OH_mean_out,'12+log(O/H) [N2]','12+log(O/H) [MRC]')
#sc1 = OH_plotterarino(M_N2,OH_Te,'12+log(O/H) [N2]','12+log(O/H) [Te]')
#sc1 = OH_plotterarino(np.log10((1/2.86)*N2),OH_Te,'log(N2)','12+log(O/H) [Te]')
#sc1 = OH_plotterarino(np.log10((1/2.86)*N2),OH_mean_out,'log(N2)','12+log(O/H) [out]')

