#!/usr/bin/env python
# -*- coding: utf-8 -*-


import csv
import sys
import itertools
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from matplotlib.mlab import griddata
import sympy
from polynomial2d import polyfit2d
import polyfitter as plf
import pandas as pd
import numpy.ma as ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import seaborn as sns
#sns.set(style="white")
import matplotlib.cm as cm

np.seterr(divide='ignore', invalid='ignore')     # divide: Treatment for division by zero &



#-------Command line input files--------------#

usage = 'Usage: %s infile' % sys.argv[0]

try:
    input_file = sys.argv[1]      #catalogo Marino OH_Te
    input_file = sys.argv[2]      #orden polinomial
except:
    print usage; sys.exit(1)


#---------------------HII regions catalogues----------------------#

#-----------Marino et al. (2013)---------------#
#Clipped catalogue
input_file = sys.argv[1]

g = pd.read_csv(input_file, comment = '#', header=None)


OII3727 = g[2]
OII3729 = g[3]
OIII4959 = g[5]
OIII5007 = g[6]
NII6548 = g[8]
NII6584 = g[9]
SII6717 = g[11]
SII6731 = g[12]

OH_Te = ma.masked_where(g[16] == 0,g[16])


#Proposed ratios
#Following the Pilyugin et al. (2010) notation

R2 = ma.masked_invalid(OII3727 + OII3729)
R3 = ma.masked_invalid(OIII4959 + OIII5007)
#R3 = ma.masked_invalid(OIII5007)
N2 = ma.masked_invalid(NII6548 + NII6584)
#N2 = ma.masked_invalid(NII6584)
S2 = ma.masked_invalid(SII6717 + SII6731)
#S2 = ma.masked_invalid(SII6717)

R3N2 = ma.masked_invalid(R3/N2)
R23  = ma.masked_invalid(R2+R3)
R3R2 = ma.masked_invalid(R3/R2)
N2R2 = ma.masked_invalid(N2/R2)
N2S2 = ma.masked_invalid(N2/S2)

#===============Polynomial order=================#

input_order = sys.argv[2]  #orden polinomial
order = int(sys.argv[2])

#========================================================#


#================Polynomial fitting======================#

def polyfit2D(x, y, z, order):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)

    return m
    

def polyval2D(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


#Aqui hacemos la funcion donde calculamos las abundancias y los diagramas de los 
#distintos cocientes.

def OH_poly(par1_in,par2_in,OH_in,tagx,tagy):

#Con esto eliminamos los infinitos
#enmascaramos

    par=par1_in*par2_in
    par1=par1_in[np.isfinite(par)]
    par2=par2_in[np.isfinite(par)]
    OH=OH_in[np.isfinite(par)]

#Dentro de esta funcion OH_poly para el c\'alculo de la abundancia, llamamos a las funciones polyfit2d y polyval2d
#polyfitd2 recibe como entrada los cocientes y la abundancia del catalogo, los ajusta y entrega m, que se requiere para el
#calculo de la nueva abundancia

    m = polyfit2D(par1,par2,OH,order) 

    OH_mod = polyval2D(par1, par2, m)
    #print OH_mod
    
#definimos limites para los ejes

    par1min = min(par1)
    par1max = max(par1)
    par2min = min(par2)
    par2max = max(par2)

  
#Con el nuevo valor de abundancias modificadas, hacemos la interpolacion
#definimos los ejes para la malla

    xi = np.linspace(par1min, par1max, 594)
    yi = np.linspace(par2min, par2max, 594)

#Hacemos una malla donde haremos la interpolacion (usando la tecnica del vecino natural) a los cocientes y abundancia, y aqui se calcula 
#las nuevas abundancias (OH ajustada e interpolada = OH_i_mod)

    OH_i_mod = griddata(par1, par2, OH_mod, xi, yi, interp='nn')
    

    cm = plt.cm.get_cmap('nipy_spectral')
    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')
    #CS = ax.matshow(OH_i_mod, extent=[par1min,par1max,par2max,par2min], aspect='auto')
    CS = ax.imshow(OH_i_mod, extent=[par1min,par1max,par2max,par2min], aspect='auto')
    CS1 = ax.scatter(par1, par2, marker='o', c=OH, s=OH*10, alpha = 0.65,cmap=cm)
    #CS2 = ax.scatter(PAR1_eval, PAR2_eval, marker='D', c='black', alpha = 0.85)
    ax.set_xlabel(r'log({})'.format(tagx),size=25,fontweight='bold')
    ax.set_ylabel(r'log({})'.format(tagy),size=25,fontweight='bold')
    #ax.set_title(r'log[{}]-log[{}]'.format(tagx,tagy),loc = 'center',fontweight='bold',fontsize=12)
    divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="7%", pad=0.01)
    #cb = fig.colorbar(CS,cax=cax) 
    #cb.set_label('12+log(O/H)', fontweight='bold',fontsize=25)
    #cb_num_format = "%f"
    #cb.ax.tick_params(labelsize=25)
    ax.set_xlim(par1min,par1max)
    ax.set_ylim(par2min,par2max)
    ax.minorticks_on()
    ax.tick_params(axis='both',which='minor',length=5,width=2,labelsize=25)
    ax.tick_params(axis='both',which='major',length=11,width=2,labelsize=25)
    plt.tight_layout()
    plt.show()
    plt.ion()
    #plt.savefig('/home/manga/Google Drive/negra/Maestria/sem4/seminario_Titula/Resultados_plots/ratios/{}_{}.pdf'.format(tagx,tagy),format='pdf',dpi=100)
    #plt.savefig('/home/manga/Google Drive/negra/Maestria/sem4/seminario_Titula/Resultados_plots/ratios_3/{}_{}.pdf'.format(tagx,tagy),format='pdf',dpi=100)
    
    return OH_mod

#Mostramos las graficas de la comparacion de los pares


#OH_R3N2_R23 = OH_poly(np.log10(R3N2),np.log10(R23),OH_Te,'R3N2','R23')
#OH_R3N2_N2  = OH_poly(np.log10(R3N2),np.log10(N2),OH_Te,'R3N2','N2')
#OH_R23_N2   = OH_poly(np.log10(R23),np.log10(N2),OH_Te,'R23','N2')
#OH_R23_N2R2 = OH_poly(np.log10(R23),np.log10(N2R2),OH_Te,'R23','N2R2')
#OH_N2R2_N2  = OH_poly(np.log10(N2R2),np.log10(N2),OH_Te,'N2R2','N2')
#OH_N2_R3    = OH_poly(np.log10(N2),np.log10(R3),OH_Te,'N2','R3')
#OH_R3N2_R3R2= OH_poly(np.log10(R3N2),np.log10(R3R2),OH_Te,'R3N2','R3R2')
#OH_R3R2_N2R2= OH_poly(np.log10(R3R2),np.log10(N2R2),OH_Te,'R3R2','N2R2')
#OH_R3R2_N2  = OH_poly(np.log10(R3R2),np.log10(N2),OH_Te,'R3R2','N2')
#OH_N2S2_R3N2= OH_poly(np.log10(N2S2),np.log10(R3N2),OH_Te,'N2S2','R3N2')
OH_N2S2_R23 = OH_poly(np.log10(N2S2),np.log10(R23),OH_Te,'N2S2','R23')
#OH_N2S2_N2R2= OH_poly(np.log10(N2S2),np.log10(N2R2),OH_Te,'N2S2','N2R2')
#OH_N2S2_N2  = OH_poly(np.log10(N2S2),np.log10(N2),OH_Te,'N2S2','N2')


               
'''
OH_cal = np.array([OH_R3N2_R23,OH_R3N2_N2,OH_R23_N2,OH_R23_N2R2,OH_N2R2_N2,OH_N2_R3,OH_R3N2_R3R2,OH_R3R2_N2R2,OH_R3R2_N2,OH_N2S2_R3N2,OH_N2S2_R23,OH_N2S2_N2R2,OH_N2S2_N2])



OH_cal = ma.masked_invalid(np.median(OH_cal,axis=0))    #Oxygen Abundances outcomes
Delta_mean_OH =  OH_Te - OH_cal
mean_Delta_OH = np.mean(Delta_mean_OH)
sigma_Delta_OH = np.std(Delta_mean_OH)

print 'offset O/H =', mean_Delta_OH, sigma_Delta_OH
'''
