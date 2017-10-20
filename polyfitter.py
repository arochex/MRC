#!/usr/bin/env python

import csv
import itertools
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from matplotlib.mlab import griddata
from astropy.table import Table
import sys
import random as rnd


def polyfit2d(x, y, z, order):    
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G,z)
    
    return m

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j 
        #print z
    return z


def OH_poly(par1_in,par2_in,OH_in,PAR1_eval,PAR2_eval,tagx,tagy,order):

 

    par=par1_in*par2_in
    par1=par1_in[np.isfinite(par)]
    par2=par2_in[np.isfinite(par)]
    OH=OH_in[np.isfinite(par)]


    m = polyfit2d(par1,par2,OH,order) 

    
    OH_mod = polyval2d(par1, par2, m)
    
    Delta_OH=np.abs(OH-OH_mod)/1.5 #propagacion del error

    par1min = min(par1[np.isfinite(par1)])
    par1max = max(par1[np.isfinite(par1)])
    par2min = min(par2[np.isfinite(par2)])
    par2max = max(par2[np.isfinite(par2)])

    xi = np.linspace(par1min, par1max, 500)
    yi = np.linspace(par2min, par2max, 500)

    OH_i_mod = griddata(par1, par2, OH_mod, xi, yi, interp='nn')
    OH_i_error = griddata(par1, par2, Delta_OH, xi, yi, interp='nn')
    
    #OH_i_mod = polyval2d(xi, yi, m)

    step_par1=(par1max-par1min)/200
    step_par2=(par2max-par2min)/200


    OH_eval = polyval2d(PAR1_eval,PAR2_eval,m)
    #print len(OH_eval)

    II=(PAR1_eval-par1min)/step_par1
    JJ=(PAR2_eval-par2min)/step_par2
    I=II.astype(int)
    J=JJ.astype(int)
    I[I<0]==0
    J[J<0]==0



#    print 'INDEX='+str(J)+','+str(I)
    e_OH_eval= OH_eval
     #OH_i_error[J,I]
    #OH_eval
    #OH_i_error[J,I]
  
    
    return OH_eval,e_OH_eval




def OH_poly_c(par1_in,par2_in,OH_in,PAR1_eval,PAR2_eval,tagx,tagy,order):


 

    par=par1_in*par2_in
    par1=par1_in[np.isfinite(par)]
    par2=par2_in[np.isfinite(par)]
    OH=OH_in[np.isfinite(par)]


    m = polyfit2d(par1,par2,OH,order) #entrega los coeficientes como resultados de un ajuste polinomial para cada par de cocientes, la abundancia conocida y el orden polinomial

    
    OH_mod = polyval2d(par1, par2, m)  #evalua los coeficientes con cada par de cocientes, y obtiene entonces la abundancia de salida
    
    Delta_OH=np.abs(OH-OH_mod)/1.5 #propagacion del error

    par1min = min(par1[np.isfinite(par1)])
    par1max = max(par1[np.isfinite(par1)])
    par2min = min(par2[np.isfinite(par2)])
    par2max = max(par2[np.isfinite(par2)])

    xi = np.linspace(par1min, par1max, 500)
    yi = np.linspace(par2min, par2max, 500)

    OH_i_mod = griddata(par1, par2, OH_mod, xi, yi, interp='nn')       #interpolacion
    OH_i_error = griddata(par1, par2, Delta_OH, xi, yi, interp='nn')
    
    #OH_i_mod = polyval2d(xi, yi, m)

    step_par1=(par1max-par1min)/200
    step_par2=(par2max-par2min)/200


    OH_eval = polyval2d(PAR1_eval,PAR2_eval,m)
    #print len(OH_eval)

    II=(PAR1_eval-par1min)/step_par1
    JJ=(PAR2_eval-par2min)/step_par2
    I=II.astype(int)
    J=JJ.astype(int)
    I[I<0]==0
    J[J<0]==0



#    print 'INDEX='+str(J)+','+str(I)
    e_OH_eval= OH_eval
     #OH_i_error[J,I]
    #OH_eval
    #OH_i_error[J,I]
  
    
    return OH_i_mod,OH_i_error,m

    
"""    
    fig, ax = plt.subplots()
  
    CS = ax.imshow(OH_i_mod, extent=[par1min,par1max,par2max,par2min], aspect='auto')
    CS1 = ax.scatter(par1, par2, marker='o', c=OH, s=OH*5, alpha = 0.65)
    CS2 = ax.scatter(PAR1_eval, PAR2_eval, marker='D', c='black', alpha = 0.85)
    ax.set_xlabel(r'log{}'.format(tagx),fontweight='bold')
    ax.set_ylabel(r'log{}'.format(tagy),fontweight='bold')
    ax.set_title(r'log{}-log{}'.format(tagx,tagy),loc = 'center',fontweight='bold',fontsize=12)
    #if plot_cb:
    cb = fig.colorbar(CS) 
    cb.set_label('12+log(O/H)', fontweight='bold',fontsize=12)
    ax.set_xlim(par1min,par1max)
    ax.set_ylim(par2min,par2max)
    plt.tight_layout()
    plt.grid(True)
    #plt.show()
    #fig.savefig('.png')
    plt.savefig('test.png')
""" 
