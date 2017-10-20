#!/usr/bin/env python
#-*- coding-utf-8 -*-


__autores__ = "A.A.Aroche & S.F. Sanchez"
__date__ = "2017"

#Abundances 

import sys
import csv
import argparse
import matplotlib.pyplot as plt
from scipy import stats
import numpy.ma as ma
from numpy import *
import pyfits
import itertools
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import polyfitter as plf

#=======================CONSTANTS=======================#

Re = 5.03   #Effective radius in kpc. Belley et al. (1992)
R25 = 14.1  #Optical radius in (kpc) at the B_25 mag arcsec(-2). Rosales-Ortega et al. (2010)
rho = 314   #arcsec #Size of the optical radius $\rho_{25}$ in arsec
A = 58      #arcsec
B = 5       #adimensionless
C = 44      #1/rads

#--------------------------------------------------------------------------#

#-------Command line input files--------------#

usage = 'Usage: %s infile' % sys.argv[0]

try:
    input_file = sys.argv[1]      #fluxes NGC628
    input_file = sys.argv[2]      #coordinates NGC628
    input_file = sys.argv[3]      #OH_abundances NGC628_Thesis
    input_file = sys.argv[4]      #coeficientes de orden n
except:
    print usage; sys.exit(1)


#--------------------Reading the fits file--------------------#

input_fits = sys.argv[4]
hdu_coeffs = pyfits.open(input_fits)

coeffs_out = hdu_coeffs[0].data
e_coeffs_out = hdu_coeffs[1].data

#---------------------HII regions catalogues----------------------#

#----------NGC628. Sanchez et al. (2011)-----------#
#Fluxes
#All normalized to Hb
input_file = sys.argv[1]

f = pd.read_csv(input_file, comment = '#', header=None)
#f1 = f.dropna(axis=0, how='all')

'''
Ha =  f1[9].astype(float)

R2  = f1[3].astype(float)
R3  = f1[5].astype(float)
NII  = f1[11].astype(float)
S2a  = f1[15].astype(float)
S2b  = f1[17].astype(float)
S2 = S2a+S2b
'''
Hb = ma.masked_invalid(f[1].astype(float))
#print Hb
Ha = ma.masked_invalid(f[9].astype(float))

R2  = ma.masked_invalid(f[3].astype(float))
R3  = ma.masked_invalid(f[5].astype(float))
NII  = ma.masked_invalid(f[11].astype(float))
S2a  = ma.masked_invalid(f[15].astype(float))
S2b  = ma.masked_invalid(f[17].astype(float))
S2 = S2a+S2b

N2 = NII/Ha 

R3N2 = R3/N2
R23 =  R2+R3
R3R2 = R3/R2
N2R2 = N2/R2
N2S2 = N2/S2

P = R3/(R23)

#=========Dust Correction========#
#Cardelli et al. 1989: A_lambda / A_v = a(x) + b(x)/R_v
#x = 1/lambda (10^-6 m)
#R = 3.1, Color excess, reddening of light

def Cardelli(x,R):
  y = (x-1.82)
  a=1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
  b=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
  extin = a + (b/R)
  return extin

#Coeficientes de correccion de Cardelli

beta_corr = Cardelli(1/0.4861,3.1)
alpha_corr = Cardelli(1/0.6563,3.1)


#Flujos corregidos de NGC628

def correc(flux,Lamb,H_alpha,H_beta):
  #rat_Ha_Hb = H_alpha/H_beta
  C = (log10(2.86) - log10(H_alpha/H_beta))/(alpha_corr - 1)
  lm = Lamb/1e4
  correccion = flux*10**(C*Cardelli(1/lm,3.1))

  return correccion


Hacor = correc(Ha,6562,Ha,1)
R2cor   = correc(R2,3727,Ha,1)
R3cor   = correc(R3,5007,Ha,1)
NIIcor   = correc(NII,6583,Ha,1)


S2acor  = correc(S2a,6717,Ha,1)
S2bcor  = correc(S2b,6731,Ha,1)

S2cor = S2acor+S2bcor
N2cor = NIIcor/Hacor

R3N2cor = R3cor/NIIcor
R23cor =  R2cor+R3cor
R3R2cor = R3cor/R2cor
N2R2cor = NIIcor/R2cor
N2S2cor = NIIcor/S2cor

Pcor = R3cor/(R23cor)


#------------------------------------------------------------------#

#----NGC628.HII.Coordinates-----#
#Radial Gradient and 2D maps distribution OH abundances

input_file = sys.argv[2]
h = pd.read_csv(input_file, comment = '#', header=None)

# delta right ascension from the center (arcsec)
# delta declination from the center

delta_RA = h[3]
delta_DEC =h[4]

#Galactocentric radius (kpc) 

R = h[7] 

#Normalized radius

radius = R/R25
radio = R/Re


#=========================CALIBRATORS. NGC628======================#

#--------O3N2 & N2 methods. Marino et al. (2013)-----------#

OHO3N2_Marino13sc = 8.533 - 0.214*log10(2.86*(R3/NII))
OHO3N2_Marino13c = 8.533 - 0.214*log10(2.86*(R3cor/NIIcor))
OHN2_Marino13c = 8.743 + 0.462*log10(NIIcor/2.86)
OHN2_Marino13sc = 8.743 + 0.462*log10(NII/2.86)

#OHO3N2_PP04   = 8.73 - 0.32*log10(R3cor/N2cor)

#OHN2_PP04   = 8.90 + 0.57*log10(NIIcor/2.86)

#np.savetxt('OHO3N2_Marino13c.csv',OHO3N2_Marino13c,fmt='%1.5f')
#np.savetxt('OHO3N2_Marino13sc.csv',OHO3N2_Marino13sc,fmt='%1.5f')
#np.savetxt('OHN2_Marino13c.csv',OHN2_Marino13c,fmt='%1.5f')


#-------ONS method. Pilyugin et al. (2010)-----------------#

OH_ONS = lambda a, R2, R3, NII, S2, P: (a[0] + a[1]*P + a[2]*np.log10(R3) + a[3]*np.log10(NII/R2)+ a[4]*np.log10(S2/R2))

#ONS coefficients

a_cool = (8.277, 0.657, -0.399, -0.061, 0.005)
a_warm = (8.816, -0.733, 0.454, 0.710, -0.337)
a_hot = (8.774, -1.855, 1.517, 0.304, 0.328)

#select regions (cool,warm,hot)

Pilyugin = []

for j in range(0,len(f[6])):
    
    if log10(NII[j]) > -0.1:
       ONS = OH_ONS(a_cool, R2[j], R3[j], NII[j], S2[j], P[j])

    if log10(NII[j]) < -0.1 and log10(NII[j]/S2[j]) > -0.25:
        ONS = OH_ONS(a_warm, R2[j], R3[j], NII[j], S2[j], P[j])
    
    if log10(NII[j]) < -0.1 and log10(NII[j]/S2[j]) < -0.25:
        ONS = OH_ONS(a_hot, R2[j], R3[j], NII[j], S2[j], P[j])

    Pilyugin.append(ONS)

    OHONSc = ma.masked_invalid(Pilyugin)

    OHONSc = ma.masked_where(OHONSc <= 8.0,OHONSc)
    OHONSc = ma.masked_where(OHONSc >= 8.8,OHONSc)


#MRC Abundances O/H, NGC628
#======================OH_calibrator_Tesis_NGC628=======================#

#Polyval2d evaluate a 2-D polynomial at points (x, y). This function returns the value. p(x,y) = \\sum_{i,j} c_{i,j} * x^i * y^j.


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
OH_R3N2_N2   = plf.polyval2d(log10(R3N2),log10(NII),c_R3N2_N2)
OH_R23_N2    = plf.polyval2d(log10(R23),log10(NII),c_R23_N2)
OH_R23_N2R2  = plf.polyval2d(log10(R23),log10(N2R2),c_R23_N2R2)
OH_N2R2_N2   = plf.polyval2d(log10(N2R2),log10(NII),c_N2R2_N2)
OH_N2_R3     = plf.polyval2d(log10(NII),log10(R3),c_N2_R3)
OH_R3N2_R3R2 = plf.polyval2d(log10(R3N2),log10(R3R2),c_R3N2_R3R2)
OH_R3R2_N2R2 = plf.polyval2d(log10(R3R2),log10(N2R2),c_R3R2_N2R2)
OH_R3R2_N2   = plf.polyval2d(log10(R3R2),log10(NII),c_R3R2_N2)
OH_N2S2_R3N2 = plf.polyval2d(log10(N2S2),log10(R3N2),c_N2S2_R3N2)
OH_N2S2_R23  = plf.polyval2d(log10(N2S2),log10(R23),c_N2S2_R23)
OH_N2S2_N2R2 = plf.polyval2d(log10(N2S2),log10(N2R2),c_N2S2_N2R2)
OH_N2S2_N2   = plf.polyval2d(log10(N2S2),log10(NII),c_N2S2_N2)

'''
OH_R3N2_R23  = plf.polyval2d(log10(R3N2cor),log10(R23cor),c_R3N2_R23)
OH_R3N2_N2   = plf.polyval2d(log10(R3N2cor),log10(NIIcor),c_R3N2_N2)
OH_R23_N2    = plf.polyval2d(log10(R23cor),log10(NIIcor),c_R23_N2)
OH_R23_N2R2  = plf.polyval2d(log10(R23cor),log10(N2R2cor),c_R23_N2R2)
OH_N2R2_N2   = plf.polyval2d(log10(N2R2cor),log10(NIIcor),c_N2R2_N2)
OH_N2_R3     = plf.polyval2d(log10(NIIcor),log10(R3cor),c_N2_R3)
OH_R3N2_R3R2 = plf.polyval2d(log10(R3N2cor),log10(R3R2cor),c_R3N2_R3R2)
OH_R3R2_N2R2 = plf.polyval2d(log10(R3R2cor),log10(N2R2cor),c_R3R2_N2R2)
OH_R3R2_N2   = plf.polyval2d(log10(R3R2cor),log10(NIIcor),c_R3R2_N2)
OH_N2S2_R3N2 = plf.polyval2d(log10(N2S2cor),log10(R3N2cor),c_N2S2_R3N2)
OH_N2S2_R23  = plf.polyval2d(log10(N2S2cor),log10(R23cor),c_N2S2_R23)
OH_N2S2_N2R2 = plf.polyval2d(log10(N2S2cor),log10(N2R2cor),c_N2S2_N2R2)
OH_N2S2_N2   = plf.polyval2d(log10(N2S2cor),log10(NIIcor),c_N2S2_N2)
'''
            
OH_array = np.array([OH_R3N2_R23,OH_R3N2_N2,OH_R23_N2,OH_R23_N2R2,OH_N2R2_N2,OH_N2_R3,OH_R3N2_R3R2,OH_R3R2_N2R2,OH_R3R2_N2,OH_N2S2_R3N2,OH_N2S2_R23,OH_N2S2_N2R2,OH_N2S2_N2])
          
OH_mean_out_628 = np.median(OH_array,axis=0)    #Oxygen Abundances outcomes
OHMRC628 = ma.masked_where(OH_mean_out_628 == 0,OH_mean_out_628)
OHMRC628_clip = np.clip(OHMRC628,7.8,8.9)
OHMRC628_clip = ma.masked_where(OHMRC628_clip <= 8.0,OHMRC628_clip)
OHMRC628_clip = ma.masked_where(OHMRC628_clip >= 8.9,OHMRC628_clip)



#--------------------OH_out_NGC628---------------------#

input_file = sys.argv[3]
c = pd.read_csv(input_file, comment='#',  header = None) 

n = c.dropna()

radios = n[0]
OH_MRC_628_corr = n[1]
OH_O3N2_628_corr = n[2]
OH_N2_628_corr = n[3]
OH_ONS_628_corr = n[4]



#=========================Plots===============================#
#NGC628
#RADIAL GRADIENTS O/H
#pil = ma.masked_invalid(Pilyu)
def radial(par1,par2,tag1,tag2,**args):
    
    fit_output = stats.linregress(par1,par2)   
    slope, intercept, r_value, p_value, std_err = fit_output
    
    line = slope*par1+intercept
    #print max(line)
    print 'O/H: {0:1f}, slope:{1:4f},stddev:{2:4f}'.format(intercept,slope,std_err)
    y_fit = polyval([slope,intercept],par1)
    #print intercept
    
    '''
    #=====Linear_fit=====#
    (slope,intercept) = polyfit(par1,par2,1)
    xr = polyval([slope,intercept],par1)
    print slope,intercept
    sc_fit = plt.plot(par1,xr,'r-',linewidth=3)
    '''

    alpha = 0.9
    cm = plt.cm.get_cmap('nipy_spectral')
   
    f1, ax1 = plt.subplots(figsize=(8,6))
    f1.patch.set_facecolor('white')
    
    sc1 = ax1.scatter(par1,par2,marker='o',c='k', s=50 ,label='OH_NGC628', alpha=alpha,cmap=cm)
    sc2 = ax1.plot(par1,y_fit,'r',label='Linear fit',linewidth=1)

    
    plt.xlabel(r'{}'.format(tag1),fontweight='bold',fontsize=25)
    plt.ylabel(r'{}'.format(tag2),fontweight='bold',fontsize=25)
    plt.xlim(0.0,0.8)
    plt.ylim(8.0,9.0)
    ax1.minorticks_on()
    plt.tick_params(axis='both',which='minor',length=5,width=2,labelsize=25)
    plt.tick_params(axis='both',which='major',length=11,width=2,labelsize=25)
    plt.legend(loc = 'best', numpoints=1, fancybox=True,fontsize=20)
    #plt.savefig('/home/manga/Google Drive/negra/Maestria/sem4/semi_grad/pyMRC_OH/OHNGC628/MRC_gradfit.pdf',format='pdf',dpi=100)
    
    plt.tight_layout()    
    plt.show()
    plt.ion()

   
    return sc1,sc2, slope, intercept


#sc_1 = radial(radios,OH_MRC_628_corr,'R/R$_{e}$','12 + log(O/H) [MRC]')

sc_1 = radial(radio,OHMRC628_clip,'R/R$_{e}$','12 + log(O/H) [MRC]')
#sc_1 = radial(radius,OHMRC628_clip,'R/R$_{25}$','12 + log(O/H) [MRC]')
#sc_2 = radial(radio,OHO3N2_Marino13sc,'R/R$_{e}$','12 + log(O/H) [O3N2]sc')
#sc_2 = radial(radio,OHO3N2_Marino13c,'R/R$_{25}$','12 + log(O/H) [O3N2]c')
#sc_3 = radial(radius,OHN2_Marino13c,'R/R$_{25}$','12 + log(O/H) [N2]M13')
#sc_3 = radial(radius,OHN2_Marino13sc,'R/R$_{25}$','12 + log(O/H) [N2]sc')
#sc_4 = radial(radius,Pilyugin,'R/R$_{25}$','12 + log(O/H) [ONS]')
#sc_4 = radial(radius,OH_ONS_628_corr,'R/R$_{25}$','12 + log(O/H) [ONS]')


#sc_3 = radial(radius,OHN2_PP04,'R/R$_{25}$','12 + log(O/H) [N2] PP04')
#=========================Plots===============================#

sys.exit()
#2D oxygen abundances maps
cm = plt.cm.get_cmap('nipy_spectral')
#cm = plt.cm.get_cmap('rainbow')
f1, ax1 = plt.subplots(figsize=(10, 8))
ax1 = plt.subplot(111, polar=False)
f1.patch.set_facecolor('white')
circ=plt.Circle((0,0), radius=rho, color='k',linestyle='dashed',linewidth=1,fill=False)
plt.gca().add_patch(circ)
#plt.axis('scaled')

ax1.axvline(x=0, c='k')
ax1.axhline(y=0, c='k')

ax1.set_xlim((250,-250))
ax1.set_ylim((-250,250))
#ax1.plot((6.5,9.5), (6.5,9.5), 'r',linewidth=1)
#size=60
alpha = 1

ax1 = plt.gca()
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="7%", pad=0.01)


#Analytical characterization of the spiral arms of NGC628
#S.F. Sanchez et al. (2012)

theta = np.arange(0.0,2*np.pi, 0.015)
r = A/(np.log(B*np.tan(theta/(2*C)))) 

ax1.plot(-3*r*np.cos(theta),3*r*np.sin(theta),'k', lw=1)
ax1.plot(3*r*np.cos(theta),-3*r*np.sin(theta), 'k',lw=1)

#sc = ax1.scatter(delta_RA,delta_DEC,c=OHMRC628_clip, s=200, marker='o',label='MRC', alpha=alpha,cmap=cm,vmin=8.0, vmax=8.9)
sc = ax1.scatter(delta_RA,delta_DEC,c=OHONSc, s=200, marker='o',label='ONS', alpha=alpha,cmap=cm,vmin=8.0, vmax=8.9)
#sc = ax1.scatter(delta_RA,delta_DEC,c=OHN2_Marino13sc, s=200, marker='o',label='NGC628', alpha=alpha,cmap=cm,vmin=8.0, vmax=8.9)
#sc = ax1.scatter(delta_RA,delta_DEC,c=OHMRC628_clip, s=200, marker='8',label='NGC628', alpha=alpha,cmap=cm,vmin=8.0, vmax=8.9)
cb = plt.colorbar(sc,cax=cax)
cb.ax.tick_params(labelsize=16)
#sc = ax1.scatter(delta_RA,delta_DEC,c=OH_out_628,s=20*OH_out_628,marker='8',label='NGC628', alpha=alpha,cmap=cm,vmin=8.0, vmax=8.9)


sc = ax1.scatter(50.6,46.9,c='k',s=250,marker='*',label='Croxal et al. (2013)', alpha=alpha,cmap=cm)
sc = ax1.scatter(-55.5,84,c='k',s=250,marker='*', alpha=alpha,cmap=cm)
sc = ax1.scatter(-37.6,112.7,c='k',s=250,marker='*', alpha=alpha,cmap=cm)
sc = ax1.scatter(41.7,-122.1,c='k',s=250,marker='*', alpha=alpha,cmap=cm)
sc = ax1.scatter(-59.2,-112.4,c='k',s=250,marker='*', alpha=alpha,cmap=cm)
sc = ax1.scatter(-40.2,-158.9,c='k',s=250,marker='*', alpha=alpha,cmap=cm)
sc = ax1.scatter(80.4,-138.8,c='k',s=250,marker='*', alpha=alpha,cmap=cm)





#cb.set_label(r'12+log(O/H)_Thesis', fontweight='bold',fontsize=12)
ax1.set_xlabel(r'$\Delta RA$ $[arcsec]$',fontsize=18)
ax1.set_ylabel(r'$\Delta Dec$ $[arcsec]$',fontsize=18)
ax1.text(180,203 , r'ONS', style='normal', fontsize= 20)
ax1.text(-205,-185 , r'R$_{25}$', fontsize= 20,style='oblique')
ax1.text(-10,228 , 'N', fontsize= 16,style='oblique',bbox={'facecolor':'none', 'alpha':0.4, 'pad':5})
ax1.text(240,-25 , 'E', fontsize= 16,style='oblique',bbox={'facecolor':'none', 'alpha':0.4, 'pad':5})
#ax1.text(-155,215 , 'n = ', fontsize= 16,style='oblique',bbox={'facecolor':'none', 'alpha':0.4, 'pad':5})
#ax1.set_title("MRC",fontweight='bold',fontsize=18)
ax1.minorticks_on()
ax1.tick_params(axis='both',which='minor',length=5,width=2,labelsize=18)
ax1.tick_params(axis='both',which='major',length=11,width=2,labelsize=18)


#ax1.legend(loc = 'best', numpoints=1, fancybox=False, fontsize = 12)
f1.tight_layout()
plt.show()
plt.ion()

