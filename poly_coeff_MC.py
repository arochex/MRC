#!/usr/bin/env python
# -*- coding: utf-8 -*-


#Calculo de coeficientes del polinomio y sus errores

#Calibrador de abundancia. En este script se hace una verficación del calibrador, a través de uan prueba de 
# chi cuadrada. Además de determinar el grado del polinomio que mejor se ajusta para la interpolación


#11.10.16

###___Readme_1___###
#Nuevo calibrador de abundancia basado en 
#n proyecciones de cocientes de líneas en un plano r
#m dimensional. Se realiza la verificación del calibrador
# mediante una prueba de chi2, cuando ésta empieza a fluctuar
#se escoge el orden del polinomio que ajustará al calibrador.
#Se implementa un MC (aprox se sigue una distro gaussiana)
#haciendo n iteraciones la funcion polinonial que es la que calcula
#la abundancia a través de las líneas ahora con un MC 
#Se haran 100, 1000, 10000 iteraciones, por cada iteracion se obtiene
#una tabla con los valores de rsm, std y el orden polinomial
#con el nuevo orden de polinomio obtenido después de haber hecho n iteraciones
#MC, se escoge el orden menor posible y se fija en la funcion  OH_poly para 
#hacer nuevamente n iteraciones MC y de este modo calcular la std.

# Fecha: 25.01.17
######____________README_2____________######
#Ya se tiene el orden polinomico n = 11, ahora se hacen n iteraciones MC para calcular los coeficientes de dicho polinomio para los 13 pares de combinaciones de cocientes de lineas y para los errores de los coeficientes.

#Nuevas pruebas MC, ahora ya pasandole los coeficientes y los pares de cocientes, para calcular las abundancias, las deltas de las abundancias para ver que tanto nos alejamos del valor dado, la stddev y la chi2

#---------------------Importo librerías de python -------------------------------
import sys
import csv
import string
import itertools
import pylab as P
import numpy as np
import polyfitter as plf
import matplotlib.pyplot as plt
import random as rnd
import pyfits
import argparse
from scipy import ndimage
from astropy.table import Table
from matplotlib.mlab import griddata
from scipy.stats import chisquare

#----------------------Métodos para calibrar la abundancia------------------------

def methods():
    np.seterr(divide='ignore', invalid='ignore')

    O2a = np.array(OII3727)
    O2b = np.array(OII3729)
    N2a = np.array(NII6548)
    N2b = np.array(NII6584)
    S2a = np.array(SII6717)
    S2b = np.array(SII6731)
    O3a = np.array(OIII4959)
    O3b = np.array(OIII5007)

    
    a = (O2a+O2b)                              
    b = (N2a+N2b)
    c = (S2a+S2b)
    d = (O3a+O3b)

   
    R3N2 = np.log10(d/b)
    R23  = np.log10(a+d)
    R3R2 = np.log10(d/a)
    N2S2 = np.log10(b/c)
    N2R2 = np.log10(b/a)
    N2   = np.log10(b)
    S2   = np.log10(c)
    R2   = np.log10(a)
    R3   = np.log10(d)


    return R3N2, R23, R3R2, N2S2, N2R2, N2, S2, R2, R3


#------------------- methods with MC--------------------#

def methods_MC():
    np.seterr(divide='ignore', invalid='ignore')
    O2a = np.array(OII3727)
    O2b = np.array(OII3729)
    N2a = np.array(NII6548)
    N2b = np.array(NII6584)
    S2a = np.array(SII6717)
    S2b = np.array(SII6731)
    O3a = np.array(OIII4959)
    O3b = np.array(OIII5007)

#    guess_error = 0.15

    O2a = O2a + O2a*rnd.gauss(0,guess_error)
    O2b = O2b + O2b*rnd.gauss(0,guess_error)
    N2a = N2a + N2a*rnd.gauss(0,guess_error)
    N2b = N2b + N2b*rnd.gauss(0,guess_error)
    S2a = S2a + S2a*rnd.gauss(0,guess_error)
    S2b = S2b + S2b*rnd.gauss(0,guess_error)
    O3a = O3a + O3a*rnd.gauss(0,guess_error)
    O3b = O3b + O3b*rnd.gauss(0,guess_error)
    
    

    a = (O2a+O2b)                              
    b = (N2a+N2b)
    c = (S2a+S2b)
    d = (O3a+O3b)
    

    R3N2 = np.log10(d/b)
    R23  = np.log10(a+d)
    R3R2 = np.log10(d/a)
    N2S2 = np.log10(b/c)
    N2R2 = np.log10(b/a)
    N2   = np.log10(b)
    S2   = np.log10(c)
    R2   = np.log10(a)
    R3   = np.log10(d)
    
    
    

    return R3N2, R23, R3R2, N2S2, N2R2, N2, S2, R2, R3

#-----------------------Catálogo de regiones HII (Marino et al. 2013)-----------------

#usage = 'Usage: %s infile' % sys.argv[0]
'''
try:
    input_file = sys.argv[1]
except:
    print usage; sys.exit(1)
'''
#parser = argparse.ArgumentParser()
#parser.add_argument("input_file", help="catalogo de entrada de regiones HII")

#args = parser.parse_args()

#print 'input file:',args.input_file

input_file = sys.argv[1]
entrada = open(input_file) #open the file for reading

#entrada = open('catalogo_RII_clean_Marinoetal2013.csv')

tabla = [0]

for fila in csv.reader(entrada):
    tabla.append(fila)
entrada.close()

OII3727 = []
OII3729 = []
NII6548 = []
NII6584 = []
SII6717 = []
SII6731 = []
OIII4959 = []
OIII5007 = []
AbO = []
color = []
size = []


for fila in range(1, len(tabla)):
    
    OII3727.append(float(tabla[fila][2]))
    OII3729.append(float(tabla[fila][3]))
    NII6548.append(float(tabla[fila][8]))
    NII6584.append(float(tabla[fila][9]))
    SII6717.append(float(tabla[fila][11]))
    SII6731.append(float(tabla[fila][12]))
    OIII4959.append(float(tabla[fila][5]))
    OIII5007.append(float(tabla[fila][6]))
    AbO.append(float(tabla[fila][16]))
    color.append(int(tabla[fila][14])-1)  
    size.append(float(tabla[fila][16]))  
    

abO = np.array(AbO)
categories = np.array(color)
colormap = np.array(['k','g','r'])      # COLUMN15: jT, integer, the index jT is 1 when the electron temperature 
categories2 = np.array(size)            #is derived from the auroral and nebular lines of [O III] = 'k', is 2 when the                                        
                                        #[N II] = 'g' lines are used, and 3 for the case of [S III] = 'r'



#--------------------Llamo a la función que contiene a los calibradores----------------------------


#R3N2, R23, R3R2, N2S2, N2R2, N2, S2, R2, R3 = methods_MC()
   
'''
#--------------------Líneas de emisión de Regiones HII. Galaxia UGC10043. Coba et al. 2015---------


l = np.genfromtxt('UGC10043.txt', names=True, comments='#', dtype=None, skip_header=1)


R3N2_ugc = np.log10(l['R3']/l['N2'])
R23_ugc  = np.log10(l['R2']+l['R3'])
R3R2_ugc = np.log10(l['R3']/l['R2'])
N2S2_ugc = np.log10(l['N2']/l['S2'])
N2R2_ugc = np.log10(l['N2']/l['R2'])
N2_ugc   = np.log10(l['N2'])
S2_ugc   = np.log10(l['S2'])
R2_ugc   = np.log10(l['R2'])
R3_ugc   = np.log10(l['R3'])
'''



#--------------------Abundancias calculadas con los calibradores propuestos----------------------------------------
#---------------------Interpoladas y ajustadas con un polinomio de orden ?-----------------------------------------
#                          Usando el catálogo de Marino et al. 2012

#De un test de MC y chi2 se obtuvo el orden polinomial (n=11), con esto, ahora se vuelven hacer j iteraciones de MC y fijo el orden polinomial para calcular los coeficientes de dicho polinomio para cada combinacion de pares (13) y sus respectivos errores. 


#---------------------Pido argumentos en la linea de comandos----------------------#
#--------------argumentos de orden de polinomio y numero de iteraciones------------#

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="catalogo de entrada de regiones HII")
parser.add_argument("guess_error", help="Error estimado")
parser.add_argument("order", help="orden del polinomio")
parser.add_argument("n_mc", help="numero de iteraciones MC")
#parser.add_argument("file_out.fits", help="cubo que contiene arreglos de nxnxm, donde n = order y m = 13 combinaciones de pares de cocientes ")

args = parser.parse_args()
print '\n'
#print 'input file:',args.input_file
print 'error en los cocientes:', args.guess_error
print 'orden del polinomio de 1 hasta:', args.order
print 'numero de  iteraciones:', args.n_mc

print '\n'
#print args.file_out.fits

#input_file = sys.argv[1]
guess_error = float(sys.argv[2])   #error propuesto para las lineas = 0.15
order = int(sys.argv[3])   #hasta el orden n
n_mc = int(sys.argv[4])    #numero de iteraciones MC

#file_out.fits = sys.args[3]

#-----------------------------------------------------------------------------------#

#n_mc = 26                            #numero de iteraciones MC
#order = 8                            #variamos el orden del polinomio
#c_R3N2_R23=np.zeros(n_poly)

for n in range(1,order+1):

    n_poly=(n+1)*(n+1)           #numero de coeficientes del polinomio, de esta forma porque en realidad tenemos una tira de numeros
    coeffs_mc=np.zeros((n_mc,13,n_poly)) #regresa un arreglo de la forma que se indica, 
    coeffs_out=np.zeros((13,n_poly))     #arreglo de coeficientes de salida
    e_coeffs_out=np.zeros((13,n_poly))   #arreglo de errores de los coeficientes de salida


    
    print 'order = '+str(n)+  ' ,test_final'
    '''
    #tabla_mc='chi_mc_coeffs_test5_'+str(nc)+'.txt'
    tabla_mc='chi_mc_coeffs_test5_'+str(n)+'.txt'
    tablachi2 = open(tabla_mc,'w')
    salida = csv.writer(tablachi2, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    salida.writerow(['#','order','mean_deltaOH', 'std', 'chi2'])
    '''

    for nc in range(0,n_mc):                         #ciclo para hacer las n iteraciones MC
        
        '''
        tabla_mc='chi_mc_coeffs_test1_'+str(nc)+'.txt'
    	tablachi2 = open(tabla_mc,'w')
    	salida = csv.writer(tablachi2, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    	salida.writerow(['#','order','mean_deltaOH', 'std', 'chi2'])
    	#print 'MC iteration = '+str(nc)	
	#print '\n'
        '''
        
        R3N2, R23, R3R2, N2S2, N2R2, N2, S2, R2, R3 = methods_MC()  
        #for n  in range(1,order+1):
        #con plf.OH_poly_c ahora pedimos tambien los coeficientes de cada par de cocientes
	
        OH_R3N2_R23,e_OH_R3N2_R23,c_R3N2_R23    = plf.OH_poly_c(R3N2,R23,abO,R3N2,R23,'O3N2','R23',n)     
        OH_R3N2_N2,e_OH_R3N2_N2,c_R3N2_N2       = plf.OH_poly_c(R3N2,N2,abO,R3N2,N2, 'O3N2','N2',n)	
        OH_R23_N2,e_OH_R23_N2,c_R23_N2          = plf.OH_poly_c(R23,N2,abO,R23,N2,'R23','N2_N2',n)
        OH_R23_N2R2,e_OH_R23_N2R2,c_R23_N2R2    = plf.OH_poly_c(R23,N2R2,abO,R23,N2R2,'R23','N2/O2',n)
        OH_N2R2_N2,e_OH_N2R2_N2, c_N2R2_N2      = plf.OH_poly_c(N2R2,N2,abO,N2R2,N2,'N2/O2','N2',n)
        OH_N2_R3,e_OH_N2_R3,c_N2_R3             = plf.OH_poly_c(N2,R3,abO,N2,R3,'N2','O3',n)
        OH_R3N2_R3R2,e_OH_R3N2_R3R2,c_R3N2_R3R2 = plf.OH_poly_c(R3N2,R3R2,abO,R3N2,R3R2,'O3N2','OIII/OII',n)     
        OH_R3R2_N2R2,e_OH_R3R2_N2R2,c_R3R2_N2R2 = plf.OH_poly_c(R3R2,N2R2,abO,R3R2,N2R2,'OIII/OII','N2/O2',n)
        OH_R3R2_N2,e_OH_R3R2_N2,c_R3R2_N2       = plf.OH_poly_c(R3R2,N2,abO,R3R2,N2,'OIII/OII','N2',n)
        OH_N2S2_R3N2,e_OH_N2S2_R3N2,c_N2S2_R3N2 = plf.OH_poly_c(N2S2,R3N2,abO,N2S2,R3N2,'N2/S2','O3N2',n)
        OH_N2S2_R23,e_OH_N2S2_R23,c_N2S2_R23    = plf.OH_poly_c(N2S2,R23,abO,N2S2,R23,'N2/S2','R23',n)
        OH_N2S2_N2R2,e_OH_N2S2_N2R2,c_N2S2_N2R2 = plf.OH_poly_c(N2S2,N2R2,abO,N2S2,N2R2,'N2/S2','N2/O2',n)
        OH_N2S2_N2,e_OH_N2S2_N2,c_N2S2_N2       = plf.OH_poly_c(N2S2,N2,abO,N2S2,N2,'N2/S2','N2',n)
            
        coeffs_mc[nc,0,:]=c_R3N2_R23
        coeffs_mc[nc,1,:]=c_R3N2_N2
        coeffs_mc[nc,2,:]=c_R23_N2
        coeffs_mc[nc,3,:]=c_R23_N2R2
        coeffs_mc[nc,4,:]=c_N2R2_N2
        coeffs_mc[nc,5,:]=c_N2_R3              #llenamos matrices con los coeficientes para las 13 combinaciones  
        coeffs_mc[nc,6,:]=c_R3N2_R3R2
        coeffs_mc[nc,7,:]=c_R3R2_N2R2
        coeffs_mc[nc,8,:]=c_R3R2_N2
        coeffs_mc[nc,9,:]=c_N2S2_R3N2
        coeffs_mc[nc,10,:]=c_N2S2_R23
        coeffs_mc[nc,11,:]=c_N2S2_N2R2
        coeffs_mc[nc,12,:]=c_N2S2_N2



        #OH_array = np.array([OH_R3N2_R23,OH_R3N2_N2,OH_R23_N2,OH_R23_N2R2,OH_N2R2_N2,OH_N2_R3,OH_R3N2_R3R2,OH_R3R2_N2R2,OH_R3R2_N2,OH_N2S2_R3N2,OH_N2S2_R23,OH_N2S2_N2R2,OH_N2S2_N2])
            
        #OH_mean = np.median(OH_array,axis=0)
        #Delta_mean_OH = abO-OH_mean
        #mean_Delta_OH_mean = np.mean(Delta_mean_OH)
        #sigma_Delta_OH_mean = np.std(Delta_mean_OH)
        #Chi = ((Delta_mean_OH)**2)/((abO*0.01)**2)#abO
        #chi_sum = np.sum(Chi)/(len(abO)-n-1)
        #print 'order = '+str(n)+' : mean_Delta_OH_mean = '+str(mean_Delta_OH_mean)+' +- '+str(sigma_Delta_OH_mean)+' chi_sum = '+str(chi_sum)
    
      
        #tablachi2.write('{}\t{}\t{}\t{}\n'.format(n, mean_Delta_OH_mean, sigma_Delta_OH_mean,chi_sum))
        # tablachi2.write( "%s\t%s\t%s\t%s\n" % (order, mean_Delta_OH_mean, sigma_Delta_OH_mean,chi_sum))

        # We close the loop

    coeffs_out=np.mean(coeffs_mc,axis=0)
    e_coeffs_out=np.std(coeffs_mc,axis=0)

#        print 'coeffs_mc = '+str(coeffs_mc)
#        print 'coeffs_out = '+str(coeffs_out)

            
        # Final test
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
    
    R3N2, R23, R3R2, N2S2, N2R2, N2, S2, R2, R3 = methods()    

    print 'test =',R3N2[0],R23[0],N2[0]
    #sys.exit(0)
        #for n  in range(1,order+1):
        
    OH_R3N2_R23 = plf.polyval2d(R3N2,R23,c_R3N2_R23)
    OH_R3N2_N2  = plf.polyval2d(R3N2,N2,c_R3N2_N2)
    OH_R23_N2 = plf.polyval2d(R23,N2,c_R23_N2)
    OH_R23_N2R2 = plf.polyval2d(R23,N2R2,c_R23_N2R2)
    OH_N2R2_N2 = plf.polyval2d(N2R2,N2,c_N2R2_N2)
    OH_N2_R3 = plf.polyval2d(N2,R3,c_N2_R3)
    OH_R3N2_R3R2 = plf.polyval2d(R3N2,R3R2,c_R3N2_R3R2)
    OH_R3R2_N2R2 = plf.polyval2d(R3R2,N2R2,c_R3R2_N2R2)
    OH_R3R2_N2 = plf.polyval2d(R3R2,N2,c_R3R2_N2)
    OH_N2S2_R3N2 = plf.polyval2d(N2S2,R3N2,c_N2S2_R3N2)
    OH_N2S2_R23 = plf.polyval2d(N2S2,R23,c_N2S2_R23)
    OH_N2S2_N2R2 = plf.polyval2d(N2S2,N2R2,c_N2S2_N2R2)
    OH_N2S2_N2 = plf.polyval2d(N2S2,N2,c_N2S2_N2)
    
            
    OH_array = np.array([OH_R3N2_R23,OH_R3N2_N2,OH_R23_N2,OH_R23_N2R2,OH_N2R2_N2,OH_N2_R3,OH_R3N2_R3R2,OH_R3R2_N2R2,OH_R3R2_N2,OH_N2S2_R3N2,OH_N2S2_R23,OH_N2S2_N2R2,OH_N2S2_N2])
                
    OH_mean = np.median(OH_array,axis=0)
    Delta_mean_OH = abO-OH_mean
    mean_Delta_OH_mean = np.mean(Delta_mean_OH)
    sigma_Delta_OH_mean = np.std(Delta_mean_OH)
    Chi = ((Delta_mean_OH)**2)/((abO*0.01)**2)#abO
    chi_sum = np.sum(Chi)/(len(abO)-n-1)
    print 'MC itera = '+str(nc)+' : mean_Delta_OH_mean = '+str(mean_Delta_OH_mean)+' +- '+str(sigma_Delta_OH_mean)+' chi_sum = '+str(chi_sum)

    primhdu = pyfits.PrimaryHDU(data=coeffs_out)   #a datacube for every order given (1<n<11)
    hdulist=pyfits.HDUList([primhdu])
    hdulist.append(pyfits.PrimaryHDU(e_coeffs_out))
    outputfile='coeffs_'+str(n)+'.fits'
    hdulist.writeto(outputfile,clobber=True)
    
    #sys.exit(0)

    f1, ax1 = plt.subplots(figsize=(9, 9))
    ax1.set_xlim((7,9))
    ax1.set_ylim(())
    ax1.plot((7,9), (7,9), '--r',linewidth=2)
    size=40
    alpha = 0.6
    ax1.scatter(abO,abO-OH_mean,s=15,label='OH_comparison', alpha=alpha)
    ax1.set_xlabel('12 + log(O/H)_Te')
    ax1.set_ylabel('log(O/H)_Te - log(O/H)_out')
    f1.tight_layout()
    plt.show()




        #print 'Test_Final: order = '+str(n)+' : mean_Delta_OH_mean = '+str(mean_Delta_OH_mean)+' +- '+str(sigma_Delta_OH_mean)+' chi_sum = '+str(chi_sum)
        #print '\n'
            
    #tablachi2.write('{}\t{}\t{}\t{}\n'.format(n, mean_Delta_OH_mean, sigma_Delta_OH_mean,chi_sum))  
                
#tablachi2.close()
            



#file_out = sys.argv[3]

#print str(coeffs_out)

#print str(e_coeffs_out)

