# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 18:03:03 2020

@author: Carlos Caicedo-Montoya
"""

#Este programa pretende mostrar como cambian la productividad y la 
#concentración de biomasa cuando se tiene una corriente de reciclo y una unidad
#de separación de biomasa
#Se asume la muerte celular despreciable
#no se tiene en cuenta el mantenimiento
#No hay producto o esta esta asociado al metabolismo energético


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


#Datos
Ks = 0.2     #kg/m3
umax = 0.8   #1/h
Si = 20      #kg /m3
Yxs = 0.5      #kg biomasa/kg sustrato
alpha = 0.4
beta = 3.0
Yxs = 0.5

#Cálculo de la tasa de dilución crítica
Dcrit, = fsolve(lambda x: 
                (Yxs/(1+alpha*(1-beta)))*(Si - (Ks*x*(1+alpha*(1-beta)))/(umax-x*(1+alpha*(1-beta)))), 1)

    
D = np.linspace (0, Dcrit-0.001)


#Ecuaciones para biomasa y Sustrato
S = (Ks*D*(1+alpha*(1-beta)))/(umax-D*(1+alpha*(1-beta)))
X = (Yxs/(1+alpha*(1-beta)))*(Si - S)
Productividad_Biomasa = D*X

#Calculo de la tasa de dilución óptima 
Dopt = D[np.where(Productividad_Biomasa==max(Productividad_Biomasa))]



#Graficar los resultados
fig, ax1 = plt.subplots()
ax1.plot(D, X, '--g', label= "Biomasa")
ax1.set_title(r'$\alpha = {0:0.2f},\/\beta = {1:0.2f} $'.format(alpha, beta))
   #\/ espacio entre las variables 
ax1.set_xlabel('Tasa de Dilución (h$^-$$^1$)')
ax1.set_ylabel('Concentración en estado estacionario (kg/m$^3$)')
ax1.grid()
ax1.legend(loc='right')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(D, Productividad_Biomasa, 'y', label = "Productividad")
ax2.set_ylabel('Productividad (kg/m$^3$h)')
ax2.legend(loc='center left')
