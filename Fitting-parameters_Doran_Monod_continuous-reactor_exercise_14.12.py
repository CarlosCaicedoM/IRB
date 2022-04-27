# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 09:51:33 2020

@author: Carlos Caicedo-Montoya


Solución del Ejercicio 14.12 
Tomado de 
BIOPROCESS ENGINEERING PRINCIPLES SECOND EDITION PAULINE M. DORAN. 2013 Elsevier

Análisis cinético de bacterias biorremediadoras utilizando un quimiostato

Una cepa de la bacteria Ancylobacter capaz de crecer en 1,2-dicloroetano
se aisló de sedimentos en el río Rin. Esta bacteria se utilizaráa para la
biorremediación in situ de suelos contaminados con halógenos clorados. 
Los parámetros cinéticos para el organismo son determinados utilizando datos
obtenidos de cultivo de quimiostato. Se utiliza un fermentador de 1 litro con
corriente de alimentación que contiene 100 μM de 1,2-dicloroetano.
Las c0ncentraciones de sustrato en estado estacionario
se miden en función del caudal del quimiostato

Caulda (mL/h), Sustrato(μM)
10 17.4
15 25.1
20 39.8
25 46.8
30 69.4
35 80.1
50 100

(a) Determinar μmax and KS para este organismo.
(b) Determinar el caudal de operación máximo 
"""

#Importar librerías 
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


#Datos del Problema
caudal = np.array([10, 15, 20, 25, 30, 35, 50])  #ml/hr
S = np.array ([17.4, 25.1, 39.8, 46.8, 69.4, 80.1, 100]) #μM
Volumen = 1.0 #L
Volumen = Volumen*1000    #mL
D = caudal/Volumen    #1/hr
Si = 100        #μM


#PAra todas las linealizaciones se descarta el último punto
#Con el propósito de obtener un mejor ajuste
#Linealización Lineweaver Burk 
x1 = 1/S
y1 = 1/D
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1[0:-1], y1[0:-1])

umax_1 = 1/intercept1
Ks_1 = slope1*umax_1

x_linea1 = np.linspace(x1[0], x1[-2])
linea_LB = intercept1 + x_linea1 * slope1

#Linealización Eadie-Hofstee 
x2 = D
y2 = D/S
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2[0:-1], y2[0:-1])

Ks_2 = -1/slope2  
umax_2 =  intercept2 * Ks_2

x_linea2 = np.linspace(x2[0], x2[-2])
linea_EH = intercept2 + x_linea2 * slope2


#Linealización Langmuir, Haneer-Wolf
x3 = S
y3 = S/D
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(x3[0:-1], y3[0:-1])

umax_3 =  1/slope3
Ks_3 = umax_3*intercept3

x_linea3 = np.linspace(x3[0], x3[-2])
linea_LHW = intercept3 + x_linea3 * slope3


#Tasa de Dilución optima
Dopt = umax_1*(1-np.sqrt(Ks_1/(Ks_1 + Si)))
Fopt = Dopt*Volumen  #L/hr
Fopt = Fopt*1000  #ml/hr

#Tasa de dilución crítica
Dcrit = umax_1*Si/(Ks_1+Si)
Fcrit = Dcrit*Volumen   #L/hr
Fcrit = Dcrit*1000    #mL/hr


#Obtener gráficas de las linealizaciones
fig, ax1 = plt.subplots()
ax1.plot(x1[0:-1], y1[0:-1], '*r')
ax1.plot(x_linea1, linea_LB)
ax1.set_xlabel(r'$ \frac{1}{S} \ ( \mu M^{-1})$')
ax1.set_ylabel(r'$ \frac{1}{D} (hr) $')
ax1.set_title ('Gráfica de Lineweaver-Burk')
ax1.text(x1[-2], y1[1], r'$R^2 = {0:3.2f} $'.format(r_value1**2))
ax1.text(x1[-2], y1[1] - abs((y1[-1]-y1[0]))/len(y1),
         r'$y={0:3.2f} x + {1:3.2f}$'.format(slope1, intercept1))
ax1.grid()

fig, ax2 = plt.subplots()
ax2.plot(x2[0:-1], y2[0:-1], '*r')
ax2.plot(x_linea2, linea_EH)
ax2.set_xlabel(r'$ D \ ( hr^{-1})$')
ax2.set_ylabel(r'$ \frac{D}{S} ( hr.\mu M^{-1}) $')
ax2.set_title ('Gráfica de Eadie-Hofstee')
#ax2.text(x1[-2], y1[1], 
#ax2.text(x1[-2], y1[1] - 6, 
ax2.text(x2[1], y2[2], r'$R^2 = {0:3.2f} $'.format(r_value2**2))
ax2.text(x2[1], y2[2] - abs((y2[-1]-y2[0]))/len(y2), 
         r'$y={0:3.3f} x + {1:3.3f}$'.format(slope2, intercept2))
ax2.grid()

fig, ax3 = plt.subplots()
ax3.plot(x3[0:-1], y3[0:-1], '*r')
ax3.plot(x_linea3, linea_LHW)
ax3.set_xlabel(r'$ S \ ( \mu M)$')
ax3.set_ylabel(r'$ \frac{S}{D} (\mu M.hr^{-1}) $')
ax3.set_title ('Gráfica de Langmuir')
ax3.text(x3[-4], y3[2], r'$R^2 = {0:3.2f} $'.format(r_value3**2))
ax3.text(x3[-4], y3[2]-abs((y3[-1]-y3[0]))/len(y3), 
         r'$y={0:3.2f} x + {1:3.2f}$'.format(slope3, intercept3))
ax3.grid()





