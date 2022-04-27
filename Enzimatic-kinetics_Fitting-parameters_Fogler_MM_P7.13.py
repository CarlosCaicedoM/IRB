# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 12:20:17 2020

@author: Carlos Caicedo-Montoya


Elements of Chemical Reaction Engineering 4th Edition
H.S Fogler
2006 Pearson Education, Inc.
 
Ejercicio p7-13
Los siguientes datos sobre levadura de panadería en un medio particular 
a 23.4 ° C y varias presiones parciales de oxígeno se obtuvieron:

PO2     QO2(no sulfanilamide)  QO2(20 mg sulfanilamide/ml added to medium)
0.0 0.0 0.0
0.5 23.5 17.4
1.O 33.0 25.6
1.5 37.5 30 8
2.5 42.0 36.4
3.5 43.0 39.6
5.0 43.0 40.0

PO2 : Presión Parcial de Oxígeno, mmHg
QO2: Velocidad de captación de oxígeno, uL de O2 por hora por mg de células
a. Calcular el QO2 máximo(Vmax), y la constante de Michaelis-Menten
Respuesta(Vmax = 52.63 uL O2/h. mg de células)
b. Determine el tipo de inhibición que la sulfanilamida causa sobre la
captación de oxígeno(Usando la gráfica de linweaver-Burk plot)

"""

from IPython import get_ipython
get_ipython().magic('reset -sf')


#Importar librerías
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats



#Datos del problema
PO2 = np.array([0.0, 0.5, 1.0, 1.5, 2.5, 3.5, 5.0])
QO2 = np.array([0.0, 23.5, 33.0, 37.5, 42.0, 43.0, 43.0])
QO2_sulfanilamide = np.array([0.0, 17.4, 25.6, 30.8, 36.4, 39.6, 40.0])
T = 23.4 +273   #

#Gráficas de la presión parcial de  óxigeno vs la velocidad de captación
fig, ax1 = plt.subplots()
ax1.plot(PO2, QO2_sulfanilamide, 'r', label = 'Con sulfanilamida')
ax1.plot(PO2, QO2, label = 'Sin sulfanilamida')
ax1.set_xlabel(r'$ P_{O_2} \/ ( mmHg)$')
ax1.set_ylabel(r'$ Q_{O_2} ( \mu L \/ de \/ O_2 \/ \/ hr^{-1} \/ mg\/de\/células ^{-1}) $')
ax1.legend()
ax1.grid()

#Linealización Lineweaver Burk 
x1 = 1/PO2[1:len(PO2)+1]
y1 = 1/QO2[1:len(PO2)+1]
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x1, y1)

vmax = 1/intercept1    #uLO2/hr.mg cell
Km = slope1*vmax       #mmHG

x_linea1 = np.linspace(x1[0], x1[-1])
linea_LB = intercept1 + x_linea1 * slope1

#Obtener gráficas de las linealizaciones
fig, ax1 = plt.subplots(figsize=(8, 4))
ax1.plot(x1, y1, '*r')
ax1.plot(x_linea1, linea_LB)
ax1.set_xlabel(r'$ \frac{1}{P_{O_2}} \/ ( mmHg^{-1})$', fontsize = 14)
ax1.set_ylabel(r'$ \frac{1}{Q_{O_2}} ( hr\/ mg\/de\/células \/ \mu L ^{-1}) $', 
               fontsize = 14 )
ax1.set_title ('Gráfica de Lineweaver-Burk')
ax1.text(x1[1], y1[2], r'$R^2 = {0:3.2f}$'.format(r_value1**2))
ax1.text(x1[1], y1[2] - abs((y1[-1]-y1[0]))/len(y1),
         r'$y={0:3.2f}x + {1:3.2f}$'.format(slope1, intercept1))
ax1.grid()


#Linealización Lineweaver Burk para el oxígeno con sulfanilamida
x2 = 1/PO2[1:len(PO2)+1]
y2 = 1/QO2_sulfanilamide[1:len(PO2)+1]
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2, y2)

vmax_sulfanilamida = 1/intercept2
Km_sulfanilamida = slope2*vmax_sulfanilamida

x_linea2 = np.linspace(x2[0], x2[-1])
linea_LB_sulfanilamida = intercept2 + x_linea1 * slope2


#Obtener gráficas de las linealizaciones
fig, ax1 = plt.subplots()
ax1.plot(x2, y2, '*r')
ax1.plot(x_linea1, linea_LB_sulfanilamida)
ax1.set_xlabel(r'$ \frac{1}{P_{O_2}} \/ ( mmHg^{-1})$', fontsize = 14)
ax1.set_ylabel(r'$ \frac{1}{Q_{O_2}} ( hr\/ mg\/de\/células \/ \mu L ^{-1}) $', 
               fontsize = 14 )
ax1.set_title ('Gráfica de Lineweaver-Burk')
ax1.text(x1[1], y1[2], r'$R^2 = {0:3.3f}$'.format(r_value2**2))
ax1.text(x1[1], y1[2] - abs((y1[-1]-y1[0]))/len(y1),
         r'$y={0:3.2f}x + {1:3.2f}$'.format(slope2, intercept2))
ax1.grid()

#Graficas juntas
fig, ax3 = plt.subplots()
ax3.plot(x_linea1, linea_LB, 'y', label = 'Sin sulfanilamida')
ax3.plot(x_linea1, linea_LB_sulfanilamida, 'r', label = 'Con sulfanilamida')
ax3.set_xlabel(r'$ \frac{1}{P_{O_2}} \/ ( mmHg^{-1})$', fontsize = 14)
ax3.set_ylabel(r'$ \frac{1}{Q_{O_2}} ( hr\/ mg\/de\/células \/ \mu L ^{-1}) $', 
               fontsize = 14 )
ax3.set_title ('Gráfica de Lineweaver-Burk')
ax3.grid()
ax3.legend()



