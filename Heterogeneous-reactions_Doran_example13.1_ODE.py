# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 23:46:35 2020

@author: Carlos Caicedo-Montoya

Perfil de concentración para una enzima inmovilizada
Datos tomados de 
BIOPROCESS ENGINEERING PRINCIPLES SECOND EDITION PAULINE M. DORAN. 2013 Elsevier


Una enzima se inmoviliza en perlas de agarosa de 8 mm de diámetro. 
La concentración de enzima en las perlas son 0,018 kg de proteína 
por m3 de gel. Se sumergen diez perlas en una solución bien mezclada
que contiene sustrato a una concentración de 3,2*10^3 kg m^-3.
La difusividad efectiva del sustrato en el gel de agarosa es 2,1*10^-9 m2 s^-1. 
La cinética de la reacción enzimática se puede aproximar como de primer orden 
con constante de reacción de primer orden igual a  3,11*10^5 s^-1 por kg 
de proteína. 
Los efectos de transferencia externa de masa son insignificantes. 
Trace el perfil de concentración de sustrato en estado estacionario dentro de 
las  perlas en función del radio de las partículas.
"""

#importar librerías
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

#Definir modelo
def perfil_concentracion_enzima_inmovilizada(r, C):
    #Cambiar el valor de rA dependiendo del tipo de 
    #cinética a utilizar
    z = C[0]
    rA = k1*z
    dz = C[1]
    dz2 = (-2/r)*dz + rA/DA
    return np.array([dz, dz2])

def bc(ya, yb):    
    return np.array([ya[1]-0, yb[0]-CAs])

#Datos
R = (8/2)/1000   #m
DA=2.1e-09    #m2/s
CAs = 3.2 *(10**-3)   #kg/m3
Volumen_perla = (4/3)*np.pi*(R**3)   #m3
Enzima_perla = 0.018*Volumen_perla
k1 = 3.11e5 * Enzima_perla   #s^-1

#Condiciones de frontera

r = np.linspace(0.005e-3, R)
y = np.zeros((2, r.size))
y[0] = CAs

sol = solve_bvp(perfil_concentracion_enzima_inmovilizada, bc,
                r, y)

conc = sol.y[0]
radius = sol.x

fig, ax = plt.subplots()
ax.plot(radius, conc, 'or')
ax.grid()
#ax.set_xlim(0, 0.004)
#ax.set_ylim(0, 0.0035)
ax.set_xlabel('Radio $\/r (m)$')
ax.set_ylabel('Concentración de sustrato $C_A (kg \/m^{-3})$')

