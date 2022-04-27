# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 18:44:20 2022

@author: Carlos Caicedo-Montoya
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

#Import libraries
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

"""
A:Bromuro de cianógeno
B:MetilAmina
C:Metil Bromuro
D:Cianamida
"""


#Definir el modelo
def bromuro_etilo(t, C):
    CA, CB, CD, CC, V = C  # B - D- T- H    
    dCA = -k*CA*CB - Q*CA/V
    dCB = -k*CA*CB + Q*(CBi-CB)/V
    dCC = k*CA*CB - Q*CC/V
    dCD = k*CA*CB - Q*CD/V
    dV = Q
    return [dCA, dCB, dCC, dCD, dV]

#Data
Q=0.05     #flujo volumetrico dm3/s
k=2.2   #dm3/s.mol
Vinit=5 #dm3
CBi=0.025       #mol/dm3
CA0=0.05        #mol/dm3

#Input for the numerical problem
C_init = np.array([CA0, 0, 0, 0, Vinit])  
t_lim = (0, 500)         #min
t_array = np.linspace(0, 500, num = 1000)

#solve the EDO system
sol = solve_ivp(bromuro_etilo, t_span = t_lim, 
                y0= C_init, t_eval = t_array)

conversion = (CA0*Vinit-sol.y[0]*sol.y[4])/(CA0*Vinit)
CA = sol.y[0]
CB = sol.y[1]
CC = sol.y[2]
CD = sol.y[3]

#Plots
f1, ax = plt.subplots()
ax.plot(t_array, CA, label = r"$C_A$", color = "olive")
ax.plot(t_array, CB, label = r"$C_B$", color = "darkblue")
ax.plot(t_array, CC, label = r"$C_C$", color = "chocolate" )
#I do not include the curve for CD because it is exactly the same of CC
ax.grid()
ax.legend()
ax.set_xlabel("Tiempo (s)", size = 14)
ax.set_ylabel("Concentración (mol/L)", size = 14)
ax.set_title("Curvas Concentración-tiempo")

f2, ax = plt.subplots()
reaction_rate = k*CA*CB
ax.plot(t_array, reaction_rate)
ax.grid()
ax.set_xlabel("Tiempo (s)", size = 14)
ax.set_ylabel("Velocidad de reacción (mol/L.s)", size = 14)
ax.set_title("Curva Velocidad de reacción-tiempo")

