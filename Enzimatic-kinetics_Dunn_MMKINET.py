# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 17:21:15 2020

@author: Carlos Caicedo-Montoya


Cinética enzimática
Ejemplo tomado de :

Rajiv Dutta
Fundamentals of Biochemical Engineering
Ane Books India, 2008

Modelo de la reacción enzimática
E + S --> ES; k1
ES--> E + S; k2
ES --> E + P; k3
"""

#Importar librerías
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Definir el modelo
#Modelo completo
def cinetica_enzimatica_completo(t, C):
    S, ES, P, E = C
    dS = -k1*E*S + k2*ES
    dES = k1*E*S - (k2 + k3)*ES
    dP = k3*ES
    dE = -k1*E*S + (k2+k3)*ES
    return [dS, dES, dP, dE]

#Modelo con aproximación del esatdo estacionario dES/dt = 0
def cinetica_enzimatica_EE(t, S):  #EE: aproximación del estado estacionario
    dS = -vmax*S/(Km + S)
    return [dS]


#Condiciones iniciales y entradas del método solve_ivp
S0 = 0.1      #mol/L
ES0 = 0.0      #mol/L
P0 = 0.0  
E0 = 0.01      #mol/L
t_span = (0, 120)    #s
C_init = np.array([S0, ES0, P0, E0])
t_array = np.linspace(0, 120, 1000)

#Datos del ejercicio
k1 = 40  #L/mol.s
k2 = 5  #s^-1    
k3 = 0.5 #s^-1      
vmax =  k3*E0   #mol/L.s        
Km = (k2+k3)/k1    #mol/L

#Resolver las ecuaciones diferenciales
sol = solve_ivp(cinetica_enzimatica_completo, t_span, C_init, t_eval = t_array)
conc_array = sol.y

sol_EE = solve_ivp(cinetica_enzimatica_EE, t_span,[S0], t_eval = t_array) 
                   #method = 'LSODA')
conc_array_EE = sol_EE.y
#Graficar los resultados

#Gráfica del modelo completo
fig, ax1 = plt.subplots()
ax1.plot(t_array, conc_array[1], 'm', label = 'Enzima-Sustrato')
ax1.plot(t_array, conc_array[3], 'y', label = 'Enzima')
ax1.legend(loc = 'center')
ax1.grid()

fig, ax2 = plt.subplots()
ax2.plot(t_array, conc_array[0], label = 'Sustrato')
ax2.plot(t_array, conc_array[2], label = 'Producto')
ax2.grid()
ax2.legend()


#Gráfica del modelo EE
fig, ax3 = plt.subplots()
ax3.plot(t_array, conc_array_EE[0], label = 'Sustrato')
ax3.grid()
ax3.legend()


#MM vs Modelo completo

fig, ax4 = plt.subplots()
ax4.plot(t_array, conc_array_EE[0], label = 'Sustrato-MM')
ax4.plot(t_array, conc_array[0], label = 'Sustrato-Modelo completo')
ax4.set_ylabel(r'$r_S  \/ (mol \/  L^{-1} s^{-1})$')
ax4.set_xlabel(r'Concentración $(mol \/ L^{-1})$')
ax4.grid()
ax4.legend()







    
    
    
    