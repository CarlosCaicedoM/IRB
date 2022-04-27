# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:38:23 2019

@author: Carlos Caicedo-Montoya

Resolver el caso de un reactor PFR, con un sistema de multiples reacciones
Pirólisis de etano en presencia de un óxido nítrico

Chemical Reactor
Analysis and Design
Fundamentals
James B. Rawlings
John G. Ekerdt
2002

Example 4.7: Ethane pyrolysis in the presence of NO

"""
#Etano+oxido nítrico <--> etil-redical+ HNO
#etil-radical --> H + C2H4
#H+Etano --> C2H5 + H2
#H+NO<--> HNO
#
#
#Etano=A
#Oxido nitrico=B
#Radical etil = C
#HNO=D
#H=E    
#C2H4=F
#H2=G

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

#Definir el modelo

def multiple_reactions_ethane(V, N, T):
    NA, NB, NC, ND, NE, NF, NG = N  # B - D- T- H 
    Ntot = NA + NB + NC + ND + NE + NF + NG   
    
    Ctot=P/(R1*T)
    CA = (NA/Ntot)*Ctot
    CB = (NB/Ntot)*Ctot
    CC = (NC/Ntot)*Ctot
    CD = (ND/Ntot)*Ctot
    CE = (NE/Ntot)*Ctot

    r1 = k1*CA*CB - k_1*CC*CD
    r2 = k2*CC
    r3 = k3*CA*CE
    r4 = k4*CB*CE - k_4*CD
   
    
    dNA_dV = -r1-r3
    dNB_dV = -r1-r4
    dNC_dV = r1-r2+r3
    dND_dV = r1+r4
    dNE_dV = r2-r3-r4
    dNF_dV = r2
    dNG_dV = r3
    return [dNA_dV, dNB_dV, dNC_dV, dND_dV, dNE_dV, dNF_dV, dNG_dV]


# Datos del problema
Q0 = 600   #L/s
T=1100    #K
P=1       #atm
R=8.3144   #J/K*mol    
R1 = 82.057 	#cc3-atm/mol-K
k1 = 1.0*10**(14)*math.exp(-217600/(R*T))  #[cm^3/mol.s
k_1 = 1.0*10**(12)*math.exp(0/(R*T))    #cm^3/mol.s
k2 = 3.0*10**(14)*math.exp(-165300/(R*T))  #s-1
k3 = 3.4*10**(12)*math.exp(-28500/(R*T))  #cm^3/mol.s
k4 = 1.0*10**(12)*math.exp(0/(R*T))     #cm^3/mol.s
k_4 = 1.0*10**(13)*math.exp(-200800/(R*T))   #s-1
yAf  = 0.95
yBf = 0.05
CAf = yAf*P/(R1*T)  #gmole/cm3
CBf = yBf*P/(R1*T)   #gmole/cm3
NA0 = CAf*Q0      #gmole/sec
NB0 = CBf*Q0      #gmole/sec


# Argumentos de entrada al método numérico
     
N_init = np.array([NA0, NB0, 0, 0, 0, 0, 0])                           #Conce
V_lim = (0, 1500)
V_array = np.linspace(0, 1500, num=500)

# Resolver el modelo
sol = solve_ivp(multiple_reactions_ethane, 
        V_lim,
        N_init,
        method='BDF', args=(1100,))

V_array = sol.t
conc_array = sol.y


#Graficar perfiles de concentración
figure2, axis = plt.subplots()
axis.plot(V_array, conc_array[0], label = "Etano", color = "goldenrod")
axis.plot(V_array, conc_array[1], label = "NO", color = "steelblue")
axis.plot(V_array, conc_array[6], label = "Etano", color = "maroon")
#lt.plot(V_array, conc_array[3])
axis.legend()
axis.set_xlabel(r'$Volumen (cm ^3)$')
axis.set_ylabel(r'$F_j \/ (mol \/\/\/ s^{-1})$')
axis.grid()
