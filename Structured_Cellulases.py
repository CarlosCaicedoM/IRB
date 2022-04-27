# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 13:04:17 2019

@author: Carlos Caicedo-Montoya

Mathematical Modeling and
Optimization of Cellulase Protein
Production Using Trichoderma
reesei RL-P37
Producción de celulasas por Trichoderma ressei

Arun Tholudur, W. Fred Ramirez, James D. McMillan
Biotechnol Bioeng. 1999;66(1):1-16.

ESQUEMA DE REACCIÓN
inoculo->micelio primario       uL, uZ
micelio primario --> micelio secundario       k1
micelio primario --> células muertas     kd1
micelio secundario --> Celulasas (Producto)    rP
micelio secundario --> esporas      k2
micelio secundario  --> células muertas    kd2

"""
from IPython import get_ipython
get_ipython().magic('reset -sf')


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Definir el modelo

def celulases(t, C):
    Z, L, Xp, Xs, Xo, Xt, P = C      
    uL = (umaxL*L/(KSL+L))*(KIL/(KIL+Z))
    uZ = (umaxZ*Z)/(KSZ+Z)
    
    if t>= 24:
        rp= alfa*uL+beta 
    else:
        rp = 0

    dZ = -(uZ)*Xp/YZ
    dL = (-(uL)*Xp)/YL
    dXp = (uL+uZ-k1-kd1)*Xp
    dXs = k1*Xp-(k2+kd2)*Xs
    dXo = k2*Xs
    dXt = (uL+uZ)*Xp-kd1*Xp-kd2*Xs
    dP = rp*Xs
    return np.array([dZ, dL, dXp, dXs, dXo, dXt, dP])


# Datos del problema

umaxL = 0.1008    # {%units h^-1}
KSL = 0.0194    # {%units g/L}
KIL = 6.502  #{%units g/L}
YL = 0.4106  #{%units g/g}
umaxZ = 0.2016 #{%units h^-1}
KSZ = 4.993   #{%units g/L}
YZ = 0.6418  #{%units g/g} 
k1 = 0.0610    #{% units h^-1}
k2 = 0.0069  #{% units h^-1}
kd1 = 0.0099   #{% units h^-1}
kd2 = 0.0089  #{%units h^-1}

alfa = 0.0525
beta = 0.0025 #{%units h^-1}
 

#Condiciones iniciales
Z0 = 8      #g/L
L0 = 42     #g/L
Xp0 = 1     #g/L
Xs0 = 0     #g/L
Xo0 = 0     #g/L
Xt0 = 1     #g/L
P0 = 0.55   #g/L


# Argumentos de entrada al método numérico
     
C_init = np.array([Z0, L0, Xp0, Xs0, Xo0, Xt0, P0])                           #Conce
t_lim = (0, 120)
t_array = np.linspace(0, 120, num=1000)

# Resolver el modelo
sol = solve_ivp(celulases, 
        t_lim,
        C_init,
        t_eval=t_array, 
        method ='Radau')

t_array = sol.t
conc_array = sol.y


#Graficar perfiles de concentración de biomasa
figure1, axis = plt.subplots()
plt.plot(t_array, conc_array[2], label = "Micelio primario")
plt.plot(t_array, conc_array[3], label = "Micelio secundario")
plt.plot(t_array, conc_array[4], label = "Esporas")
plt.plot(t_array, conc_array[5], label = "Biomasa total")

axis.legend()
axis.set_xlabel(r'$Tiempo (h)$')
axis.set_ylabel(r'$Concentración (g \/ L^{-1})$')
axis.grid()


#Graficar perfiles de concentración de sustratos y producto
figure2, axis = plt.subplots()
axis.plot(t_array, conc_array[0], label = "Xilosa")
axis.plot(t_array, conc_array[6], label = "Celulasas")

axis.legend()
axis.set_xlabel(r'$Tiempo (h)$')
axis.set_ylabel(r'$Xilosa; Celulasas \/ \/\/ (g \/ L^{-1})$')
ax2 = axis.twinx() 
ax2.plot(t_array, conc_array[1], label = "Lactosa", color = "darkred")
ax2.legend(loc = "center right")
ax2.set_ylabel(r'$Lactosa \/\/(g \/ L^{-1})$')


#Calcular volumen óptimo para alcanzar una productividad de 500 kg/año
Productividad_1 = (conc_array[6]-0.55)/(t_array + 0.0001)
Produc_max_1 = np.amax(Productividad_1)
Product_deseada_1 = 500*(1000)*(1/365)*(1/24) 
Volumen_reactor_direct = Product_deseada_1/Produc_max_1 

def  vol_opt (Vol_reactor, tol) :
    d = 100000
    while d >= tol:
        Vol_reactor = Vol_reactor + 1
        Productividad = (conc_array[6]-0.55)*Vol_reactor/(t_array + 0.0001)     #g/L*tiempo
        Product_deseada = 500*(1000)*(1/365)*(1/24)    #g/hr
        Produc_max = np.amax(Productividad)
        d = Product_deseada - Produc_max
    else : 
        return Vol_reactor


#Graficar Productividad
Vol_opt_reactor=vol_opt(1, 0.001)        

Productividad = (conc_array[6]-0.55)*Vol_opt_reactor/(t_array + 0.0001)     #g/L*tiempo
Figure3, axis = plt.subplots()
axis.plot(t_array, Productividad)
axis.set_xlabel(r'$tiempo (h)$')
axis.set_ylabel(r'$Productividad \/(g \/hr^{-1})$')
axis.grid()

#Gráficar Volumen como función de la productividad
Produc_start = Productividad_1[np.where(Productividad_1 > 0 )]
Volumen = Product_deseada_1/(Produc_start)
Figure5, axis = plt.subplots()
axis.plot(Produc_start, Volumen, 'y')
axis.set_xlabel('Productividad (g/L.hr)')
axis.set_ylabel('Volumen (L)')
axis.grid()



#Graficar alrededor de el punto óptimo para ver el volumen mínimo 

Produc_zoom = Produc_start[200:-1]
Volumen_zoom = Product_deseada_1/(Produc_zoom)
Volumen_zoom2 = Volumen [200:-1]
Figure6, axis = plt.subplots()
axis.plot(Produc_zoom, Volumen_zoom2, 'b*')
axis.set_xlabel('Productividad (g/L.hr)')
axis.set_ylabel('Volumen (L)')
axis.grid()

#Valor del tiempo del lote y del vollumen mínimo
t_optim = t_array[np.where(Productividad == np.amax(Productividad))]/24   #dias
print  ("El volumen del reactor es {0:0.2f} litros".format( Vol_opt_reactor))
print ("El tiempo de reacción de un lote es de {0:0.2f} días".format(float(t_optim)))
