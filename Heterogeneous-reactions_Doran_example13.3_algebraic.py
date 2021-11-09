# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 18:24:11 2020

@author: Carlos Caicedo Montoya

VELOCIDAD DE REACCIÓN PARA UNA ENZIMA LIBRE E INMOVILIZADA
Datos tomados de 
BIOPROCESS ENGINEERING PRINCIPLES SECOND EDITION PAULINE M. DORAN. 2013 Elsevier

La enzima invertasa es inmoilizada en una resina de intercambio iónico con 
un diámetro promedio de 1 mm. La cantidad de enzima en las perlas es medida
por un ensayo de proteínas como 0.05 kg m^3 de resina.
Un reactor de columna es empacado con 20 cm3 de perlas. 75 ml de una solución 
de glucosa 16mM es bombeada rapidamente a  través del reactor. 
En otro reactor, una cantidad identica de enzima libre es mezclada en el
mismo volumen de solución de  sacarosa. Los parametros cinéticos para la 
enzima inmovilizada se asumen iguales a los de la enzima libre.
Km = 8.8 mM
Número de recambio = 2.4*10^{-3} gmol genzima^{-1} s^{-1}
DA = 2*10^{-6} cm^2 s^{-1}

a) Cuaes es la velocidad de reacción de la enzima libre
b) cual es la velocidad de reacción de la enzima inmovilizada

sacarosa+ H20 ----> glucosa + fructosa 
"""

#importar librerías
import numpy as np
import matplotlib.pyplot as plt

#Datos
Km = 8.8  #mM = mol/m3
DA = 2.0e-06*((1/100)**2)    #m2/s
R = (1/2)*(1/1000)  #m
Conc_proteina =  0.05 #kg m^3 de resina.
Volumen_perlas = 20*((1/100)**3)    #m^3
Volumen_reaccion = 75*(1/1000)*(1/1000)     #m3    
masa_enzimas = Volumen_perlas*Conc_proteina      #kg
CAs = 16    #mM = mol/m3
numero_recambio = (2.4e-3)*1000    # mol Kg-enzima^{-1} s^{-1}
concentracion_enzima_libre = masa_enzimas/Volumen_reaccion   #kg m^{3}

vmax_l = numero_recambio*concentracion_enzima_libre   #molm^{-3} s^{-1}
vmax_i = numero_recambio*Conc_proteina      #mol m^{-3} de catalizador s^{-1}

#Velocidad de reacción enzima libre
rA_l = vmax_l*CAs/(Km+CAs)
total_rA_l = rA_l*Volumen_reaccion



#Factor de efectividad cinetica de orden cero
def modulo_thiele_cero_esfera(R, k0, DA, CAs):
    return (R/(3*np.sqrt(2)))*(np.sqrt(k0/(DA*CAs)))

thiele_cero_esfera = modulo_thiele_cero_esfera(R, vmax_i, DA, CAs)

def factor_efectividad_cero_esfera(thiele_cero_esfera):
    if thiele_cero_esfera  > 0 and thiele_cero_esfera <= 0.577:
        ni0 = 1
    else: 
        gamma = np.arccos(2/(3*thiele_cero_esfera**2)-1)
        ni0 = 1-(0.5 + np.cos((gamma+4*np.pi)/3))**3
    return ni0

factor_efectividad_cero_esfera = factor_efectividad_cero_esfera(thiele_cero_esfera)

#Factor de efectividad cinetica de primer cero
def modulo_thiele_primer_esfera(R, k1, DA):
    return (R/3)*(np.sqrt(k1/DA))

thiele_primer_esfera = modulo_thiele_primer_esfera(R, vmax_i/Km, DA)

def factor_efectividad_primer_esfera(thiele_primer_esfera):
    return (1/(3*thiele_primer_esfera**2))*(3*thiele_primer_esfera*(1/np.tanh(3*thiele_primer_esfera))-1)

factor_efectividad_primer_esfera = factor_efectividad_primer_esfera(thiele_primer_esfera)


#Factor de efectividad Michaelis-Menten
#Ecuación de Moo-Young and Kobayashi

factor_efectividad_mm = (factor_efectividad_cero_esfera+(Km/CAs)*factor_efectividad_primer_esfera)/(1+Km/CAs)


#Velocidad de reacción catalizador

rA_obs = factor_efectividad_mm*total_rA_l 


print('El factor de efectividad es{0:0.3f}'.format(factor_efectividad_mm))
print('La velocidad de reaccion de la enzima \
 libre es {0:0.8f}'.format(total_rA_l),  'mol/s')

print('La velocidad de reaccion de la enzima \
 inmovilizada es {0:0.8f}'.format(rA_obs),  'mol/s')
