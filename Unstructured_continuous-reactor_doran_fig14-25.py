# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 11:55:49 2020

@author: Carlos Caicedo-Montoya
"""

#Este programa es usado para observar el comportamiento de la biomasa
#para diferentes tasas de dilución.
#Aquí se asume que la alimentación no tiene biomasa, y que la velocidd de 
#muerte es mucho menor que la velocidad de crecimiento.
#Tambien se asume que el producto es asociado al metabolismo energético
#o que no hay un producto asociado y que el mantenimiento es despreciable.

#Importar librerías
import numpy as np
import matplotlib.pyplot as plt

#Datos
Ks= 0.2     #kg/m3
umax = 0.8   #1/h
Si= 20      #kg /m3
Yxs = 0.5      #kg biomasa/kg sustrato

#Cálculo de la tasa de dilución crítica  
Dcrit = umax*Si/(Ks+Si)
print('el valor de la tasa de dilución crítico al cual inicia el lavado del \
      reactor es {0:0.3f}'.format(Dcrit), '1/h')


#Cálculo de la tasa de dilución óptima
Dopt = umax*(1-np.sqrt(Ks/(Ks+Si)))
print('el valor de la tasa de dilución óptima al cual se tiene la máxima \
      productividad de biomasa es {0:03.3f}'.format(Dopt), '1/h' )

D = np.linspace (0, Dcrit-0.001)
X = (Si-D*Ks/(umax-D))*Yxs


D1 = np.linspace (0, Dcrit+0.1)
S=[]

for i in D1:
    Scon = i*Ks/(umax-i)
    if Scon > Si:
        S.append(Si)
    elif Scon < 0:
        S.append(Si)
    else:
        S.append(Scon)


S = np.asarray(S)


Productividad_Biomasa = D*X
Max_Productividad_Biomasa = max(Productividad_Biomasa)


#Graficar los resultados
fig, ((ax1, ax2)) = plt.subplots(1, 2)
ax1.plot(D, X, '--g')
ax1.set_title('Biomasa')
ax1.set_xlabel('Tasa de Dilución (h$^{-1}$)')
ax1.set_ylabel('Concentración en estado estacionario (kg/m$^3$)')

ax2.plot(D1, S, 'y')
ax2.set_title('Sustrato')
ax2.set_xlabel('Tasa de Dilución (h$^{-1}$)')

fig2, ax = plt.subplots()
ax.plot(D, Productividad_Biomasa)
ax.axis([D[0], D[-1],
         Productividad_Biomasa[0], Max_Productividad_Biomasa])
ax.set_xlabel('Tasa de Dilución (h${^-1}$)')
ax.set_ylabel('Productividad volumétrica de biomasa (kg/m$^3$h)')
ax.annotate(r'D$_o$$_p$$_t$ = {0:3.2f}'.format(Dopt), xy=(Dopt, Max_Productividad_Biomasa),
            xycoords='data',
            xytext=(0.8, 0.6), textcoords='axes fraction',
            arrowprops=dict(facecolor='green', shrink=0.05),
            horizontalalignment='right', verticalalignment='top',
            )
ax.grid()

#fig.savefig("Fig125", format = )
#fig2.savefig("Fig124", dpi=300)

