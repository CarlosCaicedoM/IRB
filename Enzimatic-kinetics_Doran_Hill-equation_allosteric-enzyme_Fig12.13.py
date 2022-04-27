# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:33:18 2020

@author: Carlos Caicedo-Montoya
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')



import numpy as np
import matplotlib.pyplot as plt

#Datos
vmax = 2.5*3600/1000    #mM/hr
Kh = 2.5       #mM^n

n = [0.5, 1, 2, 4]
S = np.linspace(0, 10)
r_S_array =np.zeros(len(n)*len(S)).reshape(len(n), len(S))

for i, j in enumerate(n):
    r_S_array[i] = vmax*(S**j)/(Kh + (S**j))

#Graficas
fig, ax = plt.subplots()
for i in range(len(r_S_array)):
    ax.plot(S, r_S_array[i], label='n={}'.format(n[i]))
ax.set_xlabel('Concentración (mM)')
ax.set_ylabel(r'$r_S \/ (mM \/ hr ^{-1})$')
ax.legend()
ax.grid()

