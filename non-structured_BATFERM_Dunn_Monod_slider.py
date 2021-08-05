#%% Cambiando solo mu
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 17:35:23 2019

@author: Carlos Caicedo-Montoya

Cinetica de producción -modelo monod-reactor batch

I. J. Dunn, E. Heinzle, J. Ingham, J. E. Pfenosil
Biological Reaction Engineering
Dynamic Modelling Fundamentals
with Simulation Examples
Second, Completely Revised Edition

Batch Fermentation (BATFERM)

"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button



#Definir el modelo
def batch_monod(t, C, UM):
    X, S, P = C  
    if S > 0 :
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X       #SUBSTRATE RATE EQUATION, kg/m3 h Producto asociado 
                                    #directamente con el metabolismo energetico
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = RS   
        dP_dt = RP   
    else:
        S = 0    
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X      #SUBSTRATE RATE EQUATION, kg/m3 h
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = 0  
        dP_dt = RP       
    return [dX_dt, dS_dt, dP_dt]


#Datos del modelo
UM=0.3     #1/h
Kd=0.0001      # 1/h
ks=0.000001    #g/L
Ypx=7.7     #gP/gX h
mp=1.1   #gP/gX h
Yxs=0.06    #gX/gS
ms= 2.2         #gX/gs
Xo=0.1      #g/L
So=12       #g/L
Po=0    #g/L


#Condiciones Iniciales
C_init = np.array([Xo, So, Po])  
t_lim = (0, 10)         #min
t_array = np.linspace(0, 10, num=100)

# Resolver el modelo
sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array, args=(0.3,))

t_array = sol.t
conc_array = sol.y
Biomasa = conc_array[0]
Sustrato = conc_array[1]
Producto = conc_array[2]

#Graficar perfiles de concentración
figure1, ((ax1, ax2, ax3)) = plt.subplots(3, 1)
line1, = ax1.plot(t_array, Biomasa, label = "Biomasa", color = "darkorange")
ax1.legend()
ax1.set_xlabel(r'$tiempo (h)$')
ax1.set_ylabel(r'$[X] \/  (g \/ L ^{-1})$')
ax1.grid()


line2, = ax2.plot(t_array, Sustrato, label = "Sustrato", color = "darkblue")
ax2.legend()
ax2.set_xlabel(r'$tiempo (h)$')
ax2.set_ylabel(r'$[S] \/ (g \/ L ^{-1})$')
ax2.grid()

line3, = ax3.plot(t_array, Producto, label = "Producto", color = "red")
ax3.legend()
ax3.set_xlabel(r'$tiempo (h)$')
ax3.set_ylabel(r'$[P] \/ (g \/ L ^{-1})$')
ax3.grid()



# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.3)
axcolor = 'lightcoral'
ax_miu = plt.axes([0.05, 0.5, 0.15, 0.05], facecolor=axcolor)

s_miu= Slider(ax = ax_miu, 
             label = r'$\mu$ ($h ^ {-1}$)',
             valmin = 0.1, 
             valmax = 1, 
             valinit=UM,
             valfmt='%1.1f')


def update(val):
    UM= s_miu.val
    sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array, args=(UM, ))
    conc_array = sol.y
    Biomasa = conc_array[0]
    Sustrato = conc_array[1]
    Producto = conc_array[2]
    
    line1.set_ydata(Biomasa)
    line2.set_ydata(Sustrato)
    line3.set_ydata(Producto)
    
    figure1.canvas.draw_idle()

# register the update function with each slider
s_miu.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.1, 0.84, 0.09, 0.04])
button = Button(resetax, 'Reset variables', color='cornflowerblue', hovercolor='0.975')

def reset(event):
    s_miu.reset()
button.on_clicked(reset)


#%% cambiando la velocidad de crecimiento y la afinidas por el sustrato



import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button



#Definir el modelo
def batch_monod(t, C, UM, ks):
    X, S, P = C  
    if S > 0 :
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X       #SUBSTRATE RATE EQUATION, kg/m3 h Producto asociado 
                                    #directamente con el metabolismo energetico
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = RS   
        dP_dt = RP   
    else:
        S = 0    
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X      #SUBSTRATE RATE EQUATION, kg/m3 h
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = 0  
        dP_dt = RP       
    return [dX_dt, dS_dt, dP_dt]


#Datos del modelo
UM=0.3     #1/h
Kd=0.0001      # 1/h
ks=0.000001    #g/L
Ypx=7.7     #gP/gX h
mp=1.1   #gP/gX h
Yxs=0.06    #gX/gS
ms= 2.2         #gX/gs
Xo=0.1      #g/L
So=12       #g/L
Po=0    #g/L


#Condiciones Iniciales
C_init = np.array([Xo, So, Po])  
t_lim = (0, 10)         #min
t_array = np.linspace(0, 10, num=100)

# Resolver el modelo
sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array, args=(0.3, 0.000001))

t_array = sol.t
conc_array = sol.y
Biomasa = conc_array[0]
Sustrato = conc_array[1]
Producto = conc_array[2]

#Graficar perfiles de concentración
figure1, ((ax1, ax2, ax3)) = plt.subplots(3, 1)
line1, = ax1.plot(t_array, Biomasa, label = "Biomasa", color = "darkorange")
ax1.legend()
ax1.set_xlabel(r'$tiempo (h)$')
ax1.set_ylabel(r'$[X] \/  (g \/ L ^{-1})$')
ax1.grid()


line2, = ax2.plot(t_array, Sustrato, label = "Sustrato", color = "darkblue")
ax2.legend()
ax2.set_xlabel(r'$tiempo (h)$')
ax2.set_ylabel(r'$[S] \/ (g \/ L ^{-1})$')
ax2.grid()

line3, = ax3.plot(t_array, Producto, label = "Producto", color = "red")
ax3.legend()
ax3.set_xlabel(r'$tiempo (h)$')
ax3.set_ylabel(r'$[P] \/ (g \/ L ^{-1})$')
ax3.grid()



# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.3)
axcolor = 'lightcoral'
ax_miu = plt.axes([0.05, 0.5, 0.15, 0.05], facecolor=axcolor)
ax_ks = plt.axes([0.05, 0.3, 0.15, 0.05], facecolor=axcolor)


s_miu= Slider(ax = ax_miu, 
             label = r'$\mu$ ($h ^ {-1}$)',
             valmin = 0.1, 
             valmax = 1, 
             valinit=UM,
             valfmt='%1.1f')
s_ks= Slider(ax = ax_ks, 
             label = r'$K_s$ ($g \/L ^ {-1}$)',
             valmin = 0.000001 , 
             valmax = 5, 
             valinit=ks,
             valfmt='%1.1f')

def update(val):
    UM= s_miu.val
    ks = s_ks.val
    sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array, args=(UM, ks))
    conc_array = sol.y
    Biomasa = conc_array[0]
    Sustrato = conc_array[1]
    Producto = conc_array[2]
    
    line1.set_ydata(Biomasa)
    line2.set_ydata(Sustrato)
    line3.set_ydata(Producto)
    
    figure1.canvas.draw_idle()

# register the update function with each slider
s_miu.on_changed(update)
s_ks.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.1, 0.84, 0.09, 0.04])
button = Button(resetax, 'Reset variables', color='cornflowerblue', hovercolor='0.975')

def reset(event):
    s_miu.reset()
    s_ks.reset()
button.on_clicked(reset)


#%% Cambiando todos los parámetros

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button



#Definir el modelo
def batch_monod(t, C, UM, Kd, ks, Ypx, mp, Yxs, ms):
    X, S, P = C  
    if S > 0 :
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X       #SUBSTRATE RATE EQUATION, kg/m3 h Producto asociado 
                                    #directamente con el metabolismo energetico
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = RS   
        dP_dt = RP   
    else:
        S = 0    
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X      #SUBSTRATE RATE EQUATION, kg/m3 h
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = 0  
        dP_dt = RP       
    return [dX_dt, dS_dt, dP_dt]


#Datos del modelo
UM=0.3     #1/h
Kd=0.0001      # 1/h
ks=0.000001    #g/L
Ypx=7.7     #gP/gX h
mp=1.1   #gP/gX h
Yxs=0.06    #gX/gS
ms= 2.2         #gX/gs
Xo=0.1      #g/L
So=12       #g/L
Po=0    #g/L


#Condiciones Iniciales
C_init = np.array([Xo, So, Po])  
t_lim = (0, 10)         #min
t_array = np.linspace(0, 10, num=100)

# Resolver el modelo
sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array,
        args=(UM, Kd, ks, Ypx, mp, Yxs, ms))

t_array = sol.t
conc_array = sol.y
Biomasa = conc_array[0]
Sustrato = conc_array[1]
Producto = conc_array[2]

#Graficar perfiles de concentración
figure1, ((ax1, ax2, ax3)) = plt.subplots(3, 1)
line1, = ax1.plot(t_array, Biomasa, label = "Biomasa", color = "darkorange")
ax1.legend()
ax1.set_xlabel(r'$tiempo (h)$')
ax1.set_ylabel(r'$[X] \/  (g \/ L ^{-1})$')
ax1.grid()


line2, = ax2.plot(t_array, Sustrato, label = "Sustrato", color = "darkblue")
ax2.legend()
ax2.set_xlabel(r'$tiempo (h)$')
ax2.set_ylabel(r'$[S] \/ (g \/ L ^{-1})$')
ax2.grid()

line3, = ax3.plot(t_array, Producto, label = "Producto", color = "red")
ax3.legend()
ax3.set_xlabel(r'$tiempo (h)$')
ax3.set_ylabel(r'$[P] \/ (g \/ L ^{-1})$')
ax3.grid()



# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.3)
axcolor = 'lightcoral'
ax_miu = plt.axes([0.05, 0.8, 0.15, 0.05], facecolor=axcolor)
ax_Kd = plt.axes([0.05, 0.7, 0.15, 0.05], facecolor=axcolor)
ax_ks = plt.axes([0.05, 0.6, 0.15, 0.05], facecolor=axcolor)
ax_Ypx = plt.axes([0.05, 0.5, 0.15, 0.05], facecolor=axcolor)
ax_mp = plt.axes([0.05, 0.4, 0.15, 0.05], facecolor=axcolor)
ax_Yxs = plt.axes([0.05, 0.3, 0.15, 0.05], facecolor=axcolor)
ax_ms = plt.axes([0.05, 0.2, 0.15, 0.05], facecolor=axcolor)

s_miu= Slider(ax = ax_miu, 
             label = r'$\mu$ ($h ^ {-1}$)',
             valmin = 0.1, 
             valmax = 1, 
             valinit=UM,
             valfmt='%1.1f')
s_Kd= Slider(ax = ax_Kd, label = r'$K_d$ ($ \/h ^ {-1}$)', valmin = 0.000001 , 
             valmax = 0.5, valinit=Kd, valfmt='%1.1f')

s_ks= Slider(ax = ax_ks, 
             label = r'$K_s$ ($g \/L ^ {-1}$)',
             valmin = 0.000001 , 
             valmax = 5, 
             valinit=ks,
             valfmt='%1.1f')

s_Ypx= Slider(ax = ax_Ypx, label = r'$Y_{PX}$ ($ gP\/gX ^ {-1}$)', valmin = 1 , 
             valmax = 10, valinit=Kd, valfmt='%1.1f')

s_mp= Slider(ax = ax_mp, label = r'$m_{p}$ ($ gP\/gX ^ {-1}$)', valmin = 0.5 , 
             valmax = 2, valinit=mp, valfmt='%1.1f')

s_Yxs= Slider(ax = ax_Yxs, label = r'$Y_{XS}$ ($ gX\/gS ^ {-1}$)', valmin = 0.001 , 
             valmax = 1, valinit=Yxs, valfmt='%1.1f')

s_ms= Slider(ax = ax_ms, label = r'$m_{S}$ ($ gS\/gS ^ {-1}$)', valmin = 0.5 , 
             valmax = 5, valinit=ms, valfmt='%1.1f')


def update(val):
    UM = s_miu.val
    Kd = s_Kd.val
    ks = s_ks.val
    Ypx = s_Ypx.val
    mp = s_mp.val
    Yxs = s_Yxs.val
    ms = s_ms.val

    sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array, 
        args=(UM, Kd, ks, Ypx, mp, Yxs, ms))
    conc_array = sol.y
    Biomasa = conc_array[0]
    Sustrato = conc_array[1]
    Producto = conc_array[2]
    
    line1.set_ydata(Biomasa)
    line2.set_ydata(Sustrato)
    line3.set_ydata(Producto)
    
    figure1.canvas.draw_idle()

# register the update function with each slider
s_miu.on_changed(update)
s_Kd.on_changed(update)
s_ks.on_changed(update)
s_Ypx.on_changed(update)
s_mp.on_changed(update)
s_Yxs.on_changed(update)
s_ms.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.1, 0.9, 0.09, 0.04])
button = Button(resetax, 'Reset variables', color='cornflowerblue', hovercolor='0.975')

def reset(event):
    s_miu.reset()
    s_Kd.reset()
    s_ks.reset()
    s_Ypx.reset()
    s_mp.reset()
    s_Yxs.rese()
    s_ms.reset()
button.on_clicked(reset)


#%% Cambiando las condiciones iniciales

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button



#Definir el modelo
def batch_monod(t, C):
    X, S, P = C  
    if S > 0 :
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X       #SUBSTRATE RATE EQUATION, kg/m3 h Producto asociado 
                                    #directamente con el metabolismo energetico
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = RS   
        dP_dt = RP   
    else:
        S = 0    
        U = UM*S/(ks+S)  #;MONOD EQUATION, 1/h
        RX = (U-Kd)*X    #BIOMASS RATE EQUATION, kg/m3 h 
        RS = -(U/Yxs + ms)*X      #SUBSTRATE RATE EQUATION, kg/m3 h
        RP = (mp+Ypx*U)*X   #PRODUCT RATE EQUATION, g/m3 h
        dX_dt = RX   
        dS_dt = 0  
        dP_dt = RP       
    return [dX_dt, dS_dt, dP_dt]


#Datos del modelo
UM=0.3     #1/h
Kd=0.0001      # 1/h
ks=0.000001    #g/L
Ypx=7.7     #gP/gX h
mp=1.1   #gP/gX h
Yxs=0.06    #gX/gS
ms= 2.2         #gX/gs
Xo=0.1      #g/L
So=12       #g/L
Po=0    #g/L


#Condiciones Iniciales
C_init = np.array([Xo, So, Po])  
t_lim = (0, 10)         #min
t_array = np.linspace(0, 10, num=100)

# Resolver el modelo
sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array)

t_array = sol.t
conc_array = sol.y
Biomasa = conc_array[0]
Sustrato = conc_array[1]
Producto = conc_array[2]

#Graficar perfiles de concentración
figure1, ((ax1, ax2, ax3)) = plt.subplots(3, 1)
line1, = ax1.plot(t_array, Biomasa, label = "Biomasa", color = "darkorange")
ax1.legend()
ax1.set_xlabel(r'$tiempo (h)$')
ax1.set_ylabel(r'$[X] \/  (g \/ L ^{-1})$')
ax1.grid()
ax1.set_ylim([0, 2])


line2, = ax2.plot(t_array, Sustrato, label = "Sustrato", color = "darkblue")
ax2.legend()
ax2.set_xlabel(r'$tiempo (h)$')
ax2.set_ylabel(r'$[S] \/ (g \/ L ^{-1})$')
ax2.grid()
ax2.set_ylim([0, 20])

line3, = ax3.plot(t_array, Producto, label = "Producto", color = "red")
ax3.legend()
ax3.set_xlabel(r'$tiempo (h)$')
ax3.set_ylabel(r'$[P] \/ (g \/ L ^{-1})$')
ax3.grid()
ax3.set_ylim([0, 20])



# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.3)
axcolor = 'azure'
ax_X0 = plt.axes([0.05, 0.7, 0.15, 0.05], facecolor=axcolor)
ax_S0 = plt.axes([0.05, 0.6, 0.15, 0.05], facecolor=axcolor)
ax_P0 = plt.axes([0.05, 0.5, 0.15, 0.05], facecolor=axcolor)


s_X0= Slider(ax = ax_X0, label = r'$X_{0}$ ($ g\/L ^ {-1}$)', valmin = 0.1 , 
             valmax = 1, valinit=Xo, valfmt='%1.1f')

s_S0= Slider(ax = ax_S0, label = r'$S_{0}$ ($ g\/L ^ {-1}$)', valmin = 1 , 
             valmax = 50, valinit=So, valfmt='%1.1f')

s_P0= Slider(ax = ax_P0, label = r'$P_{o}$ ($ g\/L ^ {-1}$)', valmin = 0 , 
             valmax = 1, valinit=Po, valfmt='%1.1f')


def update(val):
    Xo = s_X0.val
    So = s_S0.val
    Po= s_P0.val
    C_init = np.array([Xo, So, Po])

    sol = solve_ivp(batch_monod, 
        t_lim,
        C_init,
        t_eval=t_array)
    conc_array = sol.y
    Biomasa = conc_array[0]
    Sustrato = conc_array[1]
    Producto = conc_array[2]
    
    line1.set_ydata(Biomasa)
    line2.set_ydata(Sustrato)
    line3.set_ydata(Producto)
    
    figure1.canvas.draw_idle()

# register the update function with each slider
s_X0.on_changed(update)
s_S0.on_changed(update)
s_P0.on_changed(update)


# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.1, 0.9, 0.09, 0.04])
button = Button(resetax, 'Reset variables', color='cornflowerblue', hovercolor='0.975')

def reset(event):
    s_X0.reset()
    s_S0.reset()
    s_P0.reset()
   
button.on_clicked(reset)
