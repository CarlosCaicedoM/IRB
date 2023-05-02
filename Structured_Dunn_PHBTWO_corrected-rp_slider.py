# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 23:21:42 2020

@author: Carlos Caicedo-Montoya

I. J. Dunn, E. Heinzle, J. Ingham, J. E. Pfenosil
Biological Reaction Engineering
Dynamic Modelling Fundamentals
with Simulation Examples
Second, Completely Revised Edition
2003 WILEY-VCH

Módelo estructurado para la producción de polihidroxibutiratos (PHB)


Heinzle y Lafferty (1980)  han presentado un modelo  para describir el cultivo 
por lotes del microorganismo Alcaligenes eutrophus bajo condiciones  
de crecimiento quimiolitoautotróficas.
El crecimiento y el almacenamiento de los PHB estan descritas como funciones 
del sustrato limitante (NH4+), la biomasa residual R y el producto (PHBs) 


En el modelo, las células completas (X) constan de dos compartimentos, 
a saber, PHB (P) y la biomasa residual (R), donde R se calcula como la
diferencia entre el peso seco celular total y la concentración de PHB (R =
X - P). R puede considerarse como la biomasa catalíticamente activa, que 
incluye proteínas y ácidos nucleicos. Con concentraciones constantes de los
gases disueltos, se pueden reconocer dos fases distintas: crecimiento y 
almacenamiento. 
Durante la fase de crecimiento hay suficiente NH4+ para permitir la 
síntesis de proteínas. Cuando el sustrato limitante NH4 + (S) se agota,
la síntesis de proteínas cesa y aumenta la velocidad de producción de PHB. 
Durante la fase de almacenamiento, solo PHB es producido.
El sustrato limitante NH4 + (S) es esencial para producir R y limita su 
síntesis a bajas concentraciones.

En este programa se implementa dicho ejemplo para dos reactores en 
continuo en sería

"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

#Importar librerías
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#Definir Modelo
def produccion_de_PHB_continuo(t, C):
    R1, S1, P1, R2, S2, P2= C
    #Para el reactor 1
    miu_11 = umax_1*S1/(Ks1 + S1)
    miu_21 = umax_2*((S1/Ks2)**n)/(1+(S1/Ks2)**n)
    miu1 = miu_11 + miu_21
    
    dR1 = miu1 * R1 - D1*R1
    dS1 = (-1/Yrs)*miu1*R1 + D1*(Si-S1)
    dP1 = Ypr*miu1*R1 + (KI/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
    #dP1 = Ypr*miu1*R1 + (KI*R1/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
    
    #Para el reactor 2
    miu_12 = umax_1*S2/(Ks1 + S2)
    miu_22 = umax_2*((S2/Ks2)**n)/(1+(S2/Ks2)**n)
    miu2 = miu_12 + miu_22
    dR2 = miu2 * R2 + D2*(R1-R2)
    dS2 = (-1/Yrs)*miu2*R2 + D2*(S1-S2)
    dP2 = Ypr*miu2*R2 + (KI/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
    #dP = Ypr*miu*R + (KI/(KI+S))*(-k1*P+k2*R)
    #dP2 = Ypr*miu2*R2 + (KI*R2/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)

    return np.array([dR1, dS1, dP1, dR2, dS2, dP2])

#Datos
Yrs = 1.5
umax_1 = 0.13   #h{^-1}
Ks1 = 0.1     #g/L
umax_2 = 0.08     #h^-1
Ks2 = 1.0       #g/L
n = 4
Ypr = 0.105
KI = 0.036      #g/L   0.041
k1 = 0.045      #g/L
k2 = 0.18       #g/L
F = 1.0   #L/H 

VTOT=20  #Total volume
VRAT=0.5     #{V1/Vtotal volume ratio. Must be less than 1}
V1 = VTOT*VRAT     #L
V2 = VTOT-V1     #L
D1 = F/V1
D2 =  F/V2
Si= 5

#Condiciones iniciales y parámetros del método númerico
R0 = 0.22    #g/L
S0 = 2.3  #g/L 
P0 = 0.22   #g/L


C_init = np.array([R0, S0, P0, R0, S0, P0])
t_span = (0, 100)
t_array = np.linspace(0, 100, 1000)
sol = solve_ivp(produccion_de_PHB_continuo, t_span, C_init, t_eval=t_array)

conc_array = sol.y
Biomasa_residual_1 = conc_array[0]
Sustrato_1 = conc_array[1]
Producto_1 = conc_array[2]
Biomasa_total_1 =  conc_array[2] + conc_array[0]

Biomasa_residual_2 = conc_array[3]
Sustrato_2 = conc_array[4]
Producto_2 = conc_array[5]
Biomasa_total_2 =  conc_array[3] + conc_array[5]

#Gráfica de las concentraciones de los componentes del modelo

#Gráfica de las concentraciones de los componentes del modelo
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
line1, = ax1.plot(t_array, Biomasa_residual_1, label = 'Biomasa residual 1')
line2, = ax1.plot(t_array, Biomasa_residual_2, label= 'Biomasa residual 2')
ax1.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax1.legend()
ax1.grid()


line3, = ax2.plot(t_array, Sustrato_1, color = 'gold', label = 'Sustrato 1')
line4, = ax2.plot(t_array, Sustrato_2, color='darkred', label = 'Sustrato 2')
ax2.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax2.grid()
ax2.legend()

line5, = ax3.plot(t_array, Producto_1, color = 'green', label = 'Producto 1')
line6, = ax3.plot(t_array, Producto_2, color = 'red', label = 'Producto 2')
ax3.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax3.grid()
ax3.legend()

line7, = ax4.plot(t_array, Biomasa_total_1, color = 'steelblue', label = 'Biomsa total 1')
line8, = ax4.plot(t_array, Biomasa_total_2, color= 'orchid', label = 'Biomasa total 2')
ax4.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax4.grid()
ax4.legend()

#axes_limits
ax1.set_ylim([0, max(Biomasa_residual_2)+0.2])
ax2.set_ylim([0, max(Sustrato_1)+0.2])
ax3.set_ylim([0, max(Producto_2)+0.2])
ax4.set_ylim([0, max(Biomasa_total_2)+0.2])


#%%  Investigate the time to reach steady state for a range of flowrate values 
from IPython import get_ipython
get_ipython().magic('reset -sf')

#Importar librerías
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

from scipy.integrate import solve_ivp

#Definir Modelo
def produccion_de_PHB_continuo(t, C, F):
    R1, S1, P1, R2, S2, P2= C
    
    D1 = F/V1
    D2 =  F/V2
    
    #Para el reactor 1
    miu_11 = umax_1*S1/(Ks1 + S1)
    miu_21 = umax_2*((S1/Ks2)**n)/(1+(S1/Ks2)**n)
    miu1 = miu_11 + miu_21
    dR1 = miu1 * R1 - D1*R1
    dS1 = (-1/Yrs)*miu1*R1 + D1*(Si-S1)
    dP1 = Ypr*miu1*R1 + (KI/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
    #dP1 = Ypr*miu1*R1 + (KI*R1/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
    
    #Para el reactor 2
    miu_12 = umax_2*S2/(Ks1 + S2)
    miu_22 = umax_2*((S2/Ks2)**n)/(1+(S2/Ks2)**n)
    miu2 = miu_12 + miu_22
    dR2 = miu2 * R2 + D2*(R1-R2)
    dS2 = (-1/Yrs)*miu2*R2 + D2*(S1-S2)
    dP2 = Ypr*miu2*R2 + (KI/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
    #dP2 = Ypr*miu2*R2 + (KI*R2/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)

    return np.array([dR1, dS1, dP1, dR2, dS2, dP2])

#Datos
Yrs = 1.5
umax_1 = 0.13   #h{^-1}
Ks1 = 0.1     #g/L
umax_2 = 0.08     #h^-1
Ks2 = 1.0       #g/L
n = 4
Ypr = 0.105
KI = 0.036      #g/L   0.041
k1 = 0.045      #g/L
k2 = 0.18       #g/L
F1 = 1.0   #L/H 
VTOT=20  #Total volume
VRAT=0.5     #{V1/Vtotal volume ratio. Must be less than 1}
V1 = VTOT*VRAT     #L
V2 = VTOT-V1     #L

Si= 5

#Condiciones iniciales y parámetros del método númerico
R0 = 0.22    #g/L
S0 = 2.3  #g/L 
P0 = 0.22   #g/L


C_init = np.array([R0, S0, P0, R0, S0, P0])
t_span = (0, 100)
t_array = np.linspace(0, 100, 1000)
sol = solve_ivp(produccion_de_PHB_continuo, t_span, C_init,
                t_eval=t_array, args=(F1,))

conc_array = sol.y
Biomasa_residual_1 = conc_array[0]
Sustrato_1 = conc_array[1]
Producto_1 = conc_array[2]
Biomasa_total_1 =  conc_array[2] + conc_array[0]

Biomasa_residual_2 = conc_array[3]
Sustrato_2 = conc_array[4]
Producto_2 = conc_array[5]
Biomasa_total_2 =  conc_array[3] + conc_array[5]

#Gráfica de las concentraciones de los componentes del modelo
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
line1, = ax1.plot(t_array, Biomasa_residual_1, label = 'Biomasa residual 1')
line2, = ax1.plot(t_array, Biomasa_residual_2, label= 'Biomasa residual 2')
ax1.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax1.legend()
ax1.grid()


line3, = ax2.plot(t_array, Sustrato_1, color = 'gold', label = 'Sustrato 1')
line4, = ax2.plot(t_array, Sustrato_2, color='darkred', label = 'Sustrato 2')
ax2.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax2.grid()
ax2.legend()

line5, = ax3.plot(t_array, Producto_1, color = 'green', label = 'Producto 1')
line6, = ax3.plot(t_array, Producto_2, color = 'red', label = 'Producto 2')
ax3.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax3.grid()
ax3.legend()

line7, = ax4.plot(t_array, Biomasa_total_1, color = 'steelblue', label = 'Biomsa total 1')
line8, = ax4.plot(t_array, Biomasa_total_2, color= 'orchid', label = 'Biomasa total 2')
ax4.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax4.grid()
ax4.legend()

#axes_limits
ax1.set_ylim([0, max(Biomasa_residual_2)+0.2])
ax2.set_ylim([0, max(Sustrato_1)+0.2])
ax3.set_ylim([0, max(Producto_2)+0.2])
ax4.set_ylim([0, max(Biomasa_total_2)+0.2])


#Fig_Title
fig1.suptitle(r'$F_0={0:0.2f}$'.format(F1))


# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.3)
axcolor = 'lightcoral'
ax_F = plt.axes([0.05, 0.5, 0.15, 0.05], facecolor=axcolor)

s_F= Slider(ax = ax_F, 
             label = r'$F_0$ ($ m^3  h ^ {-1}$)',
             valmin = 0.1, 
             valmax = 2, 
             valinit=F1,
             valfmt='%1.1f')


def update(val):
    F1= s_F.val
    sol = solve_ivp(produccion_de_PHB_continuo, t_span, C_init,
                t_eval=t_array, args=(F1,))
    conc_array = sol.y
    Biomasa_residual_1 = conc_array[0]
    Sustrato_1 = conc_array[1]
    Producto_1 = conc_array[2]
    Biomasa_total_1 =  conc_array[2] + conc_array[0]

    Biomasa_residual_2 = conc_array[3]
    Sustrato_2 = conc_array[4]
    Producto_2 = conc_array[5]
    Biomasa_total_2 =  conc_array[3] + conc_array[5]
    
    line1.set_ydata(Biomasa_residual_1)
    line2.set_ydata(Biomasa_residual_2)
    line3.set_ydata(Sustrato_1)
    line4.set_ydata(Sustrato_2)
    line5.set_ydata(Producto_1)
    line6.set_ydata(Producto_2)
    line7.set_ydata(Biomasa_total_1)
    line8.set_ydata(Biomasa_total_2)
    
    #axes_limits
    ax1.set_ylim([0, max(Biomasa_residual_2)+0.2])
    ax2.set_ylim([0, max(Sustrato_1)+0.2])
    ax3.set_ylim([0, max(Producto_2)+0.2])
    ax4.set_ylim([0, max(Biomasa_total_2)+0.2])


    #Fig_Title
    fig1.suptitle(r'$F_0={0:0.2f}$'.format(F1))
    
    
    fig1.canvas.draw_idle()

# register the update function with each slider
s_F.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = plt.axes([0.1, 0.84, 0.09, 0.04])
button = Button(resetax, 'Reset variables', color='cornflowerblue', hovercolor='0.975')

def reset(event):
    s_F.reset()
button.on_clicked(reset)

#%% 3 
#keeping VRAT = 0.5 make a parameter plot to see the influence of 
# the flow rate on the productivities 

from IPython import get_ipython
get_ipython().magic('reset -sf')

#Importar librerías
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.integrate import solve_ivp

F1 = list(np.linspace(0.1, 2.0))
Prod1 = []
ProdFinal = []

for i in F1: 
    #Definir Modelo
    def produccion_de_PHB_continuo(t, C, F):
        R1, S1, P1, R2, S2, P2= C
        
        
    
      
        #Para el reactor 1
        miu_11 = umax_1*S1/(Ks1 + S1)
        miu_21 = umax_2*((S1/Ks2)**n)/(1+(S1/Ks2)**n)
        miu1 = miu_11 + miu_21
        dR1 = miu1 * R1 - D1*R1
        dS1 = (-1/Yrs)*miu1*R1 + D1*(Si-S1)
        dP1 = Ypr*miu1*R1 + (KI/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
        #dP1 = Ypr*miu1*R1 + (KI*R1/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
        
        #Para el reactor 2
        miu_12 = umax_2*S2/(Ks1 + S2)
        miu_22 = umax_2*((S2/Ks2)**n)/(1+(S2/Ks2)**n)
        miu2 = miu_12 + miu_22
        dR2 = miu2 * R2 + D2*(R1-R2)
        dS2 = (-1/Yrs)*miu2*R2 + D2*(S1-S2)
        dP2 = Ypr*miu2*R2 + (KI/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
        #dP2 = Ypr*miu2*R2 + (KI*R2/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
    
        return np.array([dR1, dS1, dP1, dR2, dS2, dP2])
    
    #Datos
    Yrs = 1.5
    umax_1 = 0.13   #h{^-1}
    Ks1 = 0.1     #g/L
    umax_2 = 0.08     #h^-1
    Ks2 = 1.0       #g/L
    n = 4
    Ypr = 0.105
    KI = 0.036      #g/L   0.041
    k1 = 0.045      #g/L
    k2 = 0.18       #g/L
    F = i   #L/H 
    
    Si= 5
    
    VTOT=20  #Total volume
    VRAT=0.5     #{V1/Vtotal volume ratio. Must be less than 1}
    V1 = VTOT*VRAT     #L
    V2 = VTOT-V1     #L
    D1 = F/V1
    D2 =  F/V2
    #Condiciones iniciales y parámetros del método númerico
    R0 = 0.22    #g/L
    S0 = 2.3  #g/L 
    P0 = 0.22   #g/L
    
    
    C_init = np.array([R0, S0, P0, R0, S0, P0])
    t_span = (0, 100)
    t_array = np.linspace(0, 100, 1000)
    sol = solve_ivp(produccion_de_PHB_continuo, t_span, C_init,
                    t_eval=t_array, args=(F,))

    conc_array = sol.y
    Biomasa_residual_1 = conc_array[0]
    Sustrato_1 = conc_array[1]
    Producto_1 = conc_array[2]
    Biomasa_total_1 =  conc_array[2] + conc_array[0]
    
    Biomasa_residual_2 = conc_array[3]
    Sustrato_2 = conc_array[4]
    Producto_2 = conc_array[5]
    Biomasa_total_2 =  conc_array[3] + conc_array[5]
    
   
    
    
    PROD=(Producto_2/VTOT)[-1]*i      #kg product/m3 h}
    PROD1=(Producto_1/V1)[-1]*i  #{Productivity from the first stage., kg product/m3 h}
    
    
    Prod1.append(PROD1) 
    ProdFinal.append(PROD)





#Plot concentration profiles

ProdFinal = np.asarray(ProdFinal)
Prod1 = np.asarray(Prod1)
F1 = np.asarray(F1)


figure2, axis2 = plt.subplots()
axis2.plot(F1, ProdFinal, '--', label = "Productivity final")
axis2.plot(F1, Prod1, '+', label = "Productivity tank 1" )
axis2.legend()
axis2.grid()
axis2.set_xlabel(r'F$(m^3 h^{-1}$)')
axis2.set_ylabel('Prod (Final) (kg m$^{-3}$h$^{-1}$) ')


print("Productividad máxima se da a un flujo volumétrico de {0:0.2f} m3/h".format(np.amax(ProdFinal)))

#%% 4

#Using Flow rate value at maximun in PROD make a parameter plot to see
#the influence of the volume rario, VRAT, on the productivities
#PROD1 and PROD 


#From the previous exercise 
#F0=1.23 m3/h

from IPython import get_ipython
get_ipython().magic('reset -sf')

#Importar librerías
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.integrate import solve_ivp

VRAT1 = list(np.linspace(0.3, 0.9))
Prod1 = []
ProdFinal = []

for i in VRAT1: 
    #Definir Modelo
    def produccion_de_PHB_continuo(t, C, VRAT):
        R1, S1, P1, R2, S2, P2= C
        
        
    
      
        #Para el reactor 1
        miu_11 = umax_1*S1/(Ks1 + S1)
        miu_21 = umax_2*((S1/Ks2)**n)/(1+(S1/Ks2)**n)
        miu1 = miu_11 + miu_21
        dR1 = miu1 * R1 - D1*R1
        dS1 = (-1/Yrs)*miu1*R1 + D1*(Si-S1)
        dP1 = Ypr*miu1*R1 + (KI/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
        #dP1 = Ypr*miu1*R1 + (KI*R1/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
        
        #Para el reactor 2
        miu_12 = umax_2*S2/(Ks1 + S2)
        miu_22 = umax_2*((S2/Ks2)**n)/(1+(S2/Ks2)**n)
        miu2 = miu_12 + miu_22
        dR2 = miu2 * R2 + D2*(R1-R2)
        dS2 = (-1/Yrs)*miu2*R2 + D2*(S1-S2)
        dP2 = Ypr*miu2*R2 + (KI/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
        #dP2 = Ypr*miu2*R2 + (KI*R2/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
    
        return np.array([dR1, dS1, dP1, dR2, dS2, dP2])
    
    #Datos
    Yrs = 1.5
    umax_1 = 0.13   #h{^-1}
    Ks1 = 0.1     #g/L
    umax_2 = 0.08     #h^-1
    Ks2 = 1.0       #g/L
    n = 4
    Ypr = 0.105
    KI = 0.036      #g/L   0.041
    k1 = 0.045      #g/L
    k2 = 0.18       #g/L
    F = 1.23   #L/H   #From the previous exercise
     
    Si= 5
    
    VTOT=20  #Total volume
    VRAT=i     #{V1/Vtotal volume ratio. Must be less than 1}
    V1 = VTOT*VRAT     #L
    V2 = VTOT-V1     #L
    D1 = F/V1
    D2 =  F/V2
    #Condiciones iniciales y parámetros del método númerico
    R0 = 0.22    #g/L
    S0 = 2.3  #g/L 
    P0 = 0.22   #g/L
    
    
    C_init = np.array([R0, S0, P0, R0, S0, P0])
    t_span = (0, 100)
    t_array = np.linspace(0, 100, 1000)
    sol = solve_ivp(produccion_de_PHB_continuo, t_span, C_init,
                    t_eval=t_array, args=(VRAT,))

    conc_array = sol.y
    Biomasa_residual_1 = conc_array[0]
    Sustrato_1 = conc_array[1]
    Producto_1 = conc_array[2]
    Biomasa_total_1 =  conc_array[2] + conc_array[0]
    
    Biomasa_residual_2 = conc_array[3]
    Sustrato_2 = conc_array[4]
    Producto_2 = conc_array[5]
    Biomasa_total_2 =  conc_array[3] + conc_array[5]
    
   
    
    
    PROD=(Producto_2/VTOT)[-1]*F      #kg product/m3 h}
    PROD1=(Producto_1/V1)[-1]*F  #{Productivity from the first stage., kg product/m3 h}
    
    
    Prod1.append(PROD1) 
    ProdFinal.append(PROD)





#Plot concentration profiles

ProdFinal = np.asarray(ProdFinal)
Prod1 = np.asarray(Prod1)
VRAT1 = np.asarray(VRAT1)


figure2, axis2 = plt.subplots()
axis2.plot(VRAT1, ProdFinal, '--', label = "Productivity final")
axis2.plot(VRAT1, Prod1, '+', label = "Productivity tank 1" )
axis2.legend()
axis2.grid()
axis2.set_xlabel(r'VRAT$(m^3 h^{-1}$)')
axis2.set_ylabel('Prod (Final) (kg m$^{-3}$h$^{-1}$) ')

#%% 5 Find the values of VRAT and F that maximises the productivity
from IPython import get_ipython
get_ipython().magic('reset -sf')

#Importar librerías
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.optimize import shgo
from scipy.integrate import solve_ivp


def objective (x):
    #Definir Modelo
    
    def produccion_de_PHB_continuo(t, C):
        R1, S1, P1, R2, S2, P2= C
       
        #Para el reactor 1
        miu_11 = umax_1*S1/(Ks1 + S1)
        miu_21 = umax_2*((S1/Ks2)**n)/(1+(S1/Ks2)**n)
        miu1 = miu_11 + miu_21
        dR1 = miu1 * R1 - D1*R1
        dS1 = (-1/Yrs)*miu1*R1 + D1*(Si-S1)
        dP1 = Ypr*miu1*R1 + (KI/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
        #dP1 = Ypr*miu1*R1 + (KI*R1/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
        
        #Para el reactor 2
        miu_12 = umax_2*S2/(Ks1 + S2)
        miu_22 = umax_2*((S2/Ks2)**n)/(1+(S2/Ks2)**n)
        miu2 = miu_12 + miu_22
        dR2 = miu2 * R2 + D2*(R1-R2)
        dS2 = (-1/Yrs)*miu2*R2 + D2*(S1-S2)
        dP2 = Ypr*miu2*R2 + (KI/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
        #dP2 = Ypr*miu2*R2 + (KI*R2/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
    
        return np.array([dR1, dS1, dP1, dR2, dS2, dP2])
    
    #Datos
    Yrs = 1.5
    umax_1 = 0.13   #h{^-1}
    Ks1 = 0.1     #g/L
    umax_2 = 0.08     #h^-1
    Ks2 = 1.0       #g/L
    n = 4
    Ypr = 0.105
    KI = 0.036      #g/L   0.041
    k1 = 0.045      #g/L
    k2 = 0.18       #g/L
    Si= 5
    F=x[0]
    VRAT=x[1]
    VTOT=20  #Total volume
     
    V1 = VTOT*VRAT     #L
    V2 = VTOT-V1     #L
    D1 = F/V1
    D2 =  F/V2
    #Condiciones iniciales y parámetros del método númerico
    R0 = 0.22    #g/L
    S0 = 2.3  #g/L 
    P0 = 0.22   #g/L
    
    
    C_init = np.array([R0, S0, P0, R0, S0, P0])
    t_span = (0, 100)
    t_array = np.linspace(0, 100, 1000)
    sol = solve_ivp(produccion_de_PHB_continuo, t_span, C_init,
                    t_eval=t_array)
    
    conc_array = sol.y

    Producto_2 = conc_array[5]
    PROD=(Producto_2/VTOT)[-1]*x[0]      #kg product/m3 h}
    return -PROD



result = shgo(objective, bounds = [(0.1, 2.0), (0.3, 0.9)])

print(result.x)

#Now, using the results from the previous optimization, the concentration profile 
#for the components of the model, the productivity and the quantity of product per:
#kilogram of biomass are plotted

def produccion_de_PHB_continuo(t, C):
    R1, S1, P1, R2, S2, P2= C
    #Para el reactor 1
    miu_11 = umax_1*S1/(Ks1 + S1)
    miu_21 = umax_2*((S1/Ks2)**n)/(1+(S1/Ks2)**n)
    miu1 = miu_11 + miu_21
    dR1 = miu1 * R1 - D1*R1
    dS1 = (-1/Yrs)*miu1*R1 + D1*(Si-S1)
    dP1 = Ypr*miu1*R1 + (KI/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
    #dP1 = Ypr*miu1*R1 + (KI*R1/(KI+S1))*(-k1*P1+k2*R1) - D1*P1
    
    #Para el reactor 2
    miu_12 = umax_2*S2/(Ks1 + S2)
    miu_22 = umax_2*((S2/Ks2)**n)/(1+(S2/Ks2)**n)
    miu2 = miu_12 + miu_22
    dR2 = miu2 * R2 + D2*(R1-R2)
    dS2 = (-1/Yrs)*miu2*R2 + D2*(S1-S2)
    dP2 = Ypr*miu2*R2 + (KI/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)
    #dP2 = Ypr*miu2*R2 + (KI*R2/(KI+S2))*(-k1*P2+k2*R2) + D2*(P1-P2)

    return np.array([dR1, dS1, dP1, dR2, dS2, dP2])

#Datos
Yrs = 1.5
umax_1 = 0.13   #h{^-1}
Ks1 = 0.1     #g/L
umax_2 = 0.08     #h^-1
Ks2 = 1.0       #g/L
n = 4
Ypr = 0.105
KI = 0.036      #g/L   0.041
k1 = 0.045      #g/L
k2 = 0.18       #g/L
F = result.x[0]   #L/H 

VTOT=20  #Total volume
VRAT= result.x[1]     #{V1/Vtotal volume ratio. Must be less than 1}
V1 = VTOT*VRAT     #L
V2 = VTOT-V1     #L
D1 = F/V1
D2 =  F/V2
Si= 5

#Condiciones iniciales y parámetros del método númerico
R0 = 0.22    #g/L
S0 = 2.3  #g/L 
P0 = 0.22   #g/L


C_init = np.array([R0, S0, P0, R0, S0, P0])
t_span = (0, 100)
t_array = np.linspace(0, 100, 1000)
sol = solve_ivp(produccion_de_PHB_continuo, t_span, C_init, t_eval=t_array)

conc_array = sol.y
Biomasa_residual_1 = conc_array[0]
Sustrato_1 = conc_array[1]
Producto_1 = conc_array[2]
Biomasa_total_1 =  conc_array[2] + conc_array[0]

Biomasa_residual_2 = conc_array[3]
Sustrato_2 = conc_array[4]
Producto_2 = conc_array[5]
Biomasa_total_2 =  conc_array[3] + conc_array[5]
PROD=(Producto_2/VTOT)*F


#Gráfica de las concentraciones de los componentes del modelo
fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
line1, = ax1.plot(t_array, Biomasa_residual_1, label = 'Biomasa residual 1')
line2, = ax1.plot(t_array, Biomasa_residual_2, label= 'Biomasa residual 2')
ax1.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax1.legend()
ax1.grid()


line3, = ax2.plot(t_array, Sustrato_1, color = 'gold', label = 'Sustrato 1')
line4, = ax2.plot(t_array, Sustrato_2, color='darkred', label = 'Sustrato 2')
ax2.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax2.grid()
ax2.legend()

line5, = ax3.plot(t_array, Producto_1, color = 'green', label = 'Producto 1')
line6, = ax3.plot(t_array, Producto_2, color = 'red', label = 'Producto 2')
ax3.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax3.grid()
ax3.legend()

line7, = ax4.plot(t_array, Biomasa_total_1, color = 'steelblue', label = 'Biomsa total 1')
line8, = ax4.plot(t_array, Biomasa_total_2, color= 'orchid', label = 'Biomasa total 2')
ax4.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1})$')
ax4.grid()
ax4.legend()



#axes_limits
ax1.set_ylim([0, max(Biomasa_residual_2)+0.2])
ax2.set_ylim([0, max(Sustrato_1)+0.2])
ax3.set_ylim([0, max(Producto_2)+0.2])
ax4.set_ylim([0, max(Biomasa_total_2)+0.2])

fig2, ((ax5, ax6))=plt.subplots(2,1)
line9, = ax5.plot(t_array, PROD, color = 'teal', label = 'Productivity')
line10, = ax6.plot(t_array, Producto_2/Biomasa_total_2, color= 'orchid', 
                   label = 'Producto por kg de Biomasa')
ax5.set(xlabel = 'tiempo ($hr$)', ylabel = r'Concentración  $(g \/ L^{-1}\/h^{-1})$')
ax5.grid()
ax5.legend()
ax6.set(xlabel = 'tiempo ($hr$)', ylabel = r' $(Kg Producto/Kg Biomasa)$')
ax6.grid()
ax6.legend()