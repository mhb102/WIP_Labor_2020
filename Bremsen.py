#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Projekt-Neue Bremsen für die BSAG SOSE2020

Mario Hernandez     5004880
Lukas Hempe         5001989

"""


import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import scipy . integrate as sc


# Konstanten
c1 = -1.325
c0 = -0.6

v_0_kmh = 46.8
v_0_ms = round(46.8/3.6)

t = np.arange(-0.5,5.5,0.1) #festgelegtes Zeitarray

h = 0.1
# neues Bremssystem analytische Funktionen
a_t_n = c1*t+c0                     # Beschleunigung
v_t_n = -0.6625*t**2-0.6*t+13       # Geschwindigkeit
s_t_n = -53/240*t**3-3/10*t**2+13*t # Weg

# altes Bremssystem analytische Funktionen
a_t_a = -2.73                       # Beschleunigung
x = len(t)
a_array = np.ones(x)
a_array = a_array*-2.73           # Array zum Plotten der konstanten Beschleun.
v_t_a = -273/100 *t+13              # Geschwindigkeit
s_t_a = -273/200*t**2+13*t          # Weg


# Plotten der Funktionen 
plt.figure(1)
plt.plot(t,a_t_n,'k',label='neues Bremssystem')
plt.plot(t,a_array,'k''--',label='altes Bremssystem')
plt.xlabel('Zeit t / $s$')
plt.ylabel('Beschleunigung a(t) / $m/s²$')
plt.title('Beschleunigungs-Zeit-Diagramm Gefahrenbremsung')
plt.legend()
plt.grid(True) 
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.xticks([0,1,2,3,4,5,6])
plt.yticks([-9,-8,-7,-6,-5,-4,-3,-2,-1,0])


plt.figure(2)
plt.plot(t,v_t_n,'k',label='neues Bremssystem')
plt.plot(t,v_t_a,'k''--',label='altes Bremssystem')
plt.scatter(4,0,s=60,marker = 'x', color ='k')
plt.annotate(r'4 s',
         xy=(4,0), 
         xycoords='data',
         xytext=(3.3, -2), 
         #textcoords='offset points', 
         fontsize=8,
         arrowprops=dict(facecolor='k',
                         shrink=0.2,
                         headwidth=5, 
                         headlength=4, 
                         width=1))
plt.scatter(4.76,0,s=60,marker = 'x', color ='k')
plt.annotate(r'4,76 s',
         xy=(4.76,0), 
         xycoords='data',
         xytext=(4.76, 3.5), 
         #textcoords='offset points', 
         fontsize=8,
         arrowprops=dict(facecolor='k',
                         shrink=0.2,
                         headwidth=5, 
                         headlength=4,
                         width=1))
plt.xlabel('Zeit t / $s$')
plt.ylabel('Geschwindigkeit v(t) / $m/s$')
plt.title('Geschwindigkeits-Zeit-Diagramm Gefahrenbremsung')
plt.legend()
plt.grid(True)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.xticks([0,1,2,3,4,5,6])

plt.figure(3)
plt.plot(t,s_t_n,'k',label='neues Bremssystem')
plt.plot(t,s_t_a,'k''--',label='altes Bremssystem')
plt.scatter(4,33.07,s=60,marker = 'x', color ='k')
plt.annotate(r'33,07 m',
         xy=(4,33.07), 
         xycoords='data',
         xytext=(3.5, 26), 
         #textcoords='offset points', 
         fontsize=8,
         arrowprops=dict(facecolor='k',
                         shrink=0.2,
                         headwidth=5, 
                         headlength=4, 
                         width=1))
plt.scatter(4.76,30.95,s=60,marker = 'x', color ='k')
plt.annotate(r'30.95 m',
         xy=(4.76,30.95), 
         xycoords='data',
         xytext=(4.5, 26), 
         #textcoords='offset points', 
         fontsize=8,
         arrowprops=dict(facecolor='k',
                         shrink=0.2,
                         headwidth=5, 
                         headlength=4,
                         width=1))
plt.xlabel('Zeit t / $s$')
plt.ylabel('Strecke s(t) / $m$')
plt.title('Weg-Zeit-Diagramm Gefahrenbremsung')
plt.legend()
plt.grid(True)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.xticks([0,1,2,3,4,5,6])

#####Symbolische integration############# 

#bremsverzögerung
t = sym.symbols( 't' ) 
t0 = np.arange(-0.5,4.5,0.1)
v0 = 46.8 #km/h 
v0 = v0 /3600 *1000 #m/s 

#beschleunigung
a = c1*t+c0 
a_num =sym.lambdify(t,a)

#funktion
a_t=a_num(t0) 

#geschwindigkeit
v= sym. integrate (a)

#beschleunigung
v=v+v0
v_num =sym.lambdify(t,v)

#geschwindigkeit
v_t=v_num(t0)


#bremsweg
s=sym. integrate (v)

#geschwindigkeit
s = s+0
s_num = sym.lambdify(t,s)
s_t = s_num(t0) #zeit

t0_2 = np.arange(0 ,4.1 ,0.1)

s_t_2 = s_num(t0_2)
#zeitintervall

#funktion für v(t)
def v(x):
    return (-0.6625*x**2-0.6*x+13)

#mytrapz
def myTrapz(x,y):
    xi = np.arange(x,y+h,h)
    J = h*(0.5*v(xi[0]) + np.sum(v(xi[1:xi.size-1])) + 0.5*v(xi[-1])) # numerische integration formel
    return J
trapz=myTrapz(0 ,4)
#print('Trapez Integralberechnung:') 
#print(trapz)

#mycumtrapz
def myCumTrapz(x,y,initial):
    n = int((y-x)/h)
    xi = np.arange(x,y+h,h)
    F = np.zeros(n+1)
    F[0] = initial
    for i in np.arange(1, n+1): 
        F[i] = F[i-1] + h/2*(v(xi[i-1]) + v(xi[i])) 
    return F

CumTrapz=myCumTrapz ( 0 , 4 , 0 ) 
#print('Kumultative Integralberechnung:') 
#print (CumTrapz)

scipy_integral_1= sc.quad(v,0 ,4) #errechnet integralwet und fehler print ( ’ Scipy Quad Integralberechnung : ’ )
#print ( scipy_integral_1 )

def my_scipyint(x,y):
    F = np.zeros(len(y))
    for i in range(len(y)): 
        F[i] = sc.quad(v, x, y[i])[0] 
    return F

scipy_integral_2= my_scipyint(0,t0_2) 
#print( 'Scipy Quad Integralberechnung')
#print(scipy_integral_2)

#Differenzenberechnung
diff_1 = s_t_2 - CumTrapz
diff_2=s_t_2 - scipy_integral_2
diff_3=scipy_integral_2 - CumTrapz

# #plot bremsweg vergleich
plt.figure (4)
plt.plot(t0_2,s_t_2, 'ko',label='s(t) Analytisch')
plt.plot(t0_2,CumTrapz, 'k--',label='s(t) myCumTrapz')
plt.plot(t0_2,scipy_integral_2, 'r.',label='s(t) my_scipyint')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on ()
plt.grid(b=True, which='minor', color='#999999', linestyle='-',alpha =0.2)
plt.legend(loc='lower right')
plt.title('Vergleich der Lösungsansätze')
plt.xlabel('Zeit t in s')
plt.ylabel('Strecke s in m')



# Differenzen Plot
plt.figure()
plt.plot(t0_2,diff_1, 'kx', label='Differenz Ana & Cum')
plt.plot(t0_2,diff_3, 'k-', label='Differenz Scipy & Cum')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-',alpha =0.2)
plt.minorticks_on()
plt.legend(loc='lower right')
plt.title('Differenz Ana,Scipy und MyCumTrapz')
plt.xlabel('Zeit t in s')
plt.ylabel('Streckendifferenz in m')


plt.figure()
plt.stem(t0_2,diff_2,linefmt='white',markerfmt='ko',label='Differenz Ana & Scipy')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-',alpha =0.2)
plt.minorticks_on()
plt.legend(loc='upper left')
plt.title('Differenz Analytisch und Scipy')
plt.xlabel('Zeit t in s')
plt.ylabel('Streckendifferenz in m')

