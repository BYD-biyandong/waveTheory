# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 16:19:23 2023

@author: BI-Yandong fron Tongji Univ

function: calculate the dynamic wave pressure pd on seabed surface
According to: Ye, J.H., Jeng, D.-S., 2012. Response of porous seabed to nature loadings: waves and currents. Journal of Engineering Mechanics, ASCE 138 (6), 601-613.
"""

import math
import waveLength0
import matplotlib.pylab as plt
import numpy as np


g = 9.81                    # gravity acc (m/s^2)
rho_w = 1000                # density of water (kg/m^3)
gamma_w = rho_w*g           # unit weight of water 

#-------------------------wave parameters--------------------------
H0 = float(input("Please enter wave height under wave-only H0: "))     # wave height under wave-only
T  = float(input("Please enter wave period T: "))                      # wave period under wave-only
h  = float(input("Please enter water depth h: "))                      # water depth 

#-------------------------current parameters-----------------------------
Uc = float(input("Please enter current velocity without waves Uc: "))     # average velocity of the current without waves

#--------------calculate the wave length L0 under wave-only------------
# L0 = g*math.pow(T0,2)/(2*math.pi)* math.tanh(k0*h)        
L0 = waveLength0.L0(T, h)                             # wave length under wave-only
k0 = 2 * math.pi/L0                                   # wave number under wave-only

print ('Initial wave height L0 =', L0,'m')

#---------1. calculate the wave height and wave length considering the effect of current------
c0 = L0/T;                  # wave velocity under wave-only
cc = 1 + 4*Uc/c0;           # wave velocity considering the effect of current

# wave height H considering the effect of current
H = H0 * 2 * math.pow((cc + pow(cc,0.5)) , -0.5) * math.pow((1 + pow(cc,0.5)) , -0.5)

# wave length L considering the effect of current
L = L0 * (1/4.0) * math.pow((1 + pow(cc,0.5)) , 2)
k = 2 * math.pi/L                                    # wave number with current

#---------2. calculate the dynamic wave pressure on seabed surface----------

omega_0 = Uc * k0 + math.sqrt(g*k0* math.tanh(k0*h))  # Wave frequency under wave-only
omega_2 = ( 9 + 8*math.pow(math.sinh(k0*h),2) + 8*math.pow(math.sinh(k0*h),4) ) / ( 64*math.pow(math.sinh(k0*h),4) ) * (omega_0 - Uc*k0)
omega = omega_0 + math.pow(k0*H0,2) * omega_2  # Wave frequency

P1 = rho_w*g*H0 / ( 2*math.cosh(k0*h) ) * ( 1 - ( omega_2*k0*k0*H0*H0 )/( 2*Uc*k0 - omega_0 ) )
P2 = 3.0/8*rho_w*H0*H0 * ( ( omega_0*(omega_0-Uc*H0) )/( 2*math.pow(math.sinh(k0*h),4) ) - ( g*k0 )/( 3*math.sinh(2*k0*h) ) )
P3 = 3.0/512*rho_w*k*math.pow(H0,3)*omega_0*(omega_0 - Uc*k0) * ( 9 - 4*math.pow(math.sinh(k0*h),2) ) / ( math.pow(math.sinh(k0*h),7) )

# dynamic wave pressure acting on the seabed surface
x  = float(input("Please enter position x: "))                      # point position
t = np.linspace(0, 5*T, 500)
Pb = P1*np.cos( k0*x - omega*t ) + P2*np.cos( 2*(k0*x - omega*t) ) + P3*np.cos( 3*(k0*x - omega*t) )
plt.plot(t, Pb)
plt.show()

plt.title('H' + "%.3f"%H + 'T' + "%.3f"%T + 'h' + "%.3f"%h + "_Pb at x = " + "%.3f"%x)
plt.xlabel('Time (s)')
plt.ylabel('Wave pressure on seabed surface (Pa)')

file_name = 'H' + "%.3f"%H + 'T' + "%.3f"%T + 'h' + "%.3f"%h + "_Pb at x = " + "%.3f"%x
plt.savefig(f'{file_name}.png', dpi=600)

print (f'The max dynamic wave pressure at x = {x} on the seabed surface Pb_max =', np.max(Pb), 'Pa')


