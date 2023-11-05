# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 16:19:23 2023

@author: BI-Yandong fron Tongji Univ

function: calculate the excess pore-pressure for a uniform and isotropic seabed based on the quasi-static Biot's consolidation equations (Biot, 1941), 

According to: Zhang, Y., Jeng, D.-S., Gao, F.P., Zhang, J.S., 2013. An analytical solution for response of a porous seabed to combined wave and currents loading. Ocean Eng. 57 (57), 240-247.
"""

import math
import waveLength0
import matplotlib.pylab as plt
import numpy as np


#----------------------seabed soil parameters-----------------------
g = 9.81                    # gravity acc (m/s^2)
rho_w = 1000                # density of water (kg/m^3)
gamma_w = rho_w*g           # unit weight of water 
rho_s = 2400                # density of soil (kg/m^3)
gamma_s = rho_s*g           # unit weight of soil

mu = 1.0e7                  # Shear modulus of soil G (Pa)
nu = 0.3                    # Poisson’s ratios of soil u
ks = 1.88e-4                # soil Permeability (m/s) 
n = 0.435                   # soil porosity n=e/(1+e)
beta = 4.5e-10              # compressibility of pore fluid

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
c0 = L0/T;                                            # wave velocity under wave-only
cc = 1 + 4*Uc/c0;                                     # wave velocity considering the effect of current

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
x  = float(input("Please enter position x: "))                      # x-coordinate
t = np.linspace(0, 5*T, 500)

Pb = P1*np.cos( k0*x - omega*t ) + P2*np.cos( 2*(k0*x - omega*t) ) + P3*np.cos( 3*(k0*x - omega*t) )
# plt.plot(t, Pb)
# plt.show()
print ('The max dynamic wave pressure on the seabed surface Pb_max =', np.max(Pb), 'Pa')

#---------3. calculate the excess pore-pressure within a seabed----------
z  = float(input("Please enter depth under seabed z: "))                      # z-coordinate, depth under seabed surface, negative value

alpha = ( 1-2*nu )*n*beta / ( n*beta + ( 1-2*nu )/mu )

delta = np.array([0 + 0j, 0 + 0j, 0 + 0j], dtype=complex)
C1    = np.array([0 + 0j, 0 + 0j, 0 + 0j], dtype=complex)
C2    = np.array([0 + 0j, 0 + 0j, 0 + 0j], dtype=complex)

P = np.array([P1, P2,P3])

EPP = 0 + 0j 

list = [1, 2, 3]
length = len(list)
for m in 1, 2, 3:
     
     delta[m-1] = np.sqrt( np.power(m,2) * np.power(k0,2) - (m*omega*rho_w*g/ks*( n*beta + (1-2*nu)/(2*mu*(1-nu)) ))*1j )
     C1[m-1] = ( delta[m-1] - delta[m-1]*nu + m*k0*nu ) / ( delta[m-1] - delta[m-1]*nu + m*k0*nu + m*k0*alpha )
     C2[m-1] = m*k0*alpha / ( ( delta[m-1] - m*k0 ) * ( delta[m-1] - delta[m-1]*nu + m*k0*nu + m*k0*alpha ) )
     
     EPP += P[m-1]/(1-2*nu) * ( ( 1-2*nu-alpha ) * C1[m-1]*np.exp(m*k0*z) + ( np.power(delta[m-1],2) - np.power(m,2)*np.power(k0,2) )/( m*k0 ) * ( 1-2*nu )*C2[m-1]*np.exp(delta[m-1]*z) ) * np.exp( (m*( k0*x - omega*t ))*1j )
     
plt.plot(t, EPP)
# plot initial vertical effective stress
y_position = (gamma_s - gamma_w)*z
plt.axhline(y_position, color='red', linestyle='--', label=r'$\sigma_{z0}^\prime$ (Pa)')

plt.show()

plt.title('H' + "%.3f"%H + 'T' + "%.3f"%T + 'h' + "%.3f"%h + "_EPP at z = " + "%.3f"%z + "m")

plt.xlabel('Time (s)')
plt.ylabel('EPP (Pa)')

plt.legend(loc='upper right')

file_name = 'H' + "%.3f"%H + 'T' + "%.3f"%T + 'h' + "%.3f"%h + "_EPP at z = " + "%.3f"%z
plt.savefig(f'{file_name}.png', dpi=600)

print (f'The peak excess pore pressure at z = {z} under seabed EPP_max =', np.max(EPP), 'Pa') 
     
     
