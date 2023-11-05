# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 16:19:23 2023

@author: BI-Yandong fron Tongji Univ

function: calculate the wave height and wave length considering the effect of current
"""

import math
import waveLength0

import os 
import matplotlib.pylab as plt
import numpy as np
from PIL import Image

g = 9.81                     # gravity acc (m/s^2)
rho_w = 1000                 # density of water (kg/m^3)
gamma_w = rho_w*g            # unit weight of water 

#-------------------------wave parameters--------------------------
H0 = float(input("Please enter wave height under wave-only H0: "))     # wave height under wave-only
T  = float(input("Please enter wave period T: "))                      # wave period under wave-only
h  = float(input("Please enter water depth h: "))                      # water depth 

#-------------------------current parameters-----------------------------
Uc = float(input("Please enter current velocity without waves Uc: "))     # average velocity of the current without waves

#--------------calculate the wave length L0 under wave-only------------
# L0 = g*math.pow(T0,2)/(2*math.pi)* math.tanh(k0*h)        
L0 = waveLength0.L0(T, h) # wave length under wave-only

print ('Initial wave height L0 =', L0,'m')

#---------calculate the wave height and wave length considering the effect of current------
c0 = L0/T;                  # wave velocity under wave-only
cc = 1 + 4*Uc/c0;           # wave velocity considering the effect of current

# wave height H considering the effect of current
H = H0 * 2 * math.pow((cc + pow(cc,0.5)) , -0.5) * math.pow((1 + pow(cc,0.5)) , -0.5)

# wave length L considering the effect of current
L = L0 * (1/4.0) * math.pow((1 + pow(cc,0.5)) , 2)

print ('wave height H =', H,'m')
print ('wave length L =', L,'m')

#--------------plot wave diagram by Le Méhauté (1976) ------------
a = h/(g * math.pow(T,2.0))
b = H0/(g * math.pow(T,2.0))

xMin = 5.34*math.pow(10,-4.0)
xMinPix = 402
xMax = 0.2
xMaxPix = 1977

yMin = 0.00005
yMinPix = 1865
yMax = 0.05
yMaxPix = 28

if a > xMax:
    a = xMax
elif a < xMin:
    a = xMin

if b > yMax:
    b = yMax
elif b < yMin:
    b = yMin

auxX = 612*math.log10(xMax/a)
auxY = 612*math.log10(b/yMin)

coordX = xMaxPix - auxX
coordY = yMinPix - auxY

pathname = os.path.abspath('.')
image = Image.open('teoria.png')
arr = np.asarray(image)
plt.imshow(arr)

plt.plot(coordX,coordY,'ro',markersize=12)
plt.plot(coordX,coordY,'k.')
plt.plot(coordX,coordY,'kx',markersize=12)


plt.title('H = '+"%.3f"%H+' m; T = '+"%.3f"%T+' s; h = '+"%.3f"%h+' m.');

plt.xlim(-30,2030)
plt.ylim(2200,-30)

plt.show()

file_name = 'H' + "%.3f"%H + 'T' + "%.3f"%T + 'h' + "%.3f"%h
plt.savefig(f'{file_name}.png', dpi=600)
