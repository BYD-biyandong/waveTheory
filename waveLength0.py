# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 16:19:23 2023

@author: BI-Yandong fron Tongji Univ

function: calculate the wave length L0 under wave-only
"""

import math

def L0(T, h):

    #-------------------------wave parameters--------------------------
    # T  wave period
    # h  water depth 
    
    L0 = 9.81 * math.pow(T,2) / (2 * math.pi)
    L = L0

    for i in range(100):
        Lnew = L0 * math.tanh(2 * math.pi * h / L)
        if abs(Lnew - L) < 0.0001:
            L = Lnew
            break
        else:
            L = Lnew

    return L
