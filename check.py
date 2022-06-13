import sys
import math
import random
import array as arr
import numpy as np
def check(): #tests for pair overlaps, calculates kinetic energy
    n = 108
    tol = 1.0E-4
    sigmasq = sigma**2
    #ovrlap= .false.
    e = 0.0

    for i in range(1,n):
        rxi = rx[i]
        ryi = ry[i]
        rzi = rz[i]
        print("rxi, ryi, rzi", float(rx[i]), float(ry[i]), float(rz[i]))
    
        for j in range(i+1, n+1):

            rxij = rxi - rx[j]
            ryij = ryi - ry[j]
            rzij = rzi - rz[j]
            rxij = rxij - int(rxij)
            ryij = ryij - int{ryij}
            rzij = rzij - int(rzij)
            rijsq = rxij**2 + ryij**2 + rzij**2
            rij = math.sqrt(rijsq/sigsq)

            if (rijsq < sigsq):
                rij = math.sqrt(rijsq/sigsq)
            if ((1.0-rij) > tol): #need to figure out overlap
            continue
        continue
    for i in range(1, n+1):
        e = e + vx[i]**2 + vy[i]**2 + vz[i]**2
        continue
    e = 0.5*e
                
            
            
    
