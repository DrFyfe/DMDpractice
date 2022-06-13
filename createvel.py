import sys
import math
import random
import array as arr
import numpy as np

def createvel():
    n = 108
    destemp = input("Enter desired reduced temp: ")

    rtemp = sqrt(destemp)
    vx = np.arange(0,109)
    vy = np.arange(0,109)
    vz = np.arange(0,109)

    for i in range(1,n+1):
        vx[i]= rtemp*random.gauss(0, sigma)
        vy[i]= rtemp*random.gauss(0, sigma)
        vz[i]= rtemp*random.gauss(0, sigma)
        continue
    sumx = 0.0
    sumy = 0.0
    sumz = 0.0
    for i in range(1,n+1):
        sumx = sumx + vx[i]
        sumy = sumy + vy[i]
        sumz = sumz + vz[i]
        continue
    sumx = sumx/float(n)
    sumy = sumy/float(n)
    sumz = sumz/float(n)
    for i in range(1,n+1):
        vx[i] = vx[i]-sumx
        vy[i] = vy[i]-sumy
        vz[i] = vz[i]-sumz
        continue
    return

