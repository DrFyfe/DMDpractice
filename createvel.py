import sys
import math
import random
import array as arr
import numpy as np

def createvel():
    n = 108
    kb = 1
    temp = 3
    destemp = float(input("Enter desired reduced temp: "))

    rtemp = math.sqrt(destemp)
    meanke = (1.5*kb*rtemp)/n
    sigmavel = 1/(2*meanke)
    vx = np.arange(0,109,dtype=float)
    vy = np.arange(0,109,dtype=float)
    vz = np.arange(0,109,dtype=float)

    for i in range(1,n+1):
        vx[i]= rtemp*random.gauss(0,sigmavel)
        vy[i]= rtemp*random.gauss(0,sigmavel)
        vz[i]= rtemp*random.gauss(0,sigmavel)
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
        print(vx,vy,vz)
        continue
    return vx,vy,vz
createvel()

