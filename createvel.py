import sys
import math
import random
import array as arr
import numpy as np

def createvel():
    kb = 1
    temp = 3
    destemp = float(input("Enter desired reduced temp: "))

    rtemp = math.sqrt(destemp)
    meanke = (1.5*kb*rtemp)/n
    sigmavel = 1/(2*meanke)

    for i in range(0,n):
        vx[i]= rtemp*random.gauss(0,1)
        vy[i]= rtemp*random.gauss(0,1)
        vz[i]= rtemp*random.gauss(0,1)
        continue
    sumx = 0.0
    sumy = 0.0
    sumz = 0.0
    for i in range(0,n):
        sumx = sumx + vx[i]
        sumy = sumy + vy[i]
        sumz = sumz + vz[i]
        continue
    sumx = sumx/float(n)
    sumy = sumy/float(n)
    sumz = sumz/float(n)
    for i in range(0,n):
        vx[i] = vx[i]-sumx
        vy[i] = vy[i]-sumy
        vz[i] = vz[i]-sumz
        print(vx,vy,vz)
        continue
    return

