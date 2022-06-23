import sys
import math
import random
import array as arr
import numpy as np

def bump():
    n = 108

    sigsq = sigma**2
    rxij = rx[i] - rx[j]
    ryij = ry[i] - ry[j]
    rzij = rz[i] - rz[j]
    rxij = rxij - int(rxij)
    ryij = ryij - int(ryij)
    rzij = rzij - int(rzij)

    factor = (rxij*(vx[i] - vx[j]) + ryij*(vy[i] - vy[j]) + rzij*(vz[i] - vz[j]))/sigsq

    delvx = -factor*rxij
    delvy = -factor*ryij
    delvz = -factor*rzij

    vx[i] = vx[i] + delvx
    vx[j] = vx[j] - delvx
    vy[i] = vy[i] + delvy
    vy[j] = vy[j] - delvy
    vz[i] = vz[i] + delvz
    vz[j] = vz[j] - delvz
    
    w = delvx*rxij + delvy*ryij + delvz*rzij

    return
