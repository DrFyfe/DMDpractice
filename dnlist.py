import sys
import math
import random
import array as arr
import numpy as np

def dnlist(): #looks for collisions with atoms i<j
    n - 108
    timbig = 1.0E10

    if (j == 1):
        return
    sigsq = sigma**2
    rxj = rx[j]
    ryj = ry[j]
    rzj = rz[j]
    vxj = vx[j]
    vyj = vy[j]
    vzj = vz[j]

    for i in range(1, j):
        rxij = rx[i] - rxj
        ryij = ry[i] - ryj
        rzij = rz[i] - rzj
        rxij = rxij - int(rxij)
        ryij = ryij -  int(ryij)
        rzij = rzij - int(rzij)
        vxij = vx[i] - vxj
        vyij = vy[i] - vyj
        vzij = vz[i] - vzj
        bij = rxij*vxij + ryij*vyij + rzij*vzij
        ryij = ry(i) - ryj
        rzij = rz(i) - rzj
        rxij = rxij - int(rxij)
        ryij = ryij -  int(ryij)
        rzij = rzij - int(rzij)
        vxij = vx(i) - vxj
        vyij = vy(i) - vyj
        vzij = vz(i) - vzj
        bij = rxij*vxij + ryij*vyij + rzij*vzij

        if (bij < 0.0):
            rijsq = rxij**2 + ryij**2 + rzij**2
            vijsq = vxij**2 + vyij**2 + vzij**2
            discr = bij**2 - vijsq*(rijsq-sigsq)
        if (discr > 0.0):
            tij = (-bij - SQRT(discr))/vijsq
        if (tij < coltim[i]):
            coltim[i] = tij
            partnr[i] = j

    return

            




    
