import sys
import math
import random
import array as arr
import numpy as np
n = 108
rx = np.arange(0,n,dtype=float)
ry = np.arange(0,n,dtype=float)
rz = np.arange(0,n,dtype=float)

def createcoord(): #trying to recreate the createcoord subroutine
    nc = 3
    cell = 1.0 / nc
    cell2 = 0.5 * cell
    
    rx[0] = 0.0   #subblatice a
    ry[0] = 0.0
    rz[0] = 0.0

    rx[1] = cell2 #sublattice b
    ry[1] = cell2
    rz[1] = 0.0

    rx[2] = 0.0   #sublattice c
    ry[2] = cell2
    rz[2] = cell2

    rx[3] = cell2 #sublattice d
    ry[3] = 0.0
    rz[3] = cell2
    for m in range(0,n-4,4):
        for iz in range(0,nc):
            for iy in range(0,nc):
                for ix in range(0,nc):
                    for iref in range(0,3):
                        rx[iref + m] = rx[iref] + cell*(ix-1)
                        ry[iref + m] = ry[iref] + cell*(iy-1)
                        rz[iref + m] = rz[iref] + cell*(iz-1)
                        continue
                    continue
                continue
            continue
        continue
    for i in range(0,n):
        rx[i] = rx[i]-0.5
        ry[i] = ry[i]-0.5
        rz[i] = rz[i]-0.5
        print("rx, ry, rz", float(rx[i]),float(ry[i]),float(rz[i]))
        continue
    return
createcoord()
