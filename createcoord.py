import sys
import math
import random
import array as arr
import numpy as np

def createcoord(): #trying to recreate the createcoord subroutine
    n = 108
    nc = 3
    cell = 1.0 / nc
    cell2 = 0.5 * cell
    rx = np.arange(0,109,dtype=float)
    ry = np.arange(0,109,dtype=float)
    rz = np.arange(0,109,dtype=float)
    
    rx[1] = 0.0   #subblatice a
    ry[1] = 0.0
    rz[1] = 0.0

    rx[2] = cell2 #sublattice b
    ry[2] = cell2
    rz[2] = 0.0

    rx[3] = 0.0   #sublattice c
    ry[3] = cell2
    rz[3] = cell2

    rx[4] = cell2 #sublattice d
    ry[4] = 0.0
    rz[4] = cell2
    for m in range(0,105,4):
        for iz in range(1,nc+1):
            for iy in range(1,nc+1):
                for ix in range(1,nc+1):
                    for iref in range(1,5):
                        rx[iref + m] = rx[iref] + cell*(ix-1)
                        ry[iref + m] = ry[iref] + cell*(iy-1)
                        rz[iref + m] = rz[iref] + cell*(iz-1)
                        continue
                    continue
                continue
            continue
        continue
    for i in range(1,n+1):
        rx[i] = rx[i] - 9.5
        ry[i] = ry[i] - 9.5
        rz[i] = rz[i] - 9.5
        print("rx, ry, rz", float(rx[i]),float(ry[i]),float(rz[i]))
        continue
    return(rx)
createcoord()
print(rx)
