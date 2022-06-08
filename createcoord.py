import sys
import math
import random

def createcoord(): #trying to recreate the createcoord subroutine
    nc = 3
    cell = 1.0 / nc
    cell2 = 0.5 * cell
    def rx(n): #creating the unit cell for x coordinates
        if n == 1 or n ==3:
            rx = 0.0
        elif n == 2 or n == 4:
            rx = cell2
        return rx
    def ry(n): #creating the unit cell for y coordinates
        if n == 1 or n == 4:
            ry = 0.0
        elif n == 2 or n == 3:
            ry = cell2
        return ry
    def rz(n): #creating the unit cell for z coordinates
        if n == 1 or n == 2:
            ry = 0.0
        elif n == 3 or n ==4:
            rz = cell2
        return rz      
    m = 0
    for iz in range(1, nc+1):
        for iy in range(1, nc+1):
            for ix in range(1, nc+1):
                for iref in range(1,5):
                    rx(iref + m) = rx(iref) + cell*(ix-1)
                    ry(iref + m) = ry(iref) + cell*(iy-1)
                    rz(iref + m) = rz(iref) + cell*(iz-1)
                continue
                m = m + 4
            continue
        continue
createcoord()

