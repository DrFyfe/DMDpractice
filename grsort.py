import sys
import math
import random
import array as arr
import numpy as np

def grsort():
    n = 108
    maxbin = 500
    delr = 0.001

    for i in range(i,n):
        rxij = rxi - rx[j]
        ryij = ryi - ry[j]
        rzij = rzi - rz[j]
        rxij = rxij - int(rxij)
        ryij = ryij - int(ryij)
        rzij = rzij - int(rzij)

        rijsq = rxij**2 + ryij**2 + rzij**2
        rij = math.sqrt(rijsq)
        bin = int(rij/delr)+1

        if (bin < maxbin):
            hist[bin] = hist[bin] + 2
        continue
    continue
