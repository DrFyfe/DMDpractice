import sys
import math
import random
import array as arr
import numpy as np
from check import check
from createcoord import createcoord
from createvel import createvel
from dnlist import dnlist
from uplist import uplist
from grsort import grsort
from bump_subroutine import bump
n = 108 # number of atoms
timbig = 1.0E10 # maximum time
pi = 3.14159265359
maxbin = 500
delr = 0.001

print('PROGRAM SPHERE')
print('Molecular Dynamics of Hard Spheres')
print('Results in Units kt = sigma = 1')


title = input("Enter run title: ")
f = open(title, "w") #trying to have the output file named

n = input("Enter number of spheres:")
density = input("Enter reduced density: ")

ncoll = input("Enter number of collisions required: ") #required user inputs

print("Run title", title)
print("Reduced collision Density is", density)
print("Collisions required", ncoll)

sigma = (density/float(n))**(1.0/3.0)
rx = np.arange(0,n,dtype=float) #setting up arrays for reference in createcoord,vel
ry = np.arange(0,n,dtype=float)
rz = np.arange(0,n,dtype=float)

vx = np.arange(0,n,dtype=float)
vy = np.arange(0,n,dtype=float)
vz = np.arange(0,n,dtype=float)




createcoord()

createvel()

check()

if ovrlap:
    exit()
kecentr = 50
en = e/float(n)
temp = 2.0*en/3.0
print("Temperature", temp)
enkt = en/temp
print("Initial e/nkt", enkt)

for i in range(0, n):
    coltime[i] = timbig
    partnr[i] = n
    continue
for i in range(0,n):
    uplist()
    continue
#zero virial accumulator
f.write('After 1st UPLIST i, coltim(i)',i,coltim[i])

coltime[n+1] = 5.0 #why is it 5?

acw = 0.0

print("Start of Dynamics")

# MAIN LOOP BEGINS
steps = 0
t = 0.0
uplistcntr = 10
grcount = 20

for coll in range(0,ncoll):
    tij = timbig
    for k in range(0,n):
        if(coltim[k] < tij):
            tij = coltim[i]
            i = k
        continue
    if (grcount == coll-1):
        grsort()
        grcount = grcount + 20
        coltim[n+1] = t +5.0
        steps = steps + 1
    else:
        j = partnr[i]
    continue

#**Move particles forward by time tij**
#**and reduce collision times**
#**apply periodic boundaries**
t = t + tij
for k in range(0,n):
    coltim[k] = coltim[k] - tij
    rx[k] = rx[k] + vx[k]*tij
	ry[k] = ry[k] + vy[k]*tij
	rz[k] = rz[k] + vz[k]*tij
	rx[k] = rx[k] - int(rx[k])
	ry[k] = ry[k] - int(ry[k])
	rz[k] = rz[k] - int(rz[k])
    continue

coltim[n+1] = coltim[n+1] - tij

bump()

acw = acw + w

for k in range(0,n):
    if (k == i) or (k == j) or (partnr[k] == j)):
        uplist():
    continue

dnlist() # for i
dnlist() # for j

if (coll == uplistcntr - 1):
    for i in range(1,n+1):
        uplist() # for i
        continue
    uplistcntr = uplist cntr + 10
    continue
if (coll == kecntr):
    for i in range(1,n+1):
        e = e + vx[i]**2+vy[i]**2+vz[i]**2
        e = 0.5*e
        en = e/float(n)
        enkt = en/temp
        f.write("coll and ke and enkt and  w", coll, e, enkt,  w)
        kecntr = kecntr + 50
        continue
print("end of dynamics")
print("Final Colliding Pair", i, j)
#checking for overlaps
check()

if ovrlap:
    print("PARTICLE OVERLAP IN FINAL CONFIGURATION")

pvnkt1 = acw/(float(n)*3.0*t*temp)
en=e/float(n)
enkt = en/temp
t = t*math.sqrt(temp)/sigma
rate = float(ncoll)/t
tbc = float(n)/rate/2.0

f.write('The final n is', n)
f.write('The final e is', e)
f.write('The final temp is', temp)
f.write('The final acw is', acw)
print("Final time is", t)
print("Collision rate is", rate)
print("Mean collision time", tbc)
print("Final e/nkt is",enkt)
print("PV/nkt – 1 is", pvnkt1)
f.write('The final grcount is', grcount)

#calculate g(r) page 184 tildesley

grconst = 4.0*pi*n/3

for bin in range(0,maxbin):
	rlower = REAL(bin-1)*delr
	rupper = rlower + delr
	nideal = grconst*(rupper**3 - rlower**3)
	gr[bin] = float(hist(bin))/float(steps)/float(n)/nideal
	f[bin] = rlower+delr/2
