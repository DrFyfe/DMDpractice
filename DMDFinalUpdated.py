# -*- coding: utf-8 -*-
import sys
import math
import random
import array as arr
import numpy as np
timbig = 1.0E10 # maximum time
pi = 3.14159265359
maxbin = 500
delr = 0.001
e = 0.0
w = 0
print('PROGRAM SPHERE')
print('Molecular Dynamics of Hard Spheres')
print('Results in Units kt = sigma = 1')

title = input("Enter run title: ")
file1 = open(title, "a") #trying to have the output file named
file2 = open("miscinfo.txt","a")
file3 = open("trajfile.txt","a")
n = input("Enter number of spheres:") #number of spheres
n = int(n)
density = input("Enter reduced density: ")
density = float(density)
ncoll = input("Enter number of collisions required: ") #required user inputs
ncoll = int(ncoll)

print("Run title", title)
print("Reduced collision Density is", density)
print("Collisions required", ncoll)

sigma = float((float(density)/int(n))**float(1.0/3.0))
sigsq = sigma**2
print(sigsq)
rx = np.arange(0,n,dtype=float) #setting up arrays for reference in createcoord,vel
ry = np.arange(0,n,dtype=float)
rz = np.arange(0,n,dtype=float)

vx = np.arange(0,n,dtype=float)
vy = np.arange(0,n,dtype=float)
vz = np.arange(0,n,dtype=float)

gr = np.arange(0,maxbin,dtype=float)
f = np.arange(0,maxbin,dtype=float)
hist = np.arange(0,maxbin,dtype=float)
coltim = np.arange(0,n,dtype=float)
partnr = np.arange(0,n,dtype=float)
#possibly make arrays for gr, coltime, etc?

def createcoord(rx,ry,rz): #trying to recreate the createcoord subroutine
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
#        print("rx, ry, rz", float(rx[i]),float(ry[i]),float(rz[i]))
        continue
    return
def grsort(rx,ry,rz):
    #CHange particle number, was originally 0,n and i+1, n+1
    for i in range(0,n-1):
        rxi = rx[i]
        ryi = ry[i]
        rzi = rz[i]
        for j in range(i+1,n):
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
    return

def createvel(vx,vy,vz):
    #kb = 1
    destemp = input("Enter desired reduced temp:")
    destemp = float(destemp)

    rtemp = math.sqrt(destemp)
    #meanke = (1.5*kb*rtemp)/n
    #sigmavel = 1/(2*meanke) possibly no longer needed, keeping around just in case

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
 #       print(vx,vy,vz)
        continue
    return

def check(rx,ry,rz,e): #tests for pair overlaps, calculates kinetic energy
    tol = 1.0E-4
    for i in range(0,n):
        rxi = rx[i]
        ryi = ry[i]
        rzi = rz[i]
#        print("rxi, ryi, rzi", float(rx[i]), float(ry[i]), float(rz[i]))
    
    for j in range(i+1, n):

        rxij = rxi - rx[j]
        ryij = ryi - ry[j]
        rzij = rzi - rz[j]
            #CHECK ANINT VS INT
        rxij = rxij - int(rxij)
        ryij = ryij - int(ryij)
        rzij = rzij - int(rzij)
        rijsq = rxij**2 + ryij**2 + rzij**2
        rij = math.sqrt(rijsq/sigsq)

        if (rijsq < sigsq):
                rij = math.sqrt(rijsq/sigsq)
#            else:
#                continue
        if ((1.0-rij) > tol):
                sys.exit("Overlap Detected")
#            continue
        continue
    for i in range(0, n):
        e = e + vx[i]**2 + vy[i]**2 + vz[i]**2
        print(e)
#        continue
    e = 0.5*e
    return
def dnlist(j,rx,ry,rz,vx,vy,vz): #looks for collisions with atoms i<j
    tij = timbig
    discr = 0
    if (j == 1) or (j >= n):
        return
    rxj = rx[j]
    if j >= len(rx):
        return
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
        ryij = ry[i] - ryj
        rzij = rz[i] - rzj
        rxij = rxij - int(rxij)
        ryij = ryij -  int(ryij)
        rzij = rzij - int(rzij)
        vxij = vx[i] - vxj
        vyij = vy[i] - vyj
        vzij = vz[i] - vzj
        bij = rxij*vxij + ryij*vyij + rzij*vzij

        if (bij < 0.0):
            rijsq = rxij**2 + ryij**2 + rzij**2
            vijsq = vxij**2 + vyij**2 + vzij**2
            discr = bij**2 - vijsq*(rijsq-sigsq)
        if (discr > 0.0):
            tij = (-bij - math.sqrt(discr))/vijsq
        if (tij < coltim[i]):
            coltim[i] = tij
            partnr[i] = j
        continue
    return

def uplist(i,rx,ry,rz,vx,vy,vz):
    tij = timbig
    if (i == n):
        return
    coltim[i] = timbig
    rxi = rx[i]
    ryi = ry[i]
    rzi = rz[i]
    vxi = vx[i]
    vyi = vy[i]
    vzi = vz[i]

    for j in range(i+1,n): #potentially needs range correction
        rxij = rxi - rx[j]
        ryij = ryi - ry[j]
        rzij = rzi - rz[j]
        rxij = rxij - int(rxij) #rounds its argument to the nearest whole number
        ryij = ryij - int(ryij)
        rzij = rzij - int(rzij)
        vxij = vxi - vx[j]
        vyij = vyi - vy[j]
        vzij = vzi - vz[j]
        bij = rxij*vxij + ryij*vyij + rzij*vzij

        if (bij<0.0):
            rijsq = rxij**2 + ryij**2 + rzij**2
            vijsq = vxij**2 + vyij**2 + vzij**2
            discr = bij**2 - vijsq*(rijsq - sigsq)
            if (discr>0.0):
                tij = (-bij - math.sqrt(discr))/vijsq
                if (tij<coltim[i]):
                    coltim[i] = tij
                    partnr[i] = j
        continue
    return

def bump(w,rx,ry,rz,vx,vy,vz):
    rxij = rx[i] - rx[j]
    ryij = ry[i] - ry[j]
    rzij = rz[i] - rz[j]
    rxij = rxij - int(rxij)
    ryij = ryij - int(ryij)
    rzij = rzij - int(rzij)

    factor = rxij*(vx[i]-vx[j]) + ryij*(vy[i]-vy[j]) + rzij*(vz[i]-vz[j])/sigsq

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
    return w
print("Making Coordinates")
createcoord(ry,ry,rz)
print("Done coordinates, making vels")
createvel(vx,vy,vz)
print("Done vels, checking")
check(rx,ry,rz,e)
print("Done checking, doing loop")
kecentr = 50
en = e/float(n)
temp = 2.0*en/3.0
print(e)
print("Temperature", temp)
enkt = en/temp
print("Initial e/nkt", enkt)

for i in range(0, n):
    coltim[i] = timbig
    partnr[i] = n
    continue
for i in range(0,n):
    uplist(i,rx,ry,rz,vx,vy,vz)
    continue
#zero virial accumulator
file2.write('After 1st UPLIST i, coltim(i)')
file2.write(str(i))
file2.write(str(coltim[i]))

coltim[n-1] = 5.0 #why is it 5?

acw = 0.0

print("Start of Dynamics")

# MAIN LOOP BEGINS
steps = 0
t = 0.0
uplistcntr = 10
grcount = 20

for coll in range(0,ncoll): #main loop
    print("Collision Number: ", coll)
    tij = timbig
    for k in range(0,n):
        if k >= len(coltim):
            break
        if(coltim[k] < tij):
            tij = coltim[i]
            i = k
        continue
    if (grcount == coll-1):
        grsort(rx,ry,rz)
        grcount = grcount + 20
        if n >= len(coltim):
            break
        coltim[n] = t +5.0
        steps = steps + 1
    else:
        j = partnr[i]
    steps = steps + 1 #not part of original code, test
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
#was originally n+1
    coltim[n-1] = coltim[n-1] - tij


    bump(w,rx,ry,rz,vx,vy,vz)
    print("Bump happens at coll no.: ", coll)
    acw = acw + w

    for k in range(0,n):
        if (k == i) or (partnr[k] == i) or (k == j) or (partnr[k] == j):
            uplist(k,rx,ry,rz,vx,vy,vz) #is k the correct input?

    dnlist(i,rx,ry,rz,vx,vy,vz) # for i
    dnlist(j,rx,ry,rz,vx,vy,vz) # for j

    if (coll == uplistcntr - 1):
        for i in range(0,n):
            uplist(i,rx,ry,rz,vx,vy,vz) # for i
        uplistcntr = uplistcntr + 10
    if (coll == kecentr):
        for i in range(0,n):
            e = e + vx[i]**2+vy[i]**2+vz[i]**2
        e = 0.5*e
        en = e/float(n)
        enkt = en/temp
        e = str(e)
        coll = str(coll)
        enkt = str(enkt)
        w = str(w)
        acw = str(acw)
        temp = str(temp)
        t = str(t) # converting for string for formatting reasons
        file2.write(coll + " " + e + " " + enkt + " " + w + " " + acw + " " + temp + " " + t)
        file2.write("%i\n" % en)
        file2.write("%i\n" % acw)
        file2.write("%i\n" % temp)
        file2.write("%i\n" % t) #consider indent
        file3.write("ITEM: TIMESTEP\n" + t + "\n")
        file3.write("ITEM: NUMBER OF ATOMS\n1200")
        file3.write("ITEM: BOX BOUNDS\n0	1\n0	1\n0	1")
        file3.write("ITEM: ATOMS index type x y z" + "\n")
        for n in range(0,n):
            n = str(n)
            file3.write("     " + n + "  1" + "           ")
            file3.write('%.3f' % rx[n] + "           ")
            file3.write('%.3f' % ry[n] + "           ")
            file3.write('%.3f' % rz[n] + "\n")
        e = float(e)
        coll = float(coll)
        enkt = float(enkt)
        w = float(w)
        acw = float(acw)
        temp = float(temp)
        t = float(t) # returning to float format for any potential operations
        check(rx,ry,rz,e)
    
        kecntr = kecntr + 50
 #end of main loop
print("end of dynamics")
print("Final Colliding Pair", i, j)
#checking for overlaps
check(rx,ry,rz,e)
t = float(t)
pvnkt1 = acw/(float(n)*3.0*t*temp)
en=e/float(n)
enkt = en/temp
t = t*math.sqrt(temp)/sigma
rate = float(ncoll)/t
tbc = float(n)/rate/2.0

file1.write('The final n is')
file1.write(str(n))
file1.write('The final e is')
file1.write(str(e))
file1.write('The final temp is')
file1.write(str(temp))
file1.write('The final acw is')
file1.write(str(acw))
print("Final time is", t)
print("Collision rate is", rate)
print("Mean collision time", tbc)
print("Final e/nkt is",enkt)
print("PV/nkt  1 is", pvnkt1)
file1.write('The final grcount is')
file1.write(str(grcount))

#calculate g(r) page 184 tildesley

grconst = 4.0*pi*n/3

for bin in range(0,maxbin):
    rlower = float(bin-1)*delr
    rupper = rlower + delr
    nideal = grconst*(rupper**3 - rlower**3)
#    print(steps,nideal)
    gr[bin] = float(hist[bin])/float(steps)/float(n)/nideal
    f[bin] = rlower+delr/2
    continue
file1 = open("gr.txt", "w")
for bin in range(0,maxbin):
    file1.write("%i\n" % float((f[bin])/(sigma)))
    file1.write("%i\n" % float(gr[bin]))
    continue
#end of program
