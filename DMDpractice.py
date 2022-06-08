import sys
import math
import random
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

density = input("Enter reduced density: ")

ncoll = input("Enter number of collisions required: ") #required user inputs

def createcoord(): #trying to recreate the createcoord subroutine
    n = 108
    nc = 3
    cell = 1.0 / nc
    cell2 = 0.5 * cell
    def rx(n):
        for n in range(1,5):
            if n == 1:
                rx = 0.0
            elif n == 3:
                rx = 0.0
            else:
                rx = cell2
    def ry(n):
        for n in range(1,5):
            if n == 1:
                ry = 0.0
            elif n == 4:
                ry = 0.0
            else:
                ry = cell2
    def rz(n):
        for n in range(1,5):
            if n == 1:
                ry = 0.0
            elif n == 4:
                ry = 0.0
            else:
                ry = cell2
    m = 0
    for iz in range(1, nc):
        for iy in range(1, nc):
            for ix in range(1, nc):
                for iref in range(1,4):
                    rx = rx(iref) + cell*(ix-1)
                    ry = ry(iref) + cell*(iy-1)
                    rz = rz(iref) + cell*(iz-1) #end subroutine createcoord()
def createvel(): #begin createvel subroutine
    n = 108
    destemp = input("Enter desired reduced temp: ")
    rtemp = sqrt(destemp)
    def gauss():
        a1 = 3.949846138
        a3 = 0.252408784
        a5 = 0.076542912
        a7 = 0.008355968
        a9 = 0.029899776
        rsum = 0
        for i in range(1,12): #how is i incomporated into this?
            rsum = rsum + random()
        r = (rsum - 6)/4
        r2 = r*r #need to finish this
            
    for i in range(1, n):
        vx(i) = rtemp * gauss ()
        vy(i) = rtemp * gauss ()
        vz(i) = rtemp * gauss () #very much unfinished
        
def check():
    n = 108
    tol = 1.0E-4
    sigsq = sigma**2
    ovrlap = .flase
    e = 0.0
    for i in range(1,n-1):
        rxi = rx(i)
        ryi = ry(i)
        rzi = rz(i)
    print("rxi, ryi, rzi", #need to see how to end this
        
createcoord()
createvel()
check()
if ovrlap is true:
    sys.exit("Particle overlap in final configuration")

kecntr = 50
en = e/n
temp = 2.0*en/3.0
print("Temperature", temp)
enkt = en/temp
print("Initial e/nkt", enkt)

for i in [1,n]:
          coltim = timbig
          partnr = n
for i in [1,n]:
          uplist()
f.write('After 1st UPLIST i, coltime()i' + i + coltime(i))





    
