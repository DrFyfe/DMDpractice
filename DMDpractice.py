n = 108 # number of atoms
timbig = 1.0E10 # maximum time
pi = 3.14159265359
maxbin = 500
delr = 0.001

print('PROGRAM SPHERE')
print('Molecular Dynamics of Hard Spheres')
print('Results in Units kt = sigma = 1')


title = input("Enter run title: ")
open(title, "w") #trying to have the output file named

density = input("Enter reduced density: ")

ncoll = input("Enter number of collisions required: ") #required user inputs

def createcoord(rx,ry,rz): #trying to recreate the createcoord subroutine
    n = 108
    nc = 3
    cell = 1.0 / nc
    cell2 = 0.5 * cell
    def rx(n):
        for n in range(5):
            if n == 1:
                rx(n) = 0.0
            elif n == 3:
                rx(n) = 0.0
            else:
                rx(n) = cell2
    def ry(n):
        for n in range(5):
            if n == 1:
                ry(n) = 0.0
            elif n == 4:
                ry(n) = 0.0
            else:
                ry(n) = cell2
    def rz(n):
        for n in range(5):
            if n == 1:
                ry(n) = 0.0
            elif n == 4:
                ry(n) = 0.0
            else:
                ry(n) = cell2
        rx(1)
        ry(1)
        rz(1)
    m = 0
    for iz in range(1, nc):
        for iy in range(1, nc):
            for ix in range(1, nc):
                for iref in range(1,4):
                    #define rx(iref+m) etc
                m = m +4 # I don't understand this. is it necessary?

                
createcoord()
createvel()
check()
if ovrlap is true
    return "Particle overlap in final configuration"

kecntr = 50
en = e/n
temp = 2.0*en/3.0
print("Temperature", temp)
enkt = en/temp
print("Initial e/nkt", enkt)

for i in [1, n]:
    uplist()
    





    
