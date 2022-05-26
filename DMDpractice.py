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

def createcoord(): #trying to recreate the createcoord subroutine
    n = 108
    nc = 3
    cell = 1.0 / nc
    cell2 = 0.5 * cell
    def rx(n):
        for n in range(1,4):
            if n == 1:
                rx = 0.0
            elif n == 3:
                rx = 0.0
            else:
                rx = cell2
    def ry(n):
        for n in range(1,4):
            if n == 1:
                ry = 0.0
            elif n == 4:
                ry = 0.0
            else:
                ry = cell2
    def rz(n):
        for n in range(1,4):
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
createcoord()
createvel()
check()
if ovrlap is true:
    print("Particle overlap in final configuration")

kecntr = 50
en = e/n
temp = 2.0*en/3.0
print("Temperature", temp)
enkt = en/temp
print("Initial e/nkt", enkt)

for i in [1, n]:
    uplist()
    





    
