n = 108 # number of atoms
timbig = 1.0E10 # maximum time
pi = 3.14159265359 #duh
maxbin = 500 #still not sure what this is
delr= 0.001

print('PROGRAM SPHERE')
print('Molecular Dynamics of Hard Spheres')
print('Results in Units kt = sigma = 1')


title = input("Enter run title: ")
open(title, "w") #trying to have the output file named

density = input("Enter reduced density: ")

ncoll = input("Enter number of collisions required: ") #required user inputs

def createcoord(rx,ry,rz):
    n = 108
    nc = 3
    cell = 1.0 / nc
    cell2 = 0.5 * cell
    def rx(n):
        n = range(5)
        if n == 1:
            return 0.0
        elif n == 3:
            return 0.0
        else:
            return cell2


    