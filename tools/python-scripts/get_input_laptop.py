#!/usr/bin/env python 
import sys, os
import argparse
from itertools import islice
from subprocess import call
import numpy as np
import math

# usage: 
# python file.py --atoms at1,at2 --fix (at3 or F) --down val --up val --ins val --out val
#                --batch T (or F or nothing)
#
# needs coord.xyz, XDENS, MOL in the same directory


# read input from command line
parser = argparse.ArgumentParser()
parser.add_argument('--atoms', type=str)
parser.add_argument('--fix', type=str)
parser.add_argument('--down', type=float)
parser.add_argument('--up', type=float)
parser.add_argument('--ins', type=float)
parser.add_argument('--out', type=float)
parser.add_argument('--batch', type=str)
args = parser.parse_args() 
var = args.atoms.split(',')
idx = [int(x) for x in var]
# 
if args.fix != "F":
    fixatom = args.fix
    idxfix = int(args.fix)
    print("fixatom", fixatom)
else:
    print("no fixatom set")
#
height = [args.down,args.up]
width = [args.ins,args.out]
print(height)
print(width)
logic = args.batch
print("logic batch", logic)

def mk_sbatch(string,filename):
    f1 = open(filename,'w')
    print("#!/bin/bash", file=f1)
    print(string, file=f1) 
    print(" ", file=f1)
    f1.close()


if logic == "T": 
    l0 ="./mk_slice.py"
    l1 =" --atoms "+str(idx[0])+","+str(idx[1])
    l1f =" --fix "+str(idxfix)
    l2 =" --down "+str(height[0])+" --up "+str(height[1])
    l3 =" --ins "+str(width[0])+" --out "+str(width[1])
    line = l0+l1+l1f+l2+l3

    print(line)
    mk_sbatch(line,"run_slice_gimic")

# assume input coord file in Angstroem
# assume origin for position vectors is 0/0/0
fin = "coord.xyz"
fout = "show_fixpoint.xyz"

#conversion factors
ang2bohr = 1.8897161646320724
bohr2ang = 0.52918

natom = 0
element = []
coord_x = []
coord_y = []
coord_z = []
with open(fin) as f:
    # skip first two lines
    for _ in range(2):
        next(f)
    for i in islice(f, 0, None):
        natom = natom + 1 
        l = i.strip().split()
        element.append(l[0])
        coord_x.append(float(l[1]))
        coord_y.append(float(l[2]))
        coord_z.append(float(l[3]))

print("Number of atoms:", natom)
print("Bond defined by atoms: ", idx[0], idx[1])
print("height is: ", height)
print("width is: ", width)

# A
print("coord atom",idx[0],": ", coord_x[idx[0]-1], coord_y[idx[0]-1], coord_z[idx[0]-1])
# B
print("coord atom",idx[1],": ", coord_x[idx[1]-1], coord_y[idx[1]-1], coord_z[idx[1]-1])
# get bond vector  B-A = vec A->B
dx = coord_x[idx[1]-1] - coord_x[idx[0]-1]
dy = coord_y[idx[1]-1] - coord_y[idx[0]-1]
dz = coord_z[idx[1]-1] - coord_z[idx[0]-1]
# get bond distance
dist = np.sqrt(dx**2 + dy**2 + dz**2)
print("bond distance in A", dist)
print("bond distance in bohr", dist*ang2bohr)
# save 0.5 dist as gimic input in bohr
halfdist_bohr = 0.5*dist*ang2bohr 
halfdist = 0.5*dist
# get midpoint of bond as 0.5*A->B + O->A 
fx = 0.5*dx + coord_x[idx[0]-1] 
fy = 0.5*dy + coord_y[idx[0]-1] 
fz = 0.5*dz + coord_z[idx[0]-1] 
# get radius
rdist = np.sqrt(0.5*dx**2 + 0.5*dy**2 + 0.5*dz**2) 
# get fixpoint 
# 0C = OA + 0.5*AB +/- hc*ABperp/normAB
# AB = (a1,a2) --> ABperp = (-a2,a1) 
fxx = coord_x[idx[0]-1] + 0.5*dx + rdist*(-dy/dist)  
fyy = coord_y[idx[0]-1] + 0.5*dy + rdist*( dx/dist)  
fzz = fz
# if fixpoint is set
if args.fix != "F":
    fxx = coord_x[idxfix-1]
    fyy = coord_y[idxfix-1]
    fzz = coord_z[idxfix-1]
    print("coord fixpoint atom", idxfix, coord_x[idxfix-1], coord_y[idxfix-1], coord_z[idxfix-1])
#
element.append("X")
coord_x.append(fxx)
coord_y.append(fyy)
coord_z.append(fzz)

f1 = open(fout,'w')
print(natom + 1, file=f1)
print(" ", file=f1)
for i in range(natom + 1):
    print(element[i], coord_x[i], coord_y[i], coord_z[i], file=f1)

f1.close()
# sbatch logic
if logic=="T":
    print("running calc...")
elif logic=="F":
    call(["molden", "show_fixpoint.xyz"])

# finput = "gimic.inp"
def write_gimic_input(finput,idx,fxx,fyy,fzz,height,width):
    f1 = open(finput,'w')
    print("# GIMIC INPUT ", file=f1)
    print("calc=integral", file=f1)         
    print('title=""', file=f1)
    print('basis="MOL"', file=f1)
    print('xdens="XDENS"', file=f1)
    print("debug=1", file=f1)           
    print("openshell=false", file=f1) 
    # magnet_axis=X 
    print("magnet=[0.0,0.0,1.0]", file=f1)
    print(" ", file=f1)
    print("Grid(bond) {", file=f1)                   
    print("    type=gauss", file=f1)                
    idx = list(map(str,idx))
    line = ",".join(idx)
    string = '    bond=[',line,']' 
    string = list(map(str,string))
    newline = "".join(string) 
    print(newline, file=f1)
    # convert to bohr
    lst = [fxx*ang2bohr, fyy*ang2bohr, fzz*ang2bohr]
    lst = list(map(str,lst))
    line = ",".join(lst)
    string = '    fixcoord=[',line,']'
    newline = "".join(string) 
    print(newline, file=f1)
    print("    distance=",halfdist_bohr, file=f1)           
    print("    gauss_order=9", file=f1)           
    print("    spacing=[0.02, 0.02, 0.02]", file=f1)  
    height = list(map(str,height))
    line = ",".join(height)
    string = '    height=[',line,']'
    newline = "".join(string) 
    print(newline, file=f1)
    width = list(map(str,width))
    line = ",".join(width)
    string = '    width=[',line,']'
    newline = "".join(string) 
    print(newline, file=f1)
    print("}", file=f1)
    print(" ", file=f1)
    print("Advanced {", file=f1)
    print("    lip_order=5", file=f1)   
    print("    spherical=off", file=f1)
    print("    diamag=on", file=f1)    
    print("    paramag=on", file=f1) 
    print("    GIAO=on ", file=f1)   
    print("    screening=on", file=f1)   
    print("    screening_thrs=1.d-8 ", file=f1)
    print("}", file=f1)

    f1.close()

finput = "gimic.inp"
write_gimic_input(finput,idx,fxx,fyy,fzz,height,width)
call(["gimic", "--dryrun"])
if logic=="T":
    print("running calc...")
elif logic =="F":
    call(["molden", "grid.xyz"])

# if sbatch script had been requested
if logic=="T": 
    calcname = "bond_"+str(idx[0])+"-"+str(idx[1])
    print(calcname)
    os.mkdir(calcname)
    os.chdir(calcname)
    os.system("cp ../MOL .")
    os.system("cp ../XDENS .")
    os.system("cp ../coord.xyz .")
    os.system("cp ../get_input.py .")
    os.system("cp ../mk_slice.py .")
    os.system("cp ../run_slice_gimic .")
    os.system("chmod 715 run_slice_gimic")
    os.system("./run_slice_gimic ")
    os.chdir("../")


