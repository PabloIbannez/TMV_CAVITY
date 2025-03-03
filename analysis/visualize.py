#!/home/pablo/py3/bin/python
import sys,os
import json

import numpy as np
import MDAnalysis as mda

import pymol
from pymol import cmd
from pymol.cgo import *

from tqdm import tqdm

@cmd.extend
def takePhoto(fileName:str):
    cmd.show("cartoon")
    cmd.set("ray_trace_mode",10)
    cmd.bg_color("white")
    cmd.set("antialias",5)
    cmd.ray(3840,2160,renderer=0)
    #cmd.ray(1920,1080,renderer=1)
    cmd.png(fileName)

def compute_normal(coord1, coord2, coor3):

    dr21 = coord1-coord2
    dr32 = coord3-coord2

    return np.cross(dr21,dr32)

def genCGOobject(coord1,coord2,coord3,coord4,color):

    #normal1 = compute_normal(coord1, coord2, coord3)
    #normal2 = compute_normal(coord1, coord3, coord4)
    #normal3 = compute_normal(coord2, coord3, coord4)

    normal1 = compute_normal(coord1, coord2, coord4)
    normal2 = compute_normal(coord4, coord3, coord1)

    x1,y1,z1 = coord1
    x2,y2,z2 = coord2
    x3,y3,z3 = coord3
    x4,y4,z4 = coord4

    obj = [

      BEGIN, TRIANGLE_STRIP,

      COLOR, color[0], color[1], color[2],
      NORMAL, normal1[0], normal1[1], normal1[2],
      VERTEX, x1, y1, z1,
      VERTEX, x2, y2, z2,
      VERTEX, x3, y3, z3,
      VERTEX, x4, y4, z4,

      END
    ]

    return obj


# Read input file
if len(sys.argv) != 2:
    print("Usage: visualize.py <simulation_folder>")
    sys.exit()

import __main__

print("Starting pyMol")
#__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
__main__.pymol_argv = [ 'pymol']

pymol.finish_launching()

############################################

surfacePosition = 0.0

#Load surface

coord1 = np.asarray([-1000,-1000,surfacePosition])
coord2 = np.asarray([-1000, 1000,surfacePosition])
coord3 = np.asarray([ 1000,-1000,surfacePosition])
coord4 = np.asarray([ 1000, 1000,surfacePosition])

color = [1.0,1.0,1.0]

cmd.load_cgo(genCGOobject(coord1,coord2,coord3,coord4,color),"plane")

cmd.set_view ((-0.0728457, -0.0590521,  0.9955935,\
               0.9971556,  0.0150494,  0.0738526,\
              -0.0193442,  0.9981415,  0.0577878,\
               0.000354242,    0.001425505, -4989.715332031,\
              -3.765857697,  -87.272483826,  569.281799316,\
              -23755.427734375, 33734.906250000,  -20.000000000 ))

#####################################################################

folder = sys.argv[1]

folderLast = folder.split("/")[-1]

lipid   = float(folderLast.split("_")[2])
protein = float(folderLast.split("_")[4])

input_file = os.path.join(folder, "vesicle.sp")

photoName = f"CORONAVIRUS_lipid_{lipid}_protein_{protein}.png"
photoPath = os.path.join(".", photoName)

#Check if photo already exists
if os.path.isfile(photoPath):
    print(f"Skipping {photoName}... (already exists)")
    sys.exit()

#Check if output.sp exists
if not os.path.isfile(input_file):
    print("File %s does not exist" % input_file)
    sys.exit()

print("Processing %s" % input_file)
#Read the input file
with open(input_file, 'r') as f:

    #Compute number of particles
    count = 0
    for i,line in enumerate(f):
        if '#' in line and i != 0:
            break
        else:
            count += 1
    count -= 1

    #print("Number of particles: ", count)

    #Reset file pointer
    f.seek(0)

    #Write the output file
    with open("tmp.xyz", 'w') as g:
        for i,line in enumerate(f):
            if '#' in line:
                g.write(str(count) + '\n')
                g.write('***\n')
            else:
                x,y,z,r,t = line.split()
                g.write(f"{t} {x} {y} {z} {r}\n")

type2rad = {}
with open("tmp.xyz","r") as f:
    for line in f:
        if len(line.split()) == 5:
            t,_,_,_,r   = line.split()
            type2rad[t] = float(r)

print("Temporary file created for ", folder)

############################################

cmd.load("tmp.xyz", "particles")
os.remove("tmp.xyz")
cmd.show_as("spheres")

type2color = {}
for t in type2rad.keys():
    if t == "0":
        type2color[t] = "blue"
        continue
    if t == "12":
        #type2color[t] = "green"
        type2color[t] = "red"
        continue
    type2color[t] = "red"

for t,r in type2rad.items():
    cmd.alter('name '+t, "vdw={:}".format(r))
    cmd.color(type2color[t],'name '+t)

numberOfStates = cmd.count_states("particles")
cmd.set("state",numberOfStates) #Go to last state

#print("Taking photo")
#Save image at folder
#takePhoto(photoPath)

#print("Image saved at ", photoPath)

#Delete all objects
#cmd.delete("all")

#print("Closing pyMol")
#cmd.quit()



