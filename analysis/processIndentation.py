import sys
import numpy as np

from VLMP.utils.units import KcalMol_A_force2nanonewton

data = np.loadtxt(sys.argv[1], skiprows=1)

indStart    = 100000
indInterval = 10000

saveInterval = 100000

def processStep(s):
    # s is an indentation step, we have to find the equivalent frame (saved every saveInterval)
    frame = int(s / saveInterval)
    return frame

def processIndentation(ind):
    # From A to nm
    return ind * 0.1

def processForce(force):
    # From AKMA to pN
    return KcalMol_A_force2nanonewton(force)


for d in data:
    step,ind,force = d[0],d[1],d[2]
    print(processStep(step), processIndentation(ind), processForce(force))
