import sys

import json

import numpy as np
import matplotlib.pyplot as plt


startForce = 1.3
#collapsePosition = 4.5
collapsePosition = 5.0

selectedX = [0.0,2.0,collapsePosition,7.0,15.0]

simFolder = sys.argv[1]

simulation  = simFolder + "/simulation.json"
indentation = simFolder + "/afm.dat"
output      = simFolder + "/output.sp"

with open(simulation) as f:
    sim = json.load(f)

indentationStep = sim["simulationStep"]["afmMeasurement"]["parameters"]["intervalStep"]
outputStep      = sim["simulationStep"]["saveState"]["parameters"]["intervalStep"]

print("Reading indentation data from: "+indentation,". Indentation step: ",indentationStep)
print("Reading output data from: "+output,". Output step: ",outputStep)

data = np.loadtxt(indentation)

steps = data[:,0]
X     = data[:,1]/10.0
F     = data[:,2]

# Invert data respect to y axis
X     = -X
steps = steps[::-1]

# Shift data to start at 0
for i in range(len(X)):
    if F[i] > startForce:
        break
X     = X - X[i]

# Found the closest step for the selected X
selectedStep = []
for selx in selectedX:
    for i,x in enumerate(X):
        if x > selx:
            selectedStep.append(i*indentationStep)
            break

selectedOutputStep = []
for selx in selectedStep:
    selectedOutputStep.append(int(selx/outputStep))

for sel in selectedOutputStep:
    outputPath = simFolder+"/frame"+str(sel)+".sp"
    with open(outputPath,"w") as frame:
        frame.write("#\n")

with open(output) as f:
    framesCount = 0
    for line in f:
        if framesCount in selectedOutputStep:
            with open(simFolder+"/frame"+str(framesCount)+".sp","a") as frame:
                if "#" not in line:
                    frame.write(line)
        if "#" in line:
            framesCount += 1
            continue

plt.plot(X, F, 'k-')
plt.xlabel('Indentation (nm)')
plt.ylabel('Force (nN)')
plt.title('TMV indentation collapse')
plt.axvline(x=collapsePosition, color='r', linestyle='--')
plt.savefig(simFolder+"/collapse.png")
