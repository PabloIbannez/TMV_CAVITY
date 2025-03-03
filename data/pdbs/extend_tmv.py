import sys
import glob

import MDAnalysis as mda

import numpy as np

N = 9

tmvFile = "2om3.pdb"

########################

u = mda.Universe(tmvFile)

coms = []
for ts in u.trajectory:
    coms.append(u.atoms.center_of_mass())

z = np.mean(np.diff(np.asarray(coms), axis=0), axis=0)[2]
Z = (np.asarray(coms[-1]) - np.asarray(coms[0]))[2] + z

with mda.Writer(f"mergedTMV_{N}.pdb","w") as W:
    with open(f"mergedTMV_{N}.sp","w") as sp:
        for n in range(N):
            u = mda.Universe(tmvFile)

            for ts in u.trajectory:
                u.atoms.positions += np.array([0.0,0.0,n*Z])
                for p in u.atoms.positions:
                    radius = 1.0
                    sp.write(f"{p[0]} {p[1]} {p[2]} {radius} {ts.frame}\n")
                W.write(u.atoms)


