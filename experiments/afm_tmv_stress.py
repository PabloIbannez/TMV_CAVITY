import os

MAX_ULLINT = 18446744073709551615-1

import numpy as np

import copy
import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time, nanonewton2KcalMol_A_force
ps2AKMA = picosecond2KcalMol_A_time()
nN2AKMA = nanonewton2KcalMol_A_force()

from VLMP.experiments.HighThroughputAFM import HighThroughputAFM

tmvSize = 7

kT = 0.596 # kcal/mol at 300K

#K = 0.05   #N/m = nN/nm from https://pubs.rsc.org/en/content/articlepdf/2016/nr/c6nr01007e
K = 0.10   #Used in this experiment
K = K/10.0 #nN/A
K = K*nN2AKMA #Kcal/mol/A

#
Kxy  = 100.0 #Kcal/mol/A
xy_r = np.sqrt(kT/Kxy)

tipMass =  1e6

dt = 0.01*ps2AKMA

print("xy_r (A) = ", xy_r)

KxyFixing  = 10.0 #Kcal/mol/A
xyFixing_r = np.sqrt(kT/KxyFixing)

print("xyFixing_r (A) = ", xyFixing_r)

tipRadius  = 150.0 #A
tipEpsilon = -1.0 #kcal/mol
tipSigma   =  1.0 #A

surfEpsilon = -1.0 #kcal/mol

thermoSteps = 100000
indSteps    = MAX_ULLINT-thermoSteps #Simulation will stop when the reaches the max force

# -0.001/ps2AKMA is a proper for the tip velocity
tipVelocity = -0.001/ps2AKMA

samples = {}

epsNC     = [1.2]

folderName = str(epsNC[0])+"-"+str(epsNC[-1])
folderName.replace(".","_")
folderName = f"AFM_INDENTATION_eps_{folderName}_SIZE_{tmvSize}_STRESS_T0"

#Check if the folder exists
if os.path.exists(folderName):
    print(f"[ERROR] The folder {folderName} already exists")
    raise ValueError("Folder already exists")

for eps in epsNC:
    name = "eps"+str(eps).replace(".","_")
    samples[f"tmv_{name}"] = {"models":[{"type":"FILE",
                                         "parameters":{"inputFilePath":f"../data/jsons/tmv{tmvSize}_{name}.json"}}],
                              "modelOperations":[{"type":"rotation",
                                                  "parameters":{"axis":[1.0,0.0,0.0],
                                                                "angle":3.141592/2.0,
                                                                "selection":"FILE"}}],
                              "simulationSteps":[{"type":"stressMeasurement",
                                                  "parameters":{"outputFilePath":"stress.dat",
                                                                "intervalStep":100000,
                                                                "radiusCutOff":20.0}}
                                                 ]
                             }

parameters = {
              "simulation":{"units":"KcalMol_A",
                            "types":"basic",
                            "temperature":0.0, # !!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!
                            "box":[2000.0, 2000.0, 2000.0],
                            "samples":copy.deepcopy(samples),
                            "integrator":{"type":"BBK","parameters":{"timeStep":dt,
                                                                     "frictionConstant":0.2/ps2AKMA}}},
              "AFM":{"K":K,
                     "Kxy":Kxy,
                     "epsilon":tipEpsilon,
                     "sigma":tipSigma,
                     "tipVelocity":tipVelocity,
                     "indentationSteps":indSteps,
                     "tipMass":tipMass,
                     "tipRadius":tipRadius,
                     "initialTipSampleDistance":15.0,
                     "indentationPositionX":0.0,
                     "indentationPositionY":0.0
                    },
              "indentation":{
                             "thermalizationSteps":thermoSteps,
                             "indentationSteps":indSteps,
                             "fixSampleDuringThermalization":True,
                             "fixSampleDuringIndentation":True,
                             "KxyFixing":KxyFixing
              },
              "surface":{"epsilon":surfEpsilon,
                         "absorptionHeight":6.0,
                         "absorptionK":10.0
                        },
              "maxForce":{"force":30.0*nN2AKMA,
                          "maxForceIntervalStep":10000},
              "output":{"infoIntervalStep":1000000,
                        "saveStateIntervalStep":100000,
                        "afmMeasurementIntervalStep":10000,
                        "saveStateOutputFilePath":"output",
                        "saveStateOutputFormat":"sp"}
              }

htafm = HighThroughputAFM(parameters)

htafm.generateSimulationPool()
htafm.distributeSimulationPool("one")
htafm.setUpSimulation(folderName)
