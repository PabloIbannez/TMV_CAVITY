import sys

from VLMP.experiments.HighThroughputAFM import AnalysisAFM

analysis = AnalysisAFM(sys.argv[1],"nN_nm",maxForce = 4.0,plotTime = True, frameRate = 100000)
analysis.run()


