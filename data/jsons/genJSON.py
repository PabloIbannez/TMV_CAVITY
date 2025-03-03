import pyGrained.models.AlphaCarbon as proteinModel
import pyUAMMD

import copy

sopParamsBase = {"SASA":False,
                 "centerInput":False,
                 "aggregateChains":True,
                 "parameters":{}}

SIZE = 3
epsilonNC = [1.0,1.1,1.2,1.3,1.4]

names = []
pdbs  = []
for eps in epsilonNC:
    names.append("tmv"+str(SIZE)+"_eps"+str(eps).replace(".","_"))
    pdbs.append("../pdbs/mergedTMV_"+str(SIZE)+".pdb")

for name,pdb,eps in zip(names,pdbs,epsilonNC):
    sopParams = sopParamsBase.copy()

    sopParams["parameters"]["epsilonNC"] = eps

    sop = proteinModel.SelfOrganizedPolymer(name = name,
                                            inputPDBfilePath = pdb,
                                            params = sopParams)

    types = {"type":["Types","Basic"],
             "labels":["name","mass","radius","charge"],
             "data":[]}

    for t,tinfo in sop.getTypes().items():
        types["data"].append([tinfo["name"],tinfo["mass"],tinfo["radius"],tinfo["charge"]])

    state      = sop.getState()
    structure  = sop.getStructure()
    forceField = sop.getForceField()

    sim = pyUAMMD.simulation()

    sim["global"] = {}
    sim["global"]["types"] = copy.deepcopy(types)

    sim["state"]           = copy.deepcopy(state)

    sim["topology"] = {}
    sim["topology"]["structure"]  = copy.deepcopy(structure)
    sim["topology"]["forceField"] = copy.deepcopy(forceField)

    sim.write(name+".json")
