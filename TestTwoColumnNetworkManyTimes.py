from TwoColumnSimulation import *
# import winsound
import dill
from copy import deepcopy

sys.setrecursionlimit(10000)

params = {}
# Time Parameters
params["maxTime"] = 1000 #1000
params["tau"] = 0.1

# Population and Connection Parameters
params["popCount"] = 20

params["pyramidalSelfExcitationWeight"] = 25
params["pyramidalToPyramidalWeight"] = 50
params["pyramidalToPyramidalLikelihood"] = 0.3
params["PyramidalsToFSWeight"] = 35
params["FSToPyramidalsWeight"] = -100
params["PyramidalsToLTSWeight"] = 20 #500
params["LTStoFSWeight"] = -50
params["LTStoPyramidalsWeight"] = -50

## Regular
# Input Paramters
params["inputWeightA"] = 50
params["inputWeightB"] = 50
params["rateA"] = 1
params["rateB"] = 1
params["inputWeightAB"] = 25
params["crossModalABLikelihood"] = 0.5
params["inputWeightBA"] = 25 #7500
params["crossModalBALikelihood"] = 0.5

# Diffuse Neurotransmitter Paramters
params["serotoninLevelA"] = 10
params["serotoninLevelB"] = 10
params["Somatic5HT2AWeight"] = 85
params["Somatic5HT2AWeightLTS"] = 81
params["Somatic5HT1AWeight"] = -80
params["Axonal5HT2AWeight"] = 0.4
params["Axonal5HT1AWeight"] = -0.4

for i in range(10):
    sim = TwoColumnSimulation(params)
    sim.run()
    with open("TwoColumnPickle"+str(i)+"_NP.jar", "wb") as pickleJar:
        dill.dump(deepcopy(sim), pickleJar)
    # sim.plotColumns()

print("FINISHED!")
