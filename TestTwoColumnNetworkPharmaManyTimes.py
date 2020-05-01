from TwoColumnSimulationPharma import *
# import winsound
import dill
from copy import deepcopy

sys.setrecursionlimit(10000)

params = {}
# Time Parameters
params["maxTime"] = 1000
params["tau"] = 0.1

# Population and Connection Parameters
params["popCount"] = 10
params["pyramidalSelfExcitationWeight"] = 500
params["pyramidalToPyramidalWeight"] = 500
params["pyramidalToPyramidalLikelihood"] = 0.3
params["PyramidalsToFSWeight"] = 500
params["FSToPyramidalsWeight"] = -1000
params["PyramidalsToLTSWeight"] = 100 #500
params["LTStoFSWeight"] = -500
params["LTStoPyramidalsWeight"] = -500

## Regular
# Input Paramters
params["inputWeightA"] = 500
params["inputWeightB"] = 500
params["rateA"] = 10
params["rateB"] = 10
params["inputWeightAB"] = 250
params["crossModalABLikelihood"] = 0.5
params["inputWeightBA"] = 250 #7500
params["crossModalBALikelihood"] = 0.5

# Diffuse Neurotransmitter Paramters
params["serotoninLevelA"] = 10
params["serotoninLevelB"] = 10
params["Somatic5HT2AWeight"] = 100
params["Somatic5HT2AWeightLTS"] = 85
params["Somatic5HT1AWeight"] = -80
params["Axonal5HT2AWeight"] = 0.2
params["Axonal5HT1AWeight"] = -0.2

for i in range(10):
    sim = TwoColumnSimulationPharma(params)
    sim.run()
    with open("TwoColumnPicklePharma"+str(i)+".jar", "wb") as pickleJar:
        dill.dump(deepcopy(sim), pickleJar)
    # sim.plotColumns()

print("FINISHED!")
