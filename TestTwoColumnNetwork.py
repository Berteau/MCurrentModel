from TwoColumnSimulation import *
# import winsound
import dill
from copy import deepcopy

sys.setrecursionlimit(10000)

params = {}
# Time Parameters
params["maxTime"] = 100 #1000
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


sim = TwoColumnSimulation(params)
# weightMatrixPrior = zeros([params["popCount"], params["popCount"]])
# for x in range(len(sim.network.populations["InputA"].cells)):
#     for y in range(len(sim.network.populations["pyramidalsA"].cells)):
#         for source in sim.network.populations["InputA"].cells[x].outputs:
#             # print(source.target.name, sim.network.populations["pyramidalsA"].cells[y].name)
#             if source.target == sim.network.populations["pyramidalsA"].cells[y]:
#                 weightMatrixPrior[x,y] = source.postSynapticReceptors[0].weight
#                 # print("Found a match!")
# figure()
# pcolor(weightMatrixPrior)
# colorbar()
# title('Prior Weights from Input A to Pyramidals A')
#
# weightMatrixPriorBA = zeros([params["popCount"], params["popCount"]])
# for x in range(len(sim.network.populations["InputB"].cells)):
#     for y in range(len(sim.network.populations["pyramidalsA"].cells)):
#         for source in sim.network.populations["InputB"].cells[x].outputs:
#             if source.target == sim.network.populations["pyramidalsA"].cells[y]:
#                 weightMatrixPriorBA[x,y] = source.postSynapticReceptors[0].weight
#
# figure()
# pcolor(weightMatrixPriorBA)
# colorbar()
# title('Prior Weights from Input B to Pyramidals A')

sim.run()

#
# weightMatrixPost = zeros([params["popCount"], params["popCount"]])
# for x in range(len(sim.network.populations["InputA"].cells)):
#     for y in range(len(sim.network.populations["pyramidalsA"].cells)):
#         for source in sim.network.populations["InputA"].cells[x].outputs:
#             if source.target == sim.network.populations["pyramidalsA"].cells[y]:
#                 weightMatrixPost[x,y] = source.postSynapticReceptors[0].weight
#
# figure()
# pcolor(weightMatrixPost)
# colorbar()
# title('Posterior Weights from Input A to Pyramidals A')
#
# weightMatrixPostBA = zeros([params["popCount"], params["popCount"]])
# for x in range(len(sim.network.populations["InputB"].cells)):
#     for y in range(len(sim.network.populations["pyramidalsA"].cells)):
#         for source in sim.network.populations["InputB"].cells[x].outputs:
#             if source.target == sim.network.populations["pyramidalsA"].cells[y]:
#                 weightMatrixPostBA[x,y] = source.postSynapticReceptors[0].weight
#
# figure()
# pcolor(weightMatrixPostBA)
# colorbar()
# title('Posterior Weights from Input B to Pyramidals A')
#
# figure()
# pcolor(weightMatrixPostBA - weightMatrixPriorBA)
# colorbar()
# title('Change in BA weights')

with open("TwoColumnPickleSingleRun_NP.jar", "wb") as pickleJar:
    dill.dump(deepcopy(sim), pickleJar)

# with open("TwoColumnPickle.jar", "rb") as pickleJar:
#     sim = dill.load(pickleJar)

sim.plotColumns()
print("FINISHED!")

firstSecondSpikeTotal = 0
for j in range(10):
        firstSecondSpikes = [sim.network.populations["InputA"].cells[j].spikeRecord[i]  for i in range(len(sim.network.populations["InputA"].cells[j].spikeRecord)) if sim.network.populations["InputA"].cells[j].spikeRecord[i]< 1000]
        firstSecondSpikeTotal += len(firstSecondSpikes)
print(firstSecondSpikeTotal)
