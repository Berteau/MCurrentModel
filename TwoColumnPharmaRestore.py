from TwoColumnSimulationPharma import *
# import winsound
import dill
from copy import deepcopy

sys.setrecursionlimit(10000)

# params = {}
# # Time Parameters
# params["maxTime"] = 1000
# params["tau"] = 0.1

# # Population and Connection Parameters
# params["popCount"] = 10
# params["pyramidalSelfExcitationWeight"] = 500
# params["pyramidalToPyramidalWeight"] = 500
# params["pyramidalToPyramidalLikelihood"] = 0.3
# params["PyramidalsToFSWeight"] = 500
# params["FSToPyramidalsWeight"] = -1000
# params["PyramidalsToLTSWeight"] = 500
# params["LTStoFSWeight"] = -500
# params["LTStoPyramidalsWeight"] = -500

# ## Regular
# # Input Paramters
# params["inputWeightA"] = 500
# params["inputWeightB"] = 500
# params["rateA"] = 10
# params["rateB"] = 10
# params["inputWeightAB"] = 250
# params["crossModalABLikelihood"] = 0.5
# params["inputWeightBA"] = 250 #7500
# params["crossModalBALikelihood"] = 0.5

# # Diffuse Neurotransmitter Paramters
# params["serotoninLevelA"] = 10
# params["serotoninLevelB"] = 10
# params["Somatic5HT2AWeight"] = 100
# params["Somatic5HT2AWeightLTS"] = 100
# params["Somatic5HT1AWeight"] = -80
# params["Axonal5HT2AWeight"] = 0.005
# params["Axonal5HT1AWeight"] = -0.005

# # ## Deprivation
# # # Input Paramters
# # params["inputWeightA"] = 500
# # params["inputWeightB"] = 500
# # params["rateA"] = 0
# # params["rateB"] = 10
# # params["inputWeightAB"] = 250
# # params["crossModalABLikelihood"] = 0.3
# # params["inputWeightBA"] = 250 #7500
# # params["crossModalBALikelihood"] = 0.3
# #
# # # Diffuse Neurotransmitter Paramters
# # params["serotoninLevelA"] = 10
# # params["serotoninLevelB"] = 10
# # params["Somatic5HT2AWeight"] = 100
# # params["Somatic5HT2AWeightLTS"] = 80
# # params["Somatic5HT1AWeight"] = -80
# # params["Axonal5HT2AWeight"] = 0.2
# # params["Axonal5HT1AWeight"] = -0.2

# # # 5HT
# # # Input Paramters
# # params["inputWeightA"] = 500
# # params["inputWeightB"] = 500
# # params["rateA"] = 0
# # params["rateB"] = 10
# # params["inputWeightAB"] = 250
# # params["crossModalABLikelihood"] = 0.3
# # params["inputWeightBA"] = 250 #7500
# # params["crossModalBALikelihood"] = 0.3
# #
# # # Diffuse Neurotransmitter Paramters
# # params["serotoninLevelA"] = 40
# # params["serotoninLevelB"] = 10
# # params["Somatic5HT2AWeight"] = 100
# # params["Somatic5HT2AWeightLTS"] = 80
# # params["Somatic5HT1AWeight"] = -80
# # params["Axonal5HT2AWeight"] = 0.2
# # params["Axonal5HT1AWeight"] = -0.2

# sim = TwoColumnSimulationPharma(params)
# # weightMatrixPrior = zeros([params["popCount"], params["popCount"]])
# # for x in range(len(sim.network.populations["InputA"].cells)):
# #     for y in range(len(sim.network.populations["pyramidalsA"].cells)):
# #         for source in sim.network.populations["InputA"].cells[x].outputs:
# #             # print(source.target.name, sim.network.populations["pyramidalsA"].cells[y].name)
# #             if source.target == sim.network.populations["pyramidalsA"].cells[y]:
# #                 weightMatrixPrior[x,y] = source.postSynapticReceptors[0].weight
# #                 # print("Found a match!")
# # figure()
# # pcolor(weightMatrixPrior)
# # colorbar()
# # title('Prior Weights from Input A to Pyramidals A')
# #
# # weightMatrixPriorBA = zeros([params["popCount"], params["popCount"]])
# # for x in range(len(sim.network.populations["InputB"].cells)):
# #     for y in range(len(sim.network.populations["pyramidalsA"].cells)):
# #         for source in sim.network.populations["InputB"].cells[x].outputs:
# #             if source.target == sim.network.populations["pyramidalsA"].cells[y]:
# #                 weightMatrixPriorBA[x,y] = source.postSynapticReceptors[0].weight
# #
# # figure()
# # pcolor(weightMatrixPriorBA)
# # colorbar()
# # title('Prior Weights from Input B to Pyramidals A')

# sim.run()

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

# with open("TwoColumnPickle_Pharma.jar", "wb") as pickleJar:
#     dill.dump(deepcopy(sim), pickleJar)

with open("TwoColumnPickle_Pharma.jar", "rb") as pickleJar:
    sim = dill.load(pickleJar)

sim.plotColumns()
print("FINISHED!")

# with open("pickle.jar", "wb") as pickleJar:
#     dill.dump(weightMatrixPostBA, pickleJar)

# # Turn off plasticity?
# for x in range(len(sim.network.populations["InputA"].cells)):
#     for y in range(len(sim.network.populations["pyramidalsA"].cells)):
#         for source in sim.network.populations["InputA"].cells[x].outputs:
#             # print(source.target.name, sim.network.populations["pyramidalsA"].cells[y].name)
#             if source.target == sim.network.populations["pyramidalsA"].cells[y]:
#                 source.postSynapticReceptors[0].plasticity = False
#
# # Diffuse Neurotransmitter Paramters
# params["serotoninLevelA"] = 40
# params["serotoninLevelB"] = 10
# params["Somatic5HT2AWeight"] = 100
# params["Somatic5HT2AWeightLTS"] = 80
# params["Somatic5HT1AWeight"] = -80
# params["Axonal5HT2AWeight"] = 0.2
# params["Axonal5HT1AWeight"] = -0.2

# winsound.Beep(440, 666)

