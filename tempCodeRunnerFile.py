from TwoColumnSimulation import *
import winsound

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
params["PyramidalsToLTSWeight"] = 500
params["LTStoFSWeight"] = -500
params["LTStoPyramidalsWeight"] = -1000

# Input Paramters
params["inputWeightA"] = 500
params["inputWeightB"] = 500
params["rateA"] = 0
params["rateB"] = 10
params["inputWeightAB"] = 250
params["crossModalABLikelihood"] = 0.3
params["inputWeightBA"] = 250 #7500
params["crossModalBALikelihood"] = 0.3

# Diffuse Neurotransmitter Paramters
params["serotoninLevelA"] = 40
params["serotoninLevelB"] = 10
params["Somatic5HT2AWeight"] = 100
params["Somatic5HT2AWeightLTS"] = 80
params["Somatic5HT1AWeight"] = -80
params["Axonal5HT2AWeight"] = 0.2
params["Axonal5HT1AWeight"] = -0.2

sim = TwoColumnSimulation(params)
weightMatrixPrior = zeros([params["popCount"], params["popCount"]])
for x in range(len(sim.network.populations["InputA"].cells)):
    for y in range(len(sim.network.populations["pyramidalsA"].cells)):
        for source in sim.network.populations["InputA"].cells[x].outputs:
            if source.target == sim.network.populations["InputA"].cells[y]:
                weightMatrixPrior[x,y] = source.weight
                print("Found a match!")
figure()
pcolor(weightMatrixPrior)
colorbar()
title('Prior Weights from Input A to Pyramidals A')

weightMatrixPriorBA = zeros([params["popCount"], params["popCount"]])
for x in range(len(sim.network.populations["InputB"].cells)):
    for y in range(len(sim.network.populations["pyramidalsA"].cells)):
        for source in sim.network.populations["InputB"].cells[x].outputs:
            if source.target == sim.network.populations["InputB"].cells[y]:
                weightMatrixPriorBA[x,y] = source.weight

figure()
pcolor(weightMatrixPrior)
colorbar()
title('Prior Weights from Input B to Pyramidals A')

sim.run()

weightMatrixPost = zeros([params["popCount"], params["popCount"]])
for x in range(len(sim.network.populations["InputA"].cells)):
    for y in range(len(sim.network.populations["pyramidalsA"].cells)):
        for source in sim.network.populations["InputA"].cells[x].outputs:
            if source.target == sim.network.populations["InputA"].cells[y]:
                weightMatrixPost[x,y] = source.weight

figure()
pcolor(weightMatrixPost)
colorbar()
title('Posterior Weights from Input A to Pyramidals A')

weightMatrixPostBA = zeros([params["popCount"], params["popCount"]])
for x in range(len(sim.network.populations["InputB"].cells)):
    for y in range(len(sim.network.populations["pyramidalsA"].cells)):
        for source in sim.network.populations["InputB"].cells[x].outputs:
            if source.target == sim.network.populations["InputB"].cells[y]:
                weightMatrixPostBA[x,y] = source.weight

figure()
pcolor(weightMatrixPostBA)
colorbar()
title('Posterior Weights from Input B to Pyramidals A')


figure()
sim.plotColumns()
print("FINISHED!")


# winsound.Beep(440, 666)

