from SimulationEI import *
import winsound

params = {}
# Time Parameters
params["maxTime"] = 100
params["tau"] = 0.1

# Population and Connection Parameters
params["popCount"] = 2
params["pyramidalSelfExcitationWeight"] = 100
params["pyramidalToFSWeight"] = 0
params["fsToPyramidalWeight"] = 0

# Input Paramters
params["inputWeightA"] = 300
params["rateA"] = 20

sim = SimulationEI(params)
sim.run()
print("Final Weight:", sim.network.populations["Input"].cells[0].outputs[0].postSynapticReceptors[0].weight)
figure()
plot(sim.network.populations["Input"].cells[0].outputs[0].postSynapticReceptors[0].injectionRecord)
plot(sim.network.populations["Input"].cells[0].outputs[0].postSynapticReceptors[0].targetNeuron.ii)
plot(sim.network.populations["Input"].cells[0].outputs[0].postSynapticReceptors[0].targetNeuron.vv)
sim.plot()
print("FINISHED!")

# winsound.Beep(440, 666)

