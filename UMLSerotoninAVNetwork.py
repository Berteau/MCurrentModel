from Simulation import *
import winsound


params = {}
params["pyramidalSelfExcitationWeight"] = 15000
params["PyramidalsToFSWeight"] = 40000
params["FSToPyramidalsWeight"] = -40000
params["PyramidalsToLTSWeight"] = 50000
params["popCount"] = 40
params["serotoninLevelV"] = 10
params["serotoninLevelA"] = 10
params["inputWeightA1"] = 1000
params["InputWeightA2"] = 1000
params["rateA1"] = 20
params["rateA2"] = 3
params["inputWeightV1"] = 1000
params["InputWeightV2"] = 1000
params["InputWeightV3"] = 1000
params["rateV1"] = 20
params["rateV2"] = 1
params["rateV3"] = 1
params["inputWeightVA"] = 800
params["inputWeightAV"] = 800
params["Somatic5HT2AWeight"] = 20
params["Somatic5HT1AWeight"] = -5
params["Axonal5HT2AWeight"] = 0.2
params["Axonal5HT1AWeight"] = -0.2

sim = Simulation(params)


# sim.run()
# print("FINISHED!")
#
# pyramidalVoltages = [sim.network.populations["pyramidalsV1"].cells[i].vv for i in range(len(sim.network.populations["pyramidalsV1"].cells))]
# fsVoltages = [sim.network.populations["fastSpikingsV1"].cells[i].vv for i in range(len(sim.network.populations["fastSpikingsV1"].cells))]
#
# figure()
# pcolor(pyramidalVoltages, vmin=-100, vmax=60)
# colorbar()
# title('Pyramidal Cells')
#
# figure()
# plot(pyramidalVoltages[0])
#
# show()
#
# winsound.Beep(440, 666)