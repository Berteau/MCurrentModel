from Simulation import *
import winsound

params = {}
params["pyramidalSelfExcitationWeight"] = 15000
params["PyramidalsToFSWeight"] = 50000
params["FSToPyramidalsWeight"] = -30000
params["PyramidalsToLTSWeight"] = 80000
params["popCount"] = 40
params["serotoninLevelV"] = 10
params["serotoninLevelA"] = 10
params["inputWeightA1"] = 1500
params["InputWeightA2"] = 1500
params["rateA1"] = 20
params["rateA2"] = 3
params["inputWeightV1"] = 1500
params["InputWeightV2"] = 1500
params["InputWeightV3"] = 1500
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
sim.run()
print("FINISHED!")

v1PyramidalVoltages = [sim.network.populations["pyramidalsV1"].cells[i].vv for i in range(len(sim.network.populations["pyramidalsV1"].cells))]
a1PyramidalVoltages = [sim.network.populations["pyramidalsA1"].cells[i].vv for i in range(len(sim.network.populations["pyramidalsA1"].cells))]
v1FSVoltages = [sim.network.populations["fastSpikingsV1"].cells[i].vv for i in range(len(sim.network.populations["fastSpikingsV1"].cells))]
v1LTSVoltages = [sim.network.populations["lowThresholdsV1"].cells[i].vv for i in range(len(sim.network.populations["lowThresholdsV1"].cells))]

figure()
pcolor(v1PyramidalVoltages, vmin=-100, vmax=60)
colorbar()
title('V1 Pyramidal Cells')

figure()
pcolor(v1FSVoltages, vmin=-100, vmax=60)
colorbar()
title('V1 FS Cells')

figure()
pcolor(v1LTSVoltages, vmin=-100, vmax=60)
colorbar()
title('V1 LTS Cells')

figure()
pcolor(a1PyramidalVoltages, vmin=-100, vmax=60)
colorbar()
title('A1 Pyramidal Cells')


# figure()
# plot(v1PyramidalVoltages[0])

show()

winsound.Beep(440, 666)