import dill
import sys
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('pdf')

# arrARate = np.empty()
# arrBRate = np.empty()
# arrABInfluence = np.empty()

inputARate = []
inputBRate = []
aRate = []
bRate = []
aFSRate = []
bFSRate = []
aLTSRate = []
bLTSRate = []
aInfluence = []
bInfluence = []

baSpikeFailures = []

# i = 0
# with open("TwoColumnPickle" + str(i) + ".jar", "rb") as pickleJar:
#     sim = dill.load(pickleJar)
#     np.append(arrARate, sim.network.populations["InputA"].rateRecord, axis=0)
#     np.append(arrBRate, sim.network.populations["InputB"].rateRecord, axis=0)
#
#     baInfluence = sim.network.populations["InputB"].influenceRecord[sim.network.populations["pyramidalsA"]]
#     aaInfluence = sim.network.populations["InputA"].influenceRecord[sim.network.populations["pyramidalsA"]]
#     np.append(arrABInfluence, [baInfluence[x] - aaInfluence[x] for x in range(len(baInfluence))], axis=0)


for i in range(10):
    print(i)
    with open("TwoColumnPickle" + str(i) + "_NP.jar", "rb") as pickleJar:
            sim = dill.load(pickleJar)
            # print(len(sim.network.populations["InputA"].rateRecord))
            inputARate.append(sim.network.populations["InputA"].rateRecord)
            inputBRate.append(sim.network.populations["InputB"].rateRecord)
            aRate.append(sim.network.populations["pyramidalsA"].rateRecord)
            bRate.append(sim.network.populations["pyramidalsB"].rateRecord)
            aFSRate.append(sim.network.populations["fastSpikingsA"].rateRecord)
            bFSRate.append(sim.network.populations["fastSpikingsB"].rateRecord)
            aLTSRate.append(sim.network.populations["lowThresholdsA"].rateRecord)
            bLTSRate.append(sim.network.populations["lowThresholdsB"].rateRecord)
            baInfluence = sim.network.populations["InputB"].influenceRecord[sim.network.populations["pyramidalsA"]]
            aaInfluence = sim.network.populations["InputA"].influenceRecord[sim.network.populations["pyramidalsA"]]
            abInfluence = sim.network.populations["InputA"].influenceRecord[sim.network.populations["pyramidalsB"]]
            bbInfluence = sim.network.populations["InputB"].influenceRecord[sim.network.populations["pyramidalsB"]]
            aInfluence.append([baInfluence[x] - aaInfluence[x] for x in range(len(baInfluence))])
            bInfluence.append([bbInfluence[x] - abInfluence[x] for x in range(len(abInfluence))])
            numax = len(sim.network.populations["InputB"].outboundAxons)
            timeSpan = len(sim.network.populations["InputB"].rateRecord)
            baSpikeFailuresTemp = np.zeros((numax, timeSpan))
            axNum = 0
            for ax in sim.network.populations["InputB"].outboundAxons:
                    baSpikeFailuresTemp[axNum,:] = ax.spikeFailures
                    axNum += 1
            # To Do: sum at each time point, append that as a new line to baSpikeFailures
            baSpikeFailures.append(np.sum(baSpikeFailuresTemp, axis=0))

            # np.append(arrARate, sim.network.populations["InputA"].rateRecord, axis=1)
            # np.append(arrBRate, sim.network.populations["InputB"gggg.rateRecord, axis=1)

            # np.append(arrABInfluence, [baInfluence[x] - aaInfluence[x] for x in range(len(baInfluence))], axis=1)

# print(aRate)

arrInputARate = np.array(inputARate)
arrInputBRate = np.array(inputBRate)
arrARate = np.array(aRate)
arrBRate = np.array(bRate)
arrAFSRate = np.array(aFSRate)
arrBFSRate = np.array(bFSRate)
arrALTSRate = np.array(aLTSRate)
arrBLTSRate = np.array(bLTSRate)
arrAInfluence = np.array(aInfluence)
arrBInfluence = np.array(bInfluence)
arrBASpikeFailures = np.array(baSpikeFailures)

# print(arrARate)

meanInputARate = np.mean(arrInputARate, axis=0)
meanInputBRate =np.mean(arrInputBRate, axis=0)
meanARate = np.mean(arrARate, axis=0)
meanBRate = np.mean(arrBRate, axis=0)
meanAFSRate = np.mean(arrAFSRate, axis=0)
meanBFSRate = np.mean(arrBFSRate, axis=0)
meanALTSRate = np.mean(arrALTSRate, axis=0)
meanBLTSRate = np.mean(arrBLTSRate, axis=0)
meanAInfluence = np.mean(arrAInfluence, axis=0)
meanBInfluence = np.mean(arrBInfluence, axis=0)
meanBASpikeFailures = np.mean(arrBASpikeFailures, axis=0)

sdInputARate = np.std(arrInputARate, axis=0)
sdInputBRate = np.std(arrInputBRate, axis=0)
sdARate = np.std(arrARate, axis=0)
sdBRate = np.std(arrBRate, axis=0)
sdAFSRate = np.std(arrAFSRate, axis=0)
sdBFSRate = np.std(arrBFSRate, axis=0)
sdALTSRate = np.std(arrALTSRate, axis=0)
sdBLTSRate = np.std(arrBLTSRate, axis=0)
sdAInfluence = np.std(arrAInfluence, axis=0)
sdBInfluence = np.std(arrBInfluence, axis=0)
sdBASpikeFailures = np.std(arrBASpikeFailures, axis=0)

denom = sqrt(i)

upperInputARate = meanInputARate + 1.96 * (sdInputARate / denom)
lowerInputARate = meanInputARate - 1.96 * (sdInputARate / denom)
upperInputBRate = meanInputBRate + 1.96 * (sdInputBRate / denom)
lowerInputBRate = meanInputBRate - 1.96 * (sdInputBRate / denom)

upperARate = meanARate + 1.96 * (sdARate / denom)
lowerARate = meanARate - 1.96 * (sdARate / denom)
upperBRate = meanBRate + 1.96 * (sdBRate / denom)
lowerBRate = meanBRate - 1.96 * (sdBRate / denom)

upperAFSRate = meanAFSRate + 1.96 * (sdAFSRate / denom)
lowerAFSRate = meanAFSRate - 1.96 * (sdAFSRate / denom)
upperBFSRate = meanBFSRate + 1.96 * (sdBFSRate / denom)
lowerBFSRate = meanBFSRate - 1.96 * (sdBFSRate / denom)

upperALTSRate = meanALTSRate + 1.96 * (sdALTSRate / denom)
lowerALTSRate = meanALTSRate - 1.96 * (sdALTSRate / denom)
upperBLTSRate = meanBLTSRate + 1.96 * (sdBLTSRate / denom)
lowerBLTSRate = meanBLTSRate - 1.96 * (sdBLTSRate / denom)

upperAInfluence = meanAInfluence + 1.96 * (sdAInfluence / denom)
lowerAInfluence = meanAInfluence - 1.96 * (sdAInfluence / denom)

upperBInfluence = meanBInfluence + 1.96 * (sdBInfluence / denom)
lowerBInfluence = meanBInfluence - 1.96 * (sdBInfluence / denom)

upperBaSpikeFailures = meanBASpikeFailures + 1.96 * (sdBASpikeFailures / denom)
lowerBaSpikeFailures = meanBASpikeFailures - 1.96 * (sdBASpikeFailures / denom)

plt.figure(figsize=(9.37,2.58))
plt.title("Input A Firing Rate")
plt.plot(upperInputARate, color="red")
plt.plot(lowerInputARate, color="red")
plt.plot(meanInputARate)
plt.ylim(0, 2000)
plt.savefig("A_Input_CI_NonPharma.pdf")
plt.savefig("A_Input_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("A Pyramidal Firing Rate")
plt.plot(upperARate, color="red")
plt.plot(lowerARate, color="red")
plt.plot(meanARate)
plt.ylim(0, 2000)
plt.savefig("A_Pyramidal_CI_NonPharma.pdf")
plt.savefig("A_Pyramidal_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("A FS Firing Rate")
plt.plot(upperAFSRate, color="red")
plt.plot(lowerAFSRate, color="red")
plt.plot(meanAFSRate)
plt.ylim(0, 2000)
plt.savefig("A_FS_CI_NonPharma.pdf")
plt.savefig("A_FS_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("A LTS Firing Rate")
plt.plot(upperALTSRate, color="red")
plt.plot(lowerALTSRate, color="red")
plt.plot(meanALTSRate)
plt.ylim(0, 2000)
plt.savefig("A_LTS_CI_NonPharma.pdf")
plt.savefig("A_LTS_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("Input B Firing Rate")
plt.plot(upperInputBRate, color="red")
plt.plot(lowerInputBRate, color="red")
plt.plot(meanInputBRate)
plt.ylim(0, 2000)
plt.savefig("B_Input_CI_NonPharma.pdf")
plt.savefig("B_Input_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("B Pyramidal Firing Rate")
plt.plot(upperBRate, color="red")
plt.plot(lowerBRate, color="red")
plt.plot(meanBRate)
plt.ylim(0, 2000)
plt.savefig("B_Pyramidal_CI_NonPharma.pdf")
plt.savefig("B_Pyramidal_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("B FS Firing Rate")
plt.plot(upperBFSRate, color="red")
plt.plot(lowerBFSRate, color="red")
plt.plot(meanBFSRate)
plt.ylim(0, 2000)
plt.savefig("B_FS_CI_NonPharma.pdf")
plt.savefig("B_FS_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("B LTS Firing Rate")
plt.plot(upperBLTSRate, color="red")
plt.plot(lowerBLTSRate, color="red")
plt.plot(meanBLTSRate)
plt.ylim(0, 2000)
plt.savefig("B_LTS_CI_NonPharma.pdf")
plt.savefig("B_LTS_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("Relative influence of Inputs A (negative) and B (positive) on A Pyramidal Firing")
plt.plot(upperAInfluence, color="red")
plt.plot(lowerAInfluence, color="red")
plt.plot(meanAInfluence)
plt.savefig("A_Influence_CI_NonPharma.pdf")
plt.savefig("A_Influence_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("Relative influence of Inputs A (negative) and B (positive) on B Pyramidal Firing")
plt.plot(upperBInfluence, color="red")
plt.plot(lowerBInfluence, color="red")
plt.plot(meanBInfluence)
plt.savefig("B_Influence_CI_NonPharma.pdf")
plt.savefig("B_Influence_CI_NonPharma.png")

plt.figure(figsize=(9.37,2.58))
plt.title("Number of spike transmission failures from Input B to A Pyramidal Cells")
plt.plot(upperBaSpikeFailures, color="red")
plt.plot(lowerBaSpikeFailures, color="red")
plt.plot(meanBASpikeFailures)
plt.savefig("BA_SpikeFailures_CI_NonPharma.pdf")
plt.savefig("BA_SpikeFailures_CI_NonPharma.png")

# if sys.argv[0] is not '':
#     plt.show()

# inputAVoltages = [c.vv for c in self.network.populations["InputA"].cells]
#         inputBVoltages = [c.vv for c in self.network.populations["InputB"].cells]

#         aPyramidalVoltages = [c.vv for c in self.network.populations["pyramidalsA"].cells]
#         aFSVoltages = [c.vv for c in self.network.populations["fastSpikingsA"].cells]
#         aLTSVoltages = [c.vv for c in self.network.populations["lowThresholdsA"].cells]

#         bPyramidalVoltages = [c.vv for c in self.network.populations["pyramidalsB"].cells]
#         bFSVoltages = [c.vv for c in self.network.populations["fastSpikingsB"].cells]
#         bLTSVoltages = [c.vv for c in self.network.populations["lowThresholdsB"].cells]

#         aPyramidalSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]
#         bPyramidalSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsB"].cells]

#         print("A Spikes: ", sum(aPyramidalSpikes))
#         print("B Spikes: ", sum(bPyramidalSpikes))

#         # Sliding Windowed Spike Rates
#         figure()
#         subplot(4,1,1)
#         plot(self.network.populations["InputA"].rateRecord)
#         title('A Input')

#         subplot(4, 1, 2)
#         plot(self.network.populations["pyramidalsA"].rateRecord)
#         title('A Pyramidal Cells')

#         subplot(4, 1, 3)
#         plot(self.network.populations["fastSpikingsA"].rateRecord)
#         title('A FS Cells')

#         subplot(4, 1, 4)
#         plot(self.network.populations["lowThresholdsA"].rateRecord)
#         title('A LTS Cells')

#         figure()
#         subplot(4, 1, 1)
#         plot(self.network.populations["InputB"].rateRecord)
#         title('B Input')

#         subplot(4, 1, 2)
#         plot(self.network.populations["pyramidalsB"].rateRecord)
#         title('B Pyramidal Cells')

#         subplot(4, 1, 3)
#         plot(self.network.populations["fastSpikingsB"].rateRecord)
#         title('B FS Cells')

#         subplot(4, 1, 4)
#         plot(self.network.populations["lowThresholdsB"].rateRecord)
#         title('B LTS Cells')

#         figure()
#         # title("Column A")
#         subplot(4, 1, 1)
#         pcolor(inputAVoltages, vmin=-100, vmax=60)
#         colorbar()
#         title('A Input')

#         subplot(4, 1, 2)
#         pcolor(aPyramidalVoltages, vmin=-100, vmax=60)
#         plot(self.network.populations["pyramidalsA"].rateRecord, color='White')
#         colorbar()
#         title('A Pyramidal Cells')

#         subplot(4, 1, 3)
#         pcolor(aFSVoltages, vmin=-100, vmax=60)
#         plot(self.network.populations["fastSpikingsA"].rateRecord, color='White')
#         colorbar()
#         title('A FS Cells')

#         subplot(4, 1, 4)
#         pcolor(aLTSVoltages, vmin=-100, vmax=60)
#         plot(self.network.populations["lowThresholdsA"].rateRecord, color='White')
#         colorbar()
#         title('A LTS Cells')

#         figure()
#         subplot(4, 1, 1)
#         pcolor(inputBVoltages, vmin=-100, vmax=60)
#         plot(self.network.populations["InputB"].rateRecord, color='White')
#         colorbar()
#         title('B Input')

#         subplot(4, 1, 2)
#         pcolor(bPyramidalVoltages, vmin=-100, vmax=60)
#         plot(self.network.populations["pyramidalsB"].rateRecord, color='White')
#         colorbar()
#         title('B Pyramidal Cells')

#         subplot(4, 1, 3)
#         pcolor(bFSVoltages, vmin=-100, vmax=60)
#         plot(self.network.populations["fastSpikingsB"].rateRecord, color='White')
#         colorbar()
#         title('B FS Cells')

#         subplot(4, 1, 4)
#         pcolor(bLTSVoltages, vmin=-100, vmax=60)
#         plot(self.network.populations["lowThresholdsB"].rateRecord, color='White')
#         colorbar()
#         title('B LTS Cells')

