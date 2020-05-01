from TwoColumnNetwork import *
from scipy import *
import sys

class TwoColumnSimulationPharma():
    def __init__(self, params):
        self.params = params
        self.tau = self.params["tau"]  # In ms
        self.tspan1 = arange(0, self.params["maxTime"], self.tau)
        self.tspan2 = arange(self.params["maxTime"], self.params["maxTime"]*2, self.tau)
        self.tspan3 = arange(self.params["maxTime"]*2, self.params["maxTime"]*3, self.tau)
        self.tspan4 = arange(self.params["maxTime"]*3, self.params["maxTime"]*4, self.tau)
        self.network = TwoColumnNetwork(self.tau, self, self.params, "Test Network1")

    def run(self):
        # First Portion, full input
        print("Phase One")
        for t in self.tspan1:
            # if t % 10 == 0:
            # print("Phase 1 Time: ", t)
            self.network.step()

        self.aEndOfFirstPortionSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]

        # Second portion, Serotonin and Plasticity
        print("Phase Two")
        weightMatrixPriorBA = zeros([self.params["popCount"], self.params["popCount"]])
        for x in range(len(self.network.populations["InputB"].cells)):
            for y in range(len(self.network.populations["pyramidalsA"].cells)):
                for source in self.network.populations["InputB"].cells[x].outputs:
                    if source.target == self.network.populations["pyramidalsA"].cells[y]:
                        weightMatrixPriorBA[x, y] = source.postSynapticReceptors[0].weight

        figure()
        pcolor(weightMatrixPriorBA)
        colorbar()
        title('Prior Weights from Input B to Pyramidals A')

        # Increase 5HT
        transmittersA = {}
        transmittersA["5HT2A"] = 11
        transmittersA["5HT1A"] = 10
        self.network.setSerotoninA(transmittersA)

        transmittersB = {}
        transmittersB["5HT2A"] = 11
        transmittersB["5HT1A"] = 10
        self.network.setSerotoninB(transmittersB)

        # Turn on plasticity
        for axon in self.network.populations["InputB"].outboundAxons:
            if axon.weight > 0:
                axon.postSynapticReceptors[0].plasticity = True
                axon.postSynapticReceptors[0].c_p = 0.4
        for axon in self.network.populations["InputA"].outboundAxons:
            if axon.weight > 0:
                axon.postSynapticReceptors[0].plasticity = True
                axon.postSynapticReceptors[0].c_p = 0.43

        # for x in range(len(self.network.populations["InputB"].cells)):
        #     for y in range(len(self.network.populations["pyramidalsA"].cells)):
        #         for source in self.network.populations["InputB"].cells[x].outputs:
        #             # print(source.target.name, sim.network.populations["pyramidalsA"].cells[y].name)
        #             if source.target == self.network.populations["pyramidalsA"].cells[y]:
        #                 source.postSynapticReceptors[0].plasticity = True

        for t in self.tspan3:
            # if t % 10 == 0:
            # print("Phase 3 Time: ", t)
            self.network.step()

        weightMatrixPostBA = zeros([self.params["popCount"], self.params["popCount"]])
        for x in range(len(self.network.populations["InputB"].cells)):
            for y in range(len(self.network.populations["pyramidalsA"].cells)):
                for source in self.network.populations["InputB"].cells[x].outputs:
                    if source.target == self.network.populations["pyramidalsA"].cells[y]:
                        weightMatrixPostBA[x,y] = source.postSynapticReceptors[0].weight

        figure()
        pcolor(weightMatrixPostBA)
        colorbar()
        title('Posterior Weights from Input B to Pyramidals A')

        self.aEndOfSecondPortionSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]

        # Third portion, remapped functionality.
        print("Phase Three")
        # Turn off plasticity
        for axon in self.network.populations["InputB"].outboundAxons:
            if axon.weight > 0:
                axon.postSynapticReceptors[0].plasticity = False
                axon.postSynapticReceptors[0].c_p = 0.0
        for axon in self.network.populations["InputA"].outboundAxons:
            if axon.weight > 0:
                axon.postSynapticReceptors[0].plasticity = False
                axon.postSynapticReceptors[0].c_p = 0.0
        # for x in range(len(self.network.populations["InputA"].cells)):
        #     for y in range(len(self.network.populations["pyramidalsA"].cells)):
        #         for source in self.network.populations["InputA"].cells[x].outputs:
        #             # print(source.target.name, sim.network.populations["pyramidalsA"].cells[y].name)
        #             if source.target == self.network.populations["pyramidalsA"].cells[y]:
        #                 source.postSynapticReceptors[0].plasticity = False

        transmittersA = {}
        transmittersA["5HT2A"] = 10
        transmittersA["5HT1A"] = 10
        self.network.setSerotoninA(transmittersA)
        transmittersB = {}
        transmittersB["5HT2A"] = 10
        transmittersB["5HT1A"] = 10
        self.network.setSerotoninB(transmittersB)


        for t in self.tspan4:
            # if t % 10 == 0:
            # print("Phase 4 Time: ", t)
            self.network.step()

        # for t in self.tspan:
        #     # if t % 10 == 0:
        #     print("Time: ", t)
        #     self.network.step()

        weightMatrixEndBA = zeros([self.params["popCount"], self.params["popCount"]])
        for x in range(len(self.network.populations["InputB"].cells)):
            for y in range(len(self.network.populations["pyramidalsA"].cells)):
                for source in self.network.populations["InputB"].cells[x].outputs:
                    if source.target == self.network.populations["pyramidalsA"].cells[y]:
                        weightMatrixEndBA[x, y] = source.postSynapticReceptors[0].weight

        figure()
        pcolor(weightMatrixEndBA)
        colorbar()
        title('Final Weights from Input B to Pyramidals A')
        self.aEndOfThirdPortionSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]

    def plotColumns(self):

        baInfluence = self.network.populations["InputB"].influenceRecord[self.network.populations["pyramidalsA"]]
        aaInfluence = self.network.populations["InputA"].influenceRecord[self.network.populations["pyramidalsA"]]
        influenceDiff = [baInfluence[x] - aaInfluence[x] for x in range(len(baInfluence))]

        inputAVoltages = [c.vv for c in self.network.populations["InputA"].cells]
        inputBVoltages = [c.vv for c in self.network.populations["InputB"].cells]

        aPyramidalVoltages = [c.vv for c in self.network.populations["pyramidalsA"].cells]
        aFSVoltages = [c.vv for c in self.network.populations["fastSpikingsA"].cells]
        aLTSVoltages = [c.vv for c in self.network.populations["lowThresholdsA"].cells]

        bPyramidalVoltages = [c.vv for c in self.network.populations["pyramidalsB"].cells]
        bFSVoltages = [c.vv for c in self.network.populations["fastSpikingsB"].cells]
        bLTSVoltages = [c.vv for c in self.network.populations["lowThresholdsB"].cells]

        aPyramidalSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]
        bPyramidalSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsB"].cells]

        print("A Spikes: ", sum(aPyramidalSpikes))
        print("B Spikes: ", sum(bPyramidalSpikes))

        # Sliding Windowed Spike Rates
        figure()
        subplot(5, 1, 1)
        plot(self.network.populations["InputA"].rateRecord)
        title('A Input')

        subplot(5, 1, 2)
        plot(self.network.populations["pyramidalsA"].rateRecord)
        title('A Pyramidal Cells')

        subplot(5, 1, 3)
        plot(self.network.populations["fastSpikingsA"].rateRecord)
        title('A FS Cells')

        subplot(5, 1, 4)
        plot(self.network.populations["lowThresholdsA"].rateRecord)
        title('A LTS Cells')

        subplot(5, 1, 5)
        plot(influenceDiff)
        title('A/B Influence Balance')

        # B Rates
        figure()
        subplot(4, 1, 1)
        plot(self.network.populations["InputB"].rateRecord)
        title('B Input')

        subplot(4, 1, 2)
        plot(self.network.populations["pyramidalsB"].rateRecord)
        title('B Pyramidal Cells')

        subplot(4, 1, 3)
        plot(self.network.populations["fastSpikingsB"].rateRecord)
        title('B FS Cells')

        subplot(4, 1, 4)
        plot(self.network.populations["lowThresholdsB"].rateRecord)
        title('B LTS Cells')

        # A Rasters
        figure()
        subplot(4, 1, 1)
        pcolor(inputAVoltages, vmin=-100, vmax=60)
        colorbar()
        title('A Input')

        subplot(4, 1, 2)
        pcolor(aPyramidalVoltages, vmin=-100, vmax=60)
        colorbar()
        title('A Pyramidal Cells')

        subplot(4, 1, 3)
        pcolor(aFSVoltages, vmin=-100, vmax=60)
        colorbar()
        title('A FS Cells')

        subplot(4, 1, 4)
        pcolor(aLTSVoltages, vmin=-100, vmax=60)
        colorbar()
        title('A LTS Cells')

        # B Rasters
        figure()
        subplot(4, 1, 1)
        pcolor(inputBVoltages, vmin=-100, vmax=60)
        colorbar()
        title('B Input')

        subplot(4, 1, 2)
        pcolor(bPyramidalVoltages, vmin=-100, vmax=60)
        colorbar()
        title('B Pyramidal Cells')

        subplot(4, 1, 3)
        pcolor(bFSVoltages, vmin=-100, vmax=60)
        colorbar()
        title('B FS Cells')

        subplot(4, 1, 4)
        pcolor(bLTSVoltages, vmin=-100, vmax=60)
        colorbar()
        title('B LTS Cells')

        if sys.argv[0] is not '':
            show()
