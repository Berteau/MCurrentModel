from TwoColumnNetworkRestore import *
from scipy import *
import sys

class TwoColumnSimulationRestore():
    def __init__(self, params):
        self.params = params
        self.tau = self.params["tau"]  # In ms
        self.tspan = arange(0, self.params["maxTime"], self.tau)

        self.network = TwoColumnNetworkRestore(self.tau, self, self.params, "Test Network1")

    def run(self):
        for t in self.tspan:
            if t % 10 == 0:
                print("Time: ", t)
            self.network.step()

    def plotColumns(self):
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

        figure()
        # title("Column A")
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
