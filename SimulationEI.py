from EINetwork import *
from scipy import *

class SimulationEI():
    def __init__(self, params):
        self.tau = 0.1 # In
        self.params = params
        self.tspan = arange(0, self.params["maxTime"], self.tau)

        self.network = EINetwork(self.tau, self, self.params, "Test Network1")

    def run(self):
        for t in self.tspan:
            if t % 10 == 0:
                print("Time: ", t)
            self.network.step()

    def plot(self):
        inputAVoltages = [c.vv for c in self.network.populations["Input"].cells]
        pyramidalVoltages = [c.vv for c in self.network.populations["pyramidals"].cells]
        fsVoltages = [c.vv for c in self.network.populations["fastSpikings"].cells]

        figure()
        # title("Column A")
        subplot(3, 1, 1)
        pcolor(inputAVoltages, vmin=-100, vmax=60)
        colorbar()
        title('Input')

        subplot(3, 1, 2)
        pcolor(pyramidalVoltages, vmin=-100, vmax=60)
        colorbar()
        title('Pyramidal Cells')

        subplot(3, 1, 3)
        pcolor(fsVoltages, vmin=-100, vmax=60)
        colorbar()
        title('FS Cells')

        show()

