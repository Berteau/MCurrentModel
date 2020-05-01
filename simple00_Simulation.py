from simple00 import *
from scipy import *
from pylab import *
import sys

class simple00_Simulation():
    def __init__(self, params):
        self.params = params
        self.tau = self.params["tau"]  # In ms
        self.tspan1 = arange(0, self.params["maxTime"], self.tau)
        self.tspan2 = arange(self.params["maxTime"], self.params["maxTime"]*2, self.tau)
        self.tspan3 = arange(self.params["maxTime"]*2, self.params["maxTime"]*3, self.tau)
        self.tspan4 = arange(self.params["maxTime"]*3, self.params["maxTime"]*4, self.tau)
        
        self.network = simple00(self.tau, self, self.params, "Test Network1")

    def run(self):
        # First Portion, full input
        print("Phase One")
        for t in self.tspan1:
            # if t % 10 == 0:
            # print("Phase 1 Time: ", t)
            self.network.step()

        self.aEndOfFirstPortionSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]

        # Second portion, sensory deprivation
        print("Phase Two")
        for inputCell in self.network.populations["InputA"].cells:
            inputCell.poissonLambda = 0
        for t in self.tspan2:
            # if t % 10 == 0:
            # print("Phase 2 Time: ", t)
            self.network.step()

        self.aEndOfSecondPortionSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]

 

    def plotColumns(self):
       
        inputAVoltages = [c.vv for c in self.network.populations["InputA"].cells]
    
        aPyramidalVoltages = [c.vv for c in self.network.populations["pyramidalsA"].cells]

        aPyramidalSpikes = [len(c.spikeRecord) for c in self.network.populations["pyramidalsA"].cells]

        print("A Spikes: ", sum(aPyramidalSpikes))

        # Sliding Windowed Spike Rates
        figure()
        subplot(2, 1, 1)
        plot(self.network.populations["InputA"].rateRecord)
        title('A Input')

        subplot(2, 1, 2)
        plot(self.network.populations["pyramidalsA"].rateRecord)
        title('A Pyramidal Cells')


        # A Rasters
        figure()
        subplot(2, 1, 1)
        pcolor(inputAVoltages, vmin=-100, vmax=60)
        colorbar()
        title('A Input')

        subplot(2, 1, 2)
        pcolor(aPyramidalVoltages, vmin=-100, vmax=60)
        colorbar()
        title('A Pyramidal Cells')
        if sys.argv[0] is not '':
            show()
