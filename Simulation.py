from SerotoninAVNetwork import *
from scipy import *

class Simulation():
    def __init__(self, params):
        self.tau = 0.1 # In
        self.params = params
        self.tspan = arange(0, self.params["maxTime"], self.tau)

        self.network = SerotoninAVNetwork(self.tau, self, self.params, "Test Network1")

    def run(self):
        for t in self.tspan:
            if t % 10 == 0:
                print("Time: ", t)
            self.network.step()