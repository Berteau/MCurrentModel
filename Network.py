from abc import ABC

class Network(ABC):
    def __init__(self, tau, parentRun, params):
        self.tau = tau
        self.parentRun = parentRun
        self.params = params
        self.populations = []

    def getParentRun(self):
        return self.parentRun

    def step(self):
        # First, update all the neurons
        for population in self.populations:
            population.stepCells()
        # Then, update all the axons
        for population in self.populations:
            population.stepOutputs()





