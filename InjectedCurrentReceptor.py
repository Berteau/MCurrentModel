from PostSynapticReceptor import *


class InjectedCurrentReceptor(PostSynapticReceptor):
    def __init__(self, targetNeuron, weight):
        self.targetNeuron = targetNeuron
        targetNeuron.registerPostSynapticReceptor(self)
        self.weight = weight
        self.currentlySpiking = False

    def boutonSpike(self):
        self.currentlySpiking = True

    def step(self):
        return

    def injectCurrent(self):
        if self.currentlySpiking:
            self.targetNeuron.addSynapticTransmission(self.weight)
            self.currentlySpiking = False
