from abc import ABC

class PostSynapticReceptor(ABC):
    def __init__(self, targetNeuron, weight):
        self.targetNeuron = targetNeuron
        targetNeuron.registerPostSynapticReceptor(self)
        self.weight = weight

    def boutonSpike(self):
        return

    def postSynapticSpikeFeedback(self):
        return
    
    def step(self):
        return

    def injectCurrent(self):
        return

