from PostSynapticReceptor import *


class InjectedCurrentReceptor(PostSynapticReceptor):
    def __init__(self, targetNeuron, weight):
        self.targetNeuron = targetNeuron
        targetNeuron.registerPostSynapticReceptor(self)
        self.weight = weight
        self.currentlySpiking = False

    def boutonSpike(self):
        self.currentlySpiking = True
        self.lastSpikeTime = -100
                
                
        dTracedt = (-trace + spikes[time-1]) *tauX
    #    print vOld[0] - LTPThreshold, ", ", uBarLTP - LTPThreshold, "(", uBarLTP, v, ")"
        dMyelindt = -Altd * spikes[time-1] * max(uBarLTD - LTDThreshold, 0) + Altp * trace * max(oV - LTPThreshold, 0) * max(uBarLTP - LTPThreshold, 0) # vOld[0] in place of oV     
        dMyelindt = dMyelindt * myelinTimeConstant

    def injectCurrent(self):
        if self.currentlySpiking:
            self.targetNeuron.addSynapticTransmission(self.weight)
            self.currentlySpiking = False
