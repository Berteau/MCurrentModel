from PostSynapticReceptor import *

class GABAGSDReceptor(PostSynapticReceptor):
    def __init__(self, targetNeuron, weight, tau):
        self.K = 0
        self.Q = 0
        self.g_SD = 0
        self.tau_s = 1
        self.tau_r = 8 # ms
        self.tau_f = 30 # ms

        self.targetNeuron = targetNeuron
        targetNeuron.registerPostSynapticReceptor(self)
        self.weight = weight
        self.currentlySpiking = False

    def boutonSpike(self):
        self.currentlySpiking = True
        self.K = 1

    def step(self):
        dQ = ((1-self.Q)*self.K - (self.Q/self.tau_r))*self.tau_s
        dg_SD = ((self.tau_f + self.tau_r)/self.tau_f)*((2/self.tau_r)*(1-self.g_SD)*self.Q-(self.g_SD/self.tau_f))*self.tau_s
        self.Q += dQ
        self.g_SD += dg_SD

    def injectCurrent(self):
        self.targetNeuron.addSynapticTransmission(self.weight*self.g_SD)
        if self.K == 1:
            self.K = 0