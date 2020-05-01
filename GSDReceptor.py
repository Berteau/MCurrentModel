from PostSynapticReceptor import *

class InjectedCurrentReceptor(PostSynapticReceptor):
    def __init__(self, targetNeuron, weight, tau):
        self.K = 0
        self.Q = 0
        self.g_SD = 0
        if weight <= 0:
            self.tau_s = tau
            self.tau_r = 8 # ms
            self.tau_f = 30 # ms
        else:
            self.tau_s = tau
            self.tau_r = 2 # ms
            self.tau_f = 10 # ms
            self.tau_r_nmda = 8 # ms
            self.tau_f_nmda = 100 # ms
            self.Q_nmda = 0
            self.g_SD_nmda = 0

        self.targetNeuron = targetNeuron
        targetNeuron.registerPostSynapticReceptor(self)
        self.weight = weight
        self.currentlySpiking = False

    def boutonSpike(self):
        self.currentlySpiking = True
        self.K = 1

    def step(self):
        dQ = ((1-self.Q)*self.K - (self.Q/self.tau_r))/self.tau_s
        dg_SD = ((self.tau_f + self.tau_r)/self.tau_f)*((2/self.tau_r)*(1-self.g_SD)*self.Q-(self.g_SD/self.tau_f))/self.tau_s

        dQ_nmda = ((1 - self.Q_nmda) * self.K - (self.Q / self.tau_r)) / self.tau_s
        dg_SD_nmda = ((self.tau_f + self.tau_r) / self.tau_f) * ((2 / self.tau_r) * (1 - self.g_SD) * self.Q - (self.g_SD / self.tau_f)) / self.tau_s

        self.Q += dQ
        self.g_SD += dg_SD
        self.Q_nmda += dQ_nmda
        self.g_SD_nmda += dg_SD_nmda

    def injectCurrent(self):
        self.targetNeuron.addSynapticTransmission(self.weight*(self.g_SD+self.g_SD_nmda))
        if self.K == 1:
            self.K = 0