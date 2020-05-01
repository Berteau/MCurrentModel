from PostSynapticReceptor import *
from math import exp

class GlutGSDReceptor(PostSynapticReceptor):
    def __init__(self, targetNeuron, weight, tau):
        self.K = 0
        self.Q = 0
        self.Q_nmda = 0
        self.g_SD = 0
        self.g_SD_nmda = 0
        self.tau = tau
        self.time = 0

        self.C_d = 0.01

        self.tau_s = 1
        self.tau_r = 2 # ms
        self.tau_f = 10 # ms
        self.tau_r_nmda = 8 # ms
        self.tau_f_nmda = 100 # ms

        self.plasticity = True

        # Plasticity Parameters
        self.tau_p = 30 #ms
        self.c_p = 2
        self.th = 0.9 #presumably microAmps?

        self.lastSpikeTime = 1000000
        self.lastTargetSpikeTime = 1000000

        self.targetNeuron = targetNeuron
        targetNeuron.registerPostSynapticReceptor(self)
        self.weight = weight
        # print("Weight:", self.weight)
        self.currentlySpiking = False

        self.injectionRecord = []

    def boutonSpike(self):
        # print("BoutonSpike Called, weight:", self.weight)
        self.currentlySpiking = True
        self.K = 1
        self.lastSpikeTime = self.time

        #Implement LTD
        if self.plasticity:
            self.weight -= self.C_d
            # # Are we following a target spike?
            # candidateSpikes = [spikeTime[0] for spikeTime in self.targetNeuron.spikeRecord if self.time - spikeTime[0] < 40]
            # for spikeTime in candidateSpikes:
            #     self.weight += self.c_p * halfRectify(self.g_SD_nmda - self.th)*exp((spikeTime - self.time)/self.tau_p)
            if self.weight <= 0:
                self.weight = 0
                self.plasticity = False
        return
        
    def postSynapticSpikeFeedback(self):
        # If we spiked in the last 40 ms, modify the weight
        # print(self.time, self.lastSpikeTime)
        if 0 < self.time - self.lastSpikeTime <= 40:
            # if self.c_p * halfRectify(self.g_SD_nmda - self.th)*exp((self.lastSpikeTime - self.time)/self.tau_p) > 0:
            #     print("LTP:", self.c_p * halfRectify(self.g_SD_nmda - self.th)*exp((self.lastSpikeTime - self.time)/self.tau_p))
            self.weight += self.c_p * halfRectify(self.g_SD_nmda - self.th)*exp((self.lastSpikeTime - self.time)/self.tau_p)
        return

    def step(self):
        self.time += self.tau       
        dQ = ((1-self.Q)*self.K - (self.Q/self.tau_r))*self.tau_s
        dg_SD = ((self.tau_f + self.tau_r)/self.tau_f)*((2/self.tau_r)*(1-self.g_SD)*self.Q-(self.g_SD/self.tau_f))*self.tau_s
        dQ_nmda = ((1 - self.Q_nmda) * self.K - (self.Q_nmda / self.tau_r_nmda)) * self.tau_s
        dg_SD_nmda = ((self.tau_f_nmda + self.tau_r_nmda) / self.tau_f_nmda) * ((2 / self.tau_r_nmda) * (1 - self.g_SD_nmda) * self.Q_nmda - (self.g_SD_nmda / self.tau_f_nmda)) * self.tau_s
        self.Q += dQ
        self.g_SD += dg_SD
        self.Q_nmda += dQ_nmda
        self.g_SD_nmda += dg_SD_nmda
        if self.currentlySpiking:
            self.currentlySpiking = False
            self.K = 0 

    def injectCurrent(self):
        self.injectionRecord.append(self.weight*(self.g_SD + self.g_SD_nmda))
        self.targetNeuron.addSynapticTransmission(self.weight*(self.g_SD + self.g_SD_nmda))

def halfRectify(value):
    if value <= 0:
        return 0
    else:
        return value