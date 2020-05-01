from PostSynapticReceptor import *

class ConductanceBasedPostSynapticReceptor(PostSynapticReceptor):
    def __init__(self, targetNeuron, weight):
        self.targetNeuron = targetNeuron
        targetNeuron.registerPostSynapticReceptor(self)
        self.weight = weight

        ## Synaptic Parameters and variables
        self.I_syn = 0
        self.I_ampa = 0
        self.I_gaba = 0
        self.I_nmda = 0

        self.s_ampa = 0
        self.s_gaba = 0
        self.s_nmda = 0
        self.s_nmda_x = 0
        self.s_ampa_dt = 0
        self.s_gaba_dt = 0
        self.s_nmda_dt = 0
        self.s_nmda_x_dt = 0

        self.g_ampa = 7.5e-3
        self.g_gaba = 7.5e-3
        self.g_nmda = 2e-3
        self.mg = 1e-3
        self.tau_ampa = 2  # ms
        self.tau_gaba = 10  # ms
        self.tau_nmda_rise = 2  # ms
        self.tau_nmda_decay = 100  # ms
        self.alpha = 0.5  # ms
        self.h = 1.0  # Percent Na channels open

        self.transmission = False


    def boutonSpike(self):
        self.transmission = True
        return

    def injectCurrent(self):
            if self.transmission and self.weight > 0:  # Glutamatergic input
                self.s_ampa_dt = -(self.s_ampa / self.tau_ampa) + self.weight
                self.s_gaba_dt = -(self.s_gaba / self.tau_gaba)
                self.s_nmda_x_dt = -(self.s_nmda_x / self.tau_nmda_rise) + self.weight
                self.s_nmda_dt = -(self.s_nmda / self.tau_nmda_decay) + self.alpha * self.s_nmda_x * (1 - self.s_nmda)
            elif self.transmission and self.weight < 0:  # Gabaergic input
                self.s_ampa_dt = -self.s_ampa / self.tau_ampa
                self.s_gaba_dt = -(self.s_gaba / self.tau_gaba) + self.weight
                self.s_nmda_x_dt = -(self.s_nmda_x / self.tau_nmda_rise)
                self.s_nmda_dt = -(self.s_nmda / self.tau_nmda_decay) + self.alpha * self.s_nmda_x * (1 - self.s_nmda)
            else:  # No input
                self.s_ampa_dt = -self.s_ampa / self.tau_ampa
                self.s_gaba_dt = -self.s_gaba / self.tau_gaba
                self.s_nmda_x_dt = -(self.s_nmda_x / self.tau_nmda_rise)
                self.s_nmda_dt = -(self.s_nmda / self.tau_nmda_decay) + self.alpha * self.s_nmda_x * (1 - self.s_nmda)

            # Now update the actual value of the dynamic terms
            self.s_ampa = self.s_ampa + self.s_ampa_dt * self.tau
            self.s_gaba = self.s_gaba + self.s_gaba_dt * self.tau
            if self.debug:
                print("Synapse", a, ":", self.s_gaba[a])
            self.s_nmda[a] = self.s_nmda[a] + self.s_nmda_dt[a] * self.tau
            self.s_nmda_x[a] = self.s_nmda_x[a] + self.s_nmda_x_dt[a] * self.tau

        self.I_ampa = sum([self.weights[a] * self.s_ampa[a] for a in range(len(self.inputs))])
        self.I_gaba = sum([-self.weights[a] * self.s_gaba[a] for a in range(len(self.inputs))])
        self.I_nmda = sum([self.weights[a] * self.s_nmda[a] for a in range(len(self.inputs))])
        self.ampaOpen = self.I_ampa
        self.gabaOpen = self.I_gaba
        self.nmdaOpen = self.I_nmda
        self.I_ampa = self.I_ampa * (self.g_ampa * (self.v - self.v_Na))
        self.I_gaba = self.I_gaba * (self.g_gaba * (self.v - self.v_shunt))
        self.I_nmda = self.I_nmda * ((self.g_nmda * (self.v - self.v_Na)) / 1 + self.mg * exp((-0.062 * self.v) / 3.57))
        self.I_syn = self.I_ampa + self.I_gaba + self.I_nmda

        return

