from __future__ import division

'''
Created on Sep 16, 2011

L2/3 Pyramidal Cell from Mainen and Sejnowski 1996

@author: stefan
'''

from pylab import *
from scipy import *
from numpy import *

class MonopoleNeuron:
    def __init__(self,inputs=None, outputs=None, externalInput = 0.0, position = [0.0, 0.0, 0.0], debug=False, tau=0.1, surfaceArea=1.0, tspan=range((10000)), inactivating=True):
        self.debug = debug
        self.tspan = tspan
        self.inactivating = inactivating
        
        # Constants

        # Membrane Capacitance
        self.C = 0.75  # muF/cm^-2

        # Reversal Potentials
        self.vLeak = -70 #mV
        self.vK = -90 #mV
        self.vNa = -50 #mV
        self.vCa = 140 #mV

        # Maximum Conductances
        self.gNa = 20 #mS/cm^2
        self.gCa = 0.3 #mS/cm^2
        self.gKCa = 3 #mS/cm^2
        self.gKv = 0.1 #mS/cm^2
        self.gLeak = 1.5 #mS/cm^2, still from Prescott, need to find a better source?
        self.gM = 4.0


        ## M-Current

        self.vdapt = self.vK       #mV (for potassium m-current adaptation)
        self.beta = -35          #mV  (m-current threshold, seems different from what I've read in the past)      
        self.gamma = 5           #mV
        self.alpha = 0.005       #unitless
        self.phi = 0.15          #unitless
        self.v_1 = -1.2          #mV
        self.v_2 = 23            #mV
        self.v_3 = -2            #mV
        self.v_4 = 21            #mV


        if inputs is None:
            self.inputs = []
            self.weights = []
        else:
            self.inputs = inputs
            self.weights = [inputs[a].weight for a in self.inputs]
        if outputs is None:
            self.outputs = []
        else:
            self.outputs = outputs
        self.tau = tau
        self.position = position
        self.externalInput = externalInput
        self.surfaceArea = 1.0#surfaceArea

        # Variables
        ## Dynamics
        self.v = -60
        self.w = 0
        self.z = 0.01
        self.ahp = 0.01
        self.m = 0
        self.a_Na = 0.001
        self.b_Na = 0.001
        self.a_Ca = 0.001
        self.b_Ca = 0.001
        self.a_Kv = 0.001

        ## Synaptic Parameters and variables
        self.I_syn = 0
        self.I_ampa = 0
        self.I_gaba = 0
        self.I_nmda = 0
        
        self.s_ampa = [0 for n in range(len(inputs))]
        self.s_gaba = [0 for n in range(len(inputs))]
        self.s_nmda = [0 for n in range(len(inputs))]
        self.s_nmda_x = [0 for n in range(len(inputs))]
        self.s_ampa_dt = [0 for n in range(len(inputs))]
        self.s_gaba_dt = [0 for n in range(len(inputs))]
        self.s_nmda_dt = [0 for n in range(len(inputs))]
        self.s_nmda_x_dt = [0 for n in range(len(inputs))]
        
        self.g_ampa = 7.5e-3
        self.g_gaba = 7.5e-3
        self.g_nmda = 2e-3 # Try setting this to -4 instead, to cut down NMDA contribution
        self.mg = 1e-3
        self.tau_ampa = 2 #ms
        self.tau_gaba = 10 #ms
        self.tau_nmda_rise = 2 #ms
        self.tau_nmda_decay = 100 #ms
        self.alpha = 0.5 # ms
        self.h = 1.0 #Percent Na channels open

        ##M-current related
        Cml = 1      # uF/cm2
        celsius = 36
        tau_m_peak = 1000 / 2.3**((celsius-36)/10)

        ## Instants
        self.m_inf = 0.5*(1+tanh((self.v-self.v_1)/self.v_2))
        self.w_inf = 0.5*(1+tanh((self.v-self.v_3)/self.v_4))
        self.tau_w = 1/(1*cosh((self.v-self.v_3)/(2*self.v_4)))
        self.tau_m = tau_m_peak / ( 3.3 * exp((self.v+35)/20) + exp(-(self.v+35)/20))
        
        # Records
        self.vRecord = []
#        self.dv=[-100 for f in range(len(self.tspan))]
        self.zRecord = []
        self.spikes = 0
        self.spikeRecord = []
        self.currentlySpiking = False
        
    def step(self, time, externalCurrent = 0.0):
        # Get inputs
        vars = {}
        
        # Determine Inputs
        self.currentExcitatoryInputs = []
        self.currentInhibitoryInputs = []
        self.I_syn = 0
        self.I_ampa = 0
        self.I_gaba = 0
        self.I_nmda = 0
        
        if self.debug:
            print "Gaba by synapse:"
        for a in range(len(self.inputs)):
            axon = self.inputs[a]
            tempInput = axon.getInput(self)
            if tempInput > 0 and self.weights[a] > 0: # Glutamatergic input
                self.s_ampa_dt[a] = -(self.s_ampa[a] / self.tau_ampa) + 1 
                self.s_gaba_dt[a] = -(self.s_gaba[a] / self.tau_gaba) 
                self.s_nmda_x_dt[a] = -(self.s_nmda_x[a] / self.tau_nmda_rise) + 1
                self.s_nmda_dt[a] = -(self.s_nmda[a] / self.tau_nmda_decay) + self.alpha*self.s_nmda_x[a]*(1-self.s_nmda[a])
            elif tempInput > 0 and self.weights[a] < 0: #Gabaergic input 
                self.s_ampa_dt[a] = -self.s_ampa[a] / self.tau_ampa
                self.s_gaba_dt[a] = -(self.s_gaba[a] / self.tau_gaba) + 1 
                self.s_nmda_x_dt[a] = -(self.s_nmda_x[a] / self.tau_nmda_rise)
                self.s_nmda_dt[a] = -(self.s_nmda[a] / self.tau_nmda_decay) + self.alpha*self.s_nmda_x[a]*(1-self.s_nmda[a])
            else: # No input
                self.s_ampa_dt[a] = -self.s_ampa[a] / self.tau_ampa
                self.s_gaba_dt[a] = -self.s_gaba[a] / self.tau_gaba
                self.s_nmda_x_dt[a] = -(self.s_nmda_x[a] / self.tau_nmda_rise)
                self.s_nmda_dt[a] = -(self.s_nmda[a] / self.tau_nmda_decay) + self.alpha*self.s_nmda_x[a]*(1-self.s_nmda[a])
            
            # Now update the actual value of the dynamic terms
            self.s_ampa[a] = self.s_ampa[a] + self.s_ampa_dt[a]*self.tau
            self.s_gaba[a] = self.s_gaba[a] + self.s_gaba_dt[a]*self.tau
            if self.debug:
                print "Synapse", a, ":", self.s_gaba[a]
            self.s_nmda[a] = self.s_nmda[a] + self.s_nmda_dt[a]*self.tau
            self.s_nmda_x[a] = self.s_nmda_x[a] + self.s_nmda_x_dt[a]*self.tau
                
        self.I_ampa = sum([self.weights[a] * self.s_ampa[a] for a in range(len(self.inputs))])
        self.I_gaba = sum([-self.weights[a] * self.s_gaba[a] for a in range(len(self.inputs))])
        self.I_nmda = sum([self.weights[a] * self.s_nmda[a] for a in range(len(self.inputs))])
        self.ampaOpen = self.I_ampa
        self.gabaOpen = self.I_gaba
        self.nmdaOpen = self.I_nmda
        self.I_ampa = self.I_ampa * (self.g_ampa * (self.v - self.vNa))
        self.I_gaba = self.I_gaba * (self.g_gaba * (self.v - self.vLeak))
        self.I_nmda = self.I_nmda * ((self.g_nmda * (self.v - self.vNa)) / 1 + self.mg * exp((-0.062*self.v) / 3.57))
        self.I_syn = self.I_ampa + self.I_gaba + self.I_nmda
        Itotal = externalCurrent + self.externalInput - self.I_syn

        celsius = 36

        # I_Na Activation
        alpha_Na_a = -0.182*(self.v+25) / (1 - exp(-1*(self.v+25)/9))
        beta_Na_a = -0.124*(self.v+25) / (1 - exp((self.v+25)/9))
        a_inf_Na = alpha_Na_a / (alpha_Na_a + beta_Na_a)
        tau_a_Na = 1 / alpha_Na_a + beta_Na_a
        da_Na = (a_inf_Na - self.a_Na) / tau_a_Na
        self.a_Na = self.a_Na + da_Na

        # I_Na Inactivation
        alpha_Na_b = 0.024*(self.v+40) / (1 - exp(-1*(self.v+25)/9))
        beta_Na_b = -0.0091*(self.v+65) / (1 - exp((self.v+65)/5))
        b_inf_Na = alpha_Na_b / (alpha_Na_b + beta_Na_b)
        tau_b_Na = 1 / alpha_Na_b + beta_Na_b
        db_Na = (b_inf_Na - self.b_Na) / tau_b_Na
        self.b_Na = self.b_Na + db_Na

        # I_Ca Activation
        alpha_Ca_a = 0.055*(self.v + 27) / (1 - exp(-1*(27+self.v)/3.8))
        beta_Ca_a = 0.94*exp(-1*(self.v+75)/17)
        a_inf_Ca = alpha_Ca_a / (alpha_Ca_a + beta_Ca_a)
        tau_a_Ca = 1 / alpha_Ca_a + beta_Ca_a
        da_Ca = (a_inf_Ca - self.a_Ca) / tau_a_Ca
        self.a_Ca = self.a_Ca + da_Ca

        # I_Ca Inactivation
        alpha_Ca_b = 4.57e-4*exp((self.v+13)/50)
        beta_Ca_b = -0.0065/(1+exp(-1*(self.v+15)/28))
        b_inf_Ca = alpha_Ca_b / (alpha_Ca_b + beta_Ca_b)
        tau_b_Ca = 1 / alpha_Ca_b + beta_Ca_b
        db_Ca = (b_inf_Ca - self.b_Ca) / tau_b_Ca
        self.b_Ca = self.b_Ca + db_Ca

        # I_Kv Activation
        alpha_Kv_a = 0.02*(self.v-25) / (1 - exp(-1*(self.v-25)/9))
        beta_Kv_a = 0.002*(self.v-25) / (1 - exp(-1*(self.v-25)/9))
        a_inf_Kv = alpha_Kv_a / (alpha_Kv_a + beta_Kv_a)
        tau_a_Kv = 1 / alpha_Kv_a + beta_Kv_a
        da_Kv = (a_inf_Kv - self.a_Kv) / tau_a_Kv
        self.a_Kv = self.a_Kv + da_Kv

        # M-Current
        self.tau_z_peak = 1000 / 2.3**((celsius-36)/10)
        self.tau_z = self.tau_z_peak / ( 3.3 * exp((self.v-self.beta)/20) + exp(-(self.v-self.beta)/20))
        dz = ((1 / (1 + exp((self.beta - self.v) / self.gamma)) - self.z) / (self.tau_z))
        self.z = self.z + dz
        
        # if self.debug:
        #     print "Itotal:", Itotal
        #     print "Na Passive:", -self.g_Na*self.m_inf*(self.v-self.v_Na)
        #     print "W:", self.w
        #     print "K Voltage Diff:", (self.v-self.v_K)
        #     print "K Passive:", -self.g_K*self.w*(self.v-self.v_K)
        #     print "Shunting Passive:", -self.g_shunt*(self.v-self.v_shunt)
        #     print "M-Current:", -self.gdapt*self.z*(self.v-self.v_K)

        # Currents
        self.I_Na = self.surfaceArea*self.gNa*self.a_Na**3*self.b_Na*(self.v - self.vNa)
        self.I_Ca = self.surfaceArea*self.gCa*self.a_Ca**2*self.b_Ca*(self.v - self.vCa)
        self.I_Kv = self.surfaceArea*self.gKv*self.a_Kv**1*(self.v - self.vK)
        self.I_Leak = self.surfaceArea*self.gLeak*(self.v - self.vLeak)
        self.I_M = self.surfaceArea*self.gM*self.z*(self.v-self.vK)

        self.IExternal = (self.I_syn + self.I_Na + self.I_Ca + self.I_Kv + self.I_Leak + self.I_M) / (self.C*(1/self.tau))
        self.dv_temp = self.IExternal + externalCurrent / (self.C*(1/self.tau))
        self.v = self.v + self.dv_temp


        self.vRecord.append(self.v)
        self.zRecord.append(self.z)

        if self.currentlySpiking is False and self.v > 0:
            print "SPIKE!"
            self.spikes += 1
            self.spikeRecord.append([time,1])
            self.currentlySpiking = True
            for output in self.outputs:
                output.enqueue()
        if self.currentlySpiking == True and self.v < -10:
            self.currentlySpiking = False        

    def getState(self):
        vars = {}
        vars["v"] = self.v
        vars["i"] = self.I
        return vars
    
    def addInput(self, incomingAxon, weight):
        incomingAxon.addTarget(self)
        self.inputs.append(incomingAxon)
        self.weights.append(weight)
        self.s_ampa.append(0)
        self.s_gaba.append(0)
        self.s_nmda.append(0)
        self.s_nmda_x.append(0)
        self.s_ampa_dt.append(0)
        self.s_gaba_dt.append(0)
        self.s_nmda_dt.append(0)
        self.s_nmda_x_dt.append(0)
        
    def addOutput(self, output):
        self.outputs.append(output)
        
    def getInputs(self):
        return self.inputs