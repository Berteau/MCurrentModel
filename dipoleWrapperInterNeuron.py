from __future__ import division

'''
Created on Sep 16, 2011

@author: stefan
'''
import monopoleInterNeuron
from pylab import *
from scipy import *
from numpy import *

class DipoleWrapperInterNeuron:
    def __init__(self,inputs_a=[], inputs_b=[], outputs=[], externalInput_a = 0.0, externalInput_b = 0.0, position = [0.0, 0.0, 0.0], distance=0.06, debug=False, tau=0.1, tspan=range((10000)), inactivating=True):
        self.debug = debug
        self.tspan = tspan
        self.axialResistance = 25 #25 kOhm/cm, so that's like 
        self.distance = distance #0.5      #0.6332 #from per Henze DA CW, Barrioneuvo G. Dendritic Morphology and its effects on the amplitude and rise-time of synaptic signals in hippocampal CA3 pyramidal cells. Journal of Comparative Neurology. 1996;369:331-344.

        if inputs_a is None:
            inputs_a = []
        if inputs_b is None:
            inputs_b = []
        if outputs is None:
            outputs = []
        apexPosition = position
        apexPosition[1] += distance

#        self.apex = monopoleNeuron_alt.MonopoleNeuronPassive(inputs_a, [], externalInput_a, apexPosition, debug, tau)        
        self.base = monopoleInterNeuron.MonopoleInterNeuron(inputs_b, outputs, externalInput_b, position, debug, tau, tspan=self.tspan, inactivating=inactivating)
        self.apex = self.base
        

        # Time
        self.timestep = tau
        self.tau = tau
        self.position = position

    def step(self, time, externalInputA = 0, externalInputB = 0):
        ## Calculate Bridge Current
#        self.I_Bridge = (self.apex.v - self.base.v) / (self.axialResistance * self.distance)
#        self.I_a = -self.I_Bridge
#        self.I_b = self.I_Bridge
        self.I_Bridge = 0
        I_a = externalInputA
        I_b = externalInputB
        
#        self.apex.step(time, I_a)
        self.base.step(time, I_b)

    def getState(self):
        vars = {}
        vars["v"] = self.base.v
        vars["i_a"] = self.I_a
        vars["i_b"] = self.I_b
        vars["I_Bridge"] = self.I_Bridge
        return vars 
    
    def addInput_a(self, tempInput, weight):
        self.apex.addInput(tempInput, weight)
        
    def addInput_b(self, tempInput, weight):
        self.base.addInput(tempInput, weight)
        
    def addOutput(self, output):
        self.base.addOutput(output)
        
    def getInputs(self):
        return self.inputs

    def isBursting(self, inputOnset = 500):
        onsetSpikeRate = mean(self.base.spikeRecord[inputOnset:inputOnset+100])
        sustainedSpikeRate = mean(self.basel.spikeRecord[inputOnset+101:-1])
        if onsetSpikeRate / sustainedSpikeRate > 3.0:
            return 1
        else:
            return 0