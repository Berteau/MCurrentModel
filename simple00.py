# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:08:24 2019

@author: yohan
"""

from Population import *
from PoissonPopulation import *
from Network import *
from random import gauss

class simple00(Network):

    def __init__(self, tau, parentSimulation, params, name):
        self.tau = tau
        self.parentSimulation = parentSimulation
        self.name = name
        self.populations = {}
        self.params = params
        self.popCount = self.params["popCount"]
        
        pyramidalParams = {}
        pyramidalParams["C"] = 100
        pyramidalParams["k"] = 3
        pyramidalParams["v_r"] = -60
        pyramidalParams["v_t"] = -50
        pyramidalParams["v_peak"] = 50
        pyramidalParams["a"] = 0.02
        pyramidalParams["b"] = 0.2
        pyramidalParams["c"] = -65
        pyramidalParams["d"] = 8
        pyramidalParams["type"] = "Pyramidal"
        
        
        lambdaAParams = {}
        lambdaAParams["spikeRate"] = 7.5
        lambdaAParams["type"] = "Poisson"
        
        self.populations["InputA"] = PoissonPopulation(tau, lambdaAParams, self.popCount, [], {}, self, name)

        self.populations["pyramidalsA"] = Population(self.tau, pyramidalParams, self.popCount, [], {}, self, "Area A Pyramidal Cells")
        
        self.populations["InputA"].addOutboundConnections(self.populations["pyramidalsA"], self.gaussWeightWithChanceFactory(5, 1.0), [], [])
        
        self.populations["pyramidalsA"].addOutboundConnections(self.populations["pyramidalsA"], self.gaussWeightWithChanceFactory(5, 1.0), [], [])
        
    def gaussWeightWithChanceFactory(self, inputWeight, chance):
        def weightFunction(source, target):
            if random() < chance:
                return gauss(inputWeight / self.popCount, ((inputWeight / self.popCount) / 10))
            else:
                return None
        return weightFunction
    
    
    def step(self):
        for popName, population in self.populations.items():
            population.stepCells()
        for popName, population in self.populations.items():
            population.stepOutputs()
            
            
            
