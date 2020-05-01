import logging as log
from Population import *
from PoissonPopulation import *
from Network import *
from AxonalSerotoninReceptorFactory import *
from SomaticSerotoninReceptorFactory import *
from random import random, gauss
import functools

# TODO: change imports to a "import <file>" format, then call them with <file>.<class>, as best practices

class EINetwork(Network):
    def __init__(self, tau, parentSimulation, params, name):
        self.tau = tau
        self.parentSimulation = parentSimulation
        self.name = name
        self.populations = {}
        self.params = params

        # try:
        self.popCount = params["popCount"]
        self.inputWeightA = params["inputWeightA"]
        self.rateA = params["rateA"]
        self.pyramidalSelfExcitationWeight = params["pyramidalSelfExcitationWeight"]
        self.pyramidalToFSWeight = params["pyramidalToFSWeight"]
        self.fsToPyramidalWeight = params["fsToPyramidalWeight"]

        # Set up cell parameters
        pyramidalParams = {}
        pyramidalParams["C"] = 100
        pyramidalParams["k"] = 3
        pyramidalParams["v_r"] = -60
        pyramidalParams["v_t"] = -50
        pyramidalParams["v_peak"] = 50
        pyramidalParams["a"] = 0.01
        pyramidalParams["b"] = 5
        pyramidalParams["c"] = -60
        pyramidalParams["d"] = 400
        pyramidalParams["type"] = "Pyramidal"

        fsParams = {}
        fsParams["C"] = 20
        fsParams["k"] = 1
        fsParams["v_r"] = -55
        fsParams["v_t"] = -40
        fsParams["v_peak"] = 25
        fsParams["a"] = 0.15
        fsParams["b"] = 8
        fsParams["c"] = -55
        fsParams["d"] = 200
        fsParams["type"] = "FS"

        lambdaAParams = {}
        lambdaAParams["spikeRate"] = self.rateA
        lambdaAParams["type"] = "Poisson"

        # Build the populations
        self.populations["Input"] = PoissonPopulation(tau, lambdaAParams, self.popCount, [], {}, self, name)
        self.populations["pyramidals"] = Population(self.tau, pyramidalParams, self.popCount, [], {}, self, "Pyramidal Cells")
        self.populations["fastSpikings"] = Population(self.tau, fsParams, self.popCount, [], {}, self, "Fast Spiking Cells")

        # Build connections

        # First, unimodal input
        self.populations["Input"].addOutboundConnections(self.populations["pyramidals"], self.gaussWeightWithChanceFactory(self.inputWeightA, 1.0), [], [])

        # Now we do connections that run within each column

        # Pyramidal self-excitation
        self.populations["pyramidals"].addOutboundConnections(self.populations["pyramidals"], self.selfTargetingGaussWeightFactory(self.params["pyramidalSelfExcitationWeight"]), [], [])

        # Now Pyramidal to Fast-Spiking
        self.populations["pyramidals"].addOutboundConnections(self.populations["fastSpikings"], self.gaussWeightWithChanceFactory(self.params["pyramidalToFSWeight"], 1.0), [], [])

        # Fast Spiking to Pyramidal
        self.populations["fastSpikings"].addOutboundConnections(self.populations["pyramidals"], self.gaussWeightWithChanceFactory(self.params["fsToPyramidalWeight"], 1.0), [], [])

    # TODO: Separate out connection function factories to another library file somewhere

    def gaussWeightWithChanceFactory(self, inputWeight, chance):
        def weightFunction(source, target):
            if random() < chance:
                return gauss(inputWeight / self.popCount, ((inputWeight / self.popCount) / 10))
            else:
                return None
        return weightFunction

    def selfTargetingGaussWeightFactory(self, inputWeight):
        def weightFunction(source, target):
            if source is target:
                return gauss(inputWeight / self.popCount, ((inputWeight / self.popCount) / 10))
            else:
                return None

        return weightFunction

    def step(self):
        for popName, population in self.populations.items():
            population.stepCells()
        for popName, population in self.populations.items():
            population.stepOutputs()
