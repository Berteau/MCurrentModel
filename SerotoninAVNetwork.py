import logging as log
from Population import *
from PoissonPopulation import *
from Network import *
from AxonalSerotoninReceptorFactory import *
from SomaticSerotoninReceptorFactory import *
from random import random, gauss
import functools

# TODO: change imports to a "import <file>" format, then call them with <file>.<class>, as best practices

class SerotoninAVNetwork(Network):
    def __init__(self, tau, parentSimulation, params, name):
        self.tau = tau
        self.parentSimulation = parentSimulation
        self.name = name
        self.populations = {}
        self.params = params

        # try:
        self.popCount = params["popCount"]
        self.serotoninLevelV = params["serotoninLevelV"]
        self.serotoninLevelA = params["serotoninLevelA"]
        [self.inputWeightA1, self.inputWeightA2] = [params["inputWeightA1"], params["InputWeightA2"]]
        [self.rateA1, self.rateA2] = [params["rateA1"], params["rateA2"]]
        [self.inputWeightV1, self.inputWeightV2, self.inputWeightV3] = [params["inputWeightV1"], params["InputWeightV2"], params["InputWeightV3"]]
        [self.rateV1, self.rateV2, self.rateV3] = [params["rateV1"], params["rateV2"], params["rateV3"]]
        self.crossModalWeightVA = params["inputWeightVA"]
        self.crossModalWeightAV = params["inputWeightAV"]
        self.somaticSerotonin2AReceptorWeight = params["Somatic5HT2AWeight"]
        self.somaticSerotonin1AReceptorWeight = params["Somatic5HT1AWeight"]
        self.axonalSerotonin2AReceptorWeight = params["Axonal5HT2AWeight"]
        self.axonalSerotonin1AReceptorWeight = params["Axonal5HT1AWeight"]
        self.pyramidalSelfExcitationWeight = params["pyramidalSelfExcitationWeight"]
        # except KeyError:
        #     er.
        #     log.CRITICAL("Network " + self.name + " cannot retrieve correct parameters from dictionary.  Aborting network construction and returning Nonetype")
        #     raise ValueError("Error: Incorrect Network Parameters!")
        #     return

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

        ltsParams = {}
        ltsParams["C"] = 100
        ltsParams["k"] = 1
        ltsParams["v_r"] = -56
        ltsParams["v_t"] = -42
        ltsParams["v_peak"] = 40
        ltsParams["a"] = 0.03
        ltsParams["b"] = 8
        ltsParams["c"] = -50
        ltsParams["d"] = 200
        ltsParams["type"] = "LTS"

        inputV1Params = {}
        inputV1Params["spikeRate"] = self.rateV1
        inputV1Params["type"] = "Poisson"

        inputV2Params = {}
        inputV2Params["spikeRate"] = self.rateV2
        inputV2Params["type"] = "Poisson"

        inputV3Params = {}
        inputV3Params["spikeRate"] = self.rateV3
        inputV3Params["type"] = "Poisson"

        inputA1Params = {}
        inputA1Params["spikeRate"] = self.rateA1
        inputA1Params["type"] = "Poisson"

        inputA2Params = {}
        inputA2Params["spikeRate"] = self.rateA2
        inputA2Params["type"] = "Poisson"

        # Define defuse receptor factories
        factorySomatic5HT2A = SomaticSerotoninDiffuseReceptorFactory("5HT2A", lambda x: self.somaticSerotonin2AReceptorWeight, [])
        factorySomatic5HT1A = SomaticSerotoninDiffuseReceptorFactory("5HT1A", lambda x: self.somaticSerotonin1AReceptorWeight, [])
        factoryAxonal5HT2A = AxonalSerotoninDiffuseReceptorFactory("5HT2A", lambda x: self.axonalSerotonin2AReceptorWeight, [])
        factoryAxonal5HT1A = AxonalSerotoninDiffuseReceptorFactory("5HT1A", lambda x: self.axonalSerotonin1AReceptorWeight, [])

        # Define diffuse Transmitters for Visual Area
        transmittersV = {}
        transmittersV["5HT2A"] = self.serotoninLevelV
        transmittersV["5HT1A"] = self.serotoninLevelV

        # Define diffuse transmitters for Auditory Area
        transmittersA = {}
        transmittersA["5HT2A"] = self.serotoninLevelA
        transmittersA["5HT1A"] = self.serotoninLevelA

        # Build the populations
        self.populations["InputV1"] = PoissonPopulation(tau, inputV1Params, self.popCount, [], {}, self, name)
        self.populations["InputV2"] = PoissonPopulation(tau, inputV2Params, self.popCount, [], {}, self, name)
        self.populations["InputV3"] = PoissonPopulation(tau, inputV3Params, self.popCount, [], {}, self, name)
        self.populations["InputA1"] = PoissonPopulation(tau, inputA1Params, self.popCount, [], {}, self, name)
        self.populations["InputA2"] = PoissonPopulation(tau, inputA2Params, self.popCount, [], {}, self, name)

        self.populations["pyramidalsV1"] = Population(self.tau, pyramidalParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V1 Pyramidal Cells")
        self.populations["fastSpikingsV1"] = Population(self.tau, fsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V1 Fast Spiking Cells")
        self.populations["lowThresholdsV1"] = Population(self.tau, ltsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V1 Low Threshold Cells")

        self.populations["pyramidalsV2"] = Population(self.tau, pyramidalParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V2 Pyramidal Cells")
        self.populations["fastSpikingsV2"] = Population(self.tau, fsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V2 Fast Spiking Cells")
        self.populations["lowThresholdsV2"] = Population(self.tau, ltsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V2 Low Threshold Cells")

        self.populations["pyramidalsV3"] = Population(self.tau, pyramidalParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V3 Pyramidal Cells")
        self.populations["fastSpikingsV3"] = Population(self.tau, fsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V3 Fast Spiking Cells")
        self.populations["lowThresholdsV3"] = Population(self.tau, ltsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "V3 Low Threshold Cells")

        self.populations["pyramidalsA1"] = Population(self.tau, pyramidalParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "A1 Pyramidal Cells")
        self.populations["fastSpikingsA1"] = Population(self.tau, fsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "A1 Fast Spiking Cells")
        self.populations["lowThresholdsA1"] = Population(self.tau, ltsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "A1 Low Threshold Cells")

        self.populations["pyramidalsA2"] = Population(self.tau, pyramidalParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "A2 Pyramidal Cells")
        self.populations["fastSpikingsA2"] = Population(self.tau, fsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "A2 Fast Spiking Cells")
        self.populations["lowThresholdsA2"] = Population(self.tau, ltsParams, self.popCount, [factorySomatic5HT2A], transmittersV, self, "A2 Low Threshold Cells")

        # Build connections

        # First, visual input
        self.populations["InputV1"].addOutboundConnections(self.populations["pyramidalsV1"], self.gaussWeightWithChanceFactory(self.inputWeightV1, 1.0), [], [factoryAxonal5HT2A])
        self.populations["InputV2"].addOutboundConnections(self.populations["pyramidalsV2"], self.gaussWeightWithChanceFactory(self.inputWeightV2, 1.0), [], [factoryAxonal5HT2A])
        self.populations["InputV3"].addOutboundConnections(self.populations["pyramidalsV3"], self.gaussWeightWithChanceFactory(self.inputWeightV3, 1.0), [], [factoryAxonal5HT2A])

        # Then auditory input
        self.populations["InputA1"].addOutboundConnections(self.populations["pyramidalsA1"], self.gaussWeightWithChanceFactory(self.inputWeightA1, 1.0), [], [factoryAxonal5HT2A])
        self.populations["InputA2"].addOutboundConnections(self.populations["pyramidalsA2"], self.gaussWeightWithChanceFactory(self.inputWeightA2, 1.0), [], [factoryAxonal5HT2A])

        # Then cross-modal input
        self.populations["InputV1"].addOutboundConnections(self.populations["pyramidalsA1"], self.gaussWeightWithChanceFactory(self.crossModalWeightVA, 0.3), [], [factoryAxonal5HT2A])
        self.populations["InputV3"].addOutboundConnections(self.populations["pyramidalsA2"], self.gaussWeightWithChanceFactory(self.crossModalWeightVA, 0.3), [], [factoryAxonal5HT2A])

        self.populations["InputA1"].addOutboundConnections(self.populations["pyramidalsV1"], self.gaussWeightWithChanceFactory(self.crossModalWeightAV, 0.3), [], [factoryAxonal5HT2A])
        self.populations["InputA2"].addOutboundConnections(self.populations["pyramidalsV3"], self.gaussWeightWithChanceFactory(self.crossModalWeightAV, 0.3), [], [factoryAxonal5HT2A])

        # Now we do connections that run within each group (column, if you will permit me to use the term loosly) of populations

        # Pyramidal self-excitation
        self.populations["pyramidalsV1"].addOutboundConnections(self.populations["pyramidalsV1"], self.selfTargetingGaussWeightFactory(self.params["pyramidalSelfExcitationWeight"]), [], [])
        self.populations["pyramidalsV2"].addOutboundConnections(self.populations["pyramidalsV2"], self.selfTargetingGaussWeightFactory(self.params["pyramidalSelfExcitationWeight"]), [], [])
        self.populations["pyramidalsV3"].addOutboundConnections(self.populations["pyramidalsV3"], self.selfTargetingGaussWeightFactory(self.params["pyramidalSelfExcitationWeight"]), [], [])
        self.populations["pyramidalsA1"].addOutboundConnections(self.populations["pyramidalsA1"], self.selfTargetingGaussWeightFactory(self.params["pyramidalSelfExcitationWeight"]), [], [])
        self.populations["pyramidalsA2"].addOutboundConnections(self.populations["pyramidalsA2"], self.selfTargetingGaussWeightFactory(self.params["pyramidalSelfExcitationWeight"]), [], [])

        # Now Pyramidal to Fast-Spiking
        self.populations["pyramidalsV1"].addOutboundConnections(self.populations["fastSpikingsV1"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToFSWeight"], 1.0), [], [])
        self.populations["pyramidalsV2"].addOutboundConnections(self.populations["fastSpikingsV2"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToFSWeight"], 1.0), [], [])
        self.populations["pyramidalsV3"].addOutboundConnections(self.populations["fastSpikingsV3"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToFSWeight"], 1.0), [], [])
        self.populations["pyramidalsA1"].addOutboundConnections(self.populations["fastSpikingsA1"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToFSWeight"], 1.0), [], [])
        self.populations["pyramidalsA2"].addOutboundConnections(self.populations["fastSpikingsA2"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToFSWeight"], 1.0), [], [])

        # Fast Spiking to Pyramidal
        self.populations["fastSpikingsV1"].addOutboundConnections(self.populations["pyramidalsV1"], self.gaussWeightWithChanceFactory(self.params["FSToPyramidalsWeight"], 1.0), [], [])
        self.populations["fastSpikingsV2"].addOutboundConnections(self.populations["pyramidalsV2"], self.gaussWeightWithChanceFactory(self.params["FSToPyramidalsWeight"], 1.0), [], [])
        self.populations["fastSpikingsV3"].addOutboundConnections(self.populations["pyramidalsV3"], self.gaussWeightWithChanceFactory(self.params["FSToPyramidalsWeight"], 1.0), [], [])
        self.populations["fastSpikingsA1"].addOutboundConnections(self.populations["pyramidalsA1"], self.gaussWeightWithChanceFactory(self.params["FSToPyramidalsWeight"], 1.0), [], [])
        self.populations["fastSpikingsA2"].addOutboundConnections(self.populations["pyramidalsA2"], self.gaussWeightWithChanceFactory(self.params["FSToPyramidalsWeight"], 1.0), [], [])

        # Pyramidal to Low-Threshold
        self.populations["pyramidalsV1"].addOutboundConnections(self.populations["lowThresholdsV1"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToLTSWeight"], 1.0), [], [])
        self.populations["pyramidalsV2"].addOutboundConnections(self.populations["lowThresholdsV2"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToLTSWeight"], 1.0), [], [])
        self.populations["pyramidalsV3"].addOutboundConnections(self.populations["lowThresholdsV3"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToLTSWeight"], 1.0), [], [])
        self.populations["pyramidalsA1"].addOutboundConnections(self.populations["lowThresholdsA1"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToLTSWeight"], 1.0), [], [])
        self.populations["pyramidalsA2"].addOutboundConnections(self.populations["lowThresholdsA2"], self.gaussWeightWithChanceFactory(self.params["PyramidalsToLTSWeight"], 1.0), [], [])

        # TODO: Add cross-modal connections

        # TODO: Add competition between visual columns

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
