from __future__ import division
from InjectedCurrentReceptor import *
from GlutGSDReceptor import *
from GABAGSDReceptor import *
import logging as log
from random import random

'''
Created on September 12, 2018

A class which handles a transmission queue, weight, and STDP if needed.

@author: Stefan Berteau
'''


class Axon():
    def __init__(self, timeStep, weight, source, target, postSynapticReceptors=None, length=0, width=0.1, myelin=0.1, debug=False):
        self.timeStep = timeStep
        self.debug = debug
        self.weight = weight
        self.source = source
        self.target = target

        source.addOutput(self)
        target.addInput(self)

        # If not otherwise specified, set up a standard injection of "weight" current in picoamperes at the time of each spike.
        if postSynapticReceptors is None:
            if self.weight <= 0:
                self.postSynapticReceptors = [GABAGSDReceptor(self.target, self.weight, self.timeStep)]
            else:
                self.postSynapticReceptors = [GlutGSDReceptor(self.target, self.weight, self.timeStep)]
            # self.postSynapticReceptors = [InjectedCurrentReceptor(self.target, self.weight)]
        elif isinstance(postSynapticReceptors, list):
            self.postSynapticReceptors = postSynapticReceptors
        else:
            raise ValueError("Error: postSynapticReceptors should be a list of PostSynapticReceptor objects")

        if length == 0:
            self.fastMode = True
        else:
            self.spikes = [] # Spikes currently traveling down the axon
            self.length = length / 1000  # Convert length to (look this dimension up)
            self.width = width
            self.myelin = myelin
            self.speed = (self.width ** (2 - self.myelin)) * 1000  # determine conduction velocity

        self.failureRate = 0.1

        self.proximalDiffuseReceptors = []
        self.distalDiffuseReceptors = []

        self.proximalDiffuseTransmitterLevels = {}
        self.distalDiffuseTransmitterLevels = {}

        if self.debug:
            print(self.source.type, "connected to ", self.target.type)

    def addProximalDiffuseAxonReceptor(self, receptor):
        self.proximalDiffuseReceptors.append(receptor)

    def addDistalDiffuseAxonReceptor(self, receptor):
        self.distalDiffuseReceptors.append(receptor)

    def updateProximalDiffuseTransmitters(self, proximalDiffuseTransmitterLevels):
        self.proximalDiffuseTransmitterLevels = proximalDiffuseTransmitterLevels
        for proximalAxonReceptor in self.proximalDiffuseReceptors:  # For each proximal receptor
            try:
                print(self.proximalDiffuseTransmitterLevels[proximalAxonReceptor.getTypeString()])
                proximalAxonReceptor.setLevel(self.proximalDiffuseTransmitterLevels[proximalAxonReceptor.getTypeString()]) # Set the level to be whatever matches it in the passed dictionary
            except KeyError:
                log.warning("WARNING: There are proximal receptors of type " + proximalAxonReceptor.getTypeString() + " in axons originating in population " + self.source.parentPopulation.name + ", which does not currently have any transmitter levels set for that name.  Defaulting to a level of 0.0")
                proximalAxonReceptor.setLevel(0.0)
        return

    def updateDistalDiffuseTransmitters(self, distalDiffuseTransmitters):
        self.distalDiffuseTransmitters = distalDiffuseTransmitters
        for distalAxonReceptor in self.distalDiffuseReceptors:  # For each distal receptor
            try:
                distalAxonReceptor.setLevel(self.distalDiffuseTransmitters[distalAxonReceptor.getTypeString()]) # Set the level to be whatever matches it in the passed dictionary
            except KeyError:
                log.warning("WARNING: There are distal receptors of type " + distalAxonReceptor.getTypeString() + " in axons going to population " + self.target.parentPopulation.name + " from " + self.source.parentPopulation.name + ", which does not currently have any transmitter levels set for that name.  Defaulting to a level of 0.0")
                distalAxonReceptor.setLevel(0.0)
        return

    def enqueue(self):
        if random() > self.failureRate:
            if self.debug:
                print(self.source.type, "spiked!  Transmitting ", self.weight, "to ", self.target.type)
            if self.fastMode:   # If we are not doing timed myelinated transmission, just directly inject current into the target
                self.boutonSpike()
            else:
                self.spikes.append(0)

    def boutonSpike(self):
        for receptor in self.postSynapticReceptors:
            receptor.boutonSpike()

    def injectReceptorCurrent(self):
        for receptor in self.postSynapticReceptors:
            receptor.step()
            receptor.injectCurrent()

    def step(self):
        # Inject Current
        self.injectReceptorCurrent()

        # If we care about action potential transit times, update action potential positions
        if self.fastMode:
            return
        else:
            self.time += self.timeStep
            self.speed = (self.width ** (2 - self.myelin)) * 1000
            self.spikes = [self.spikes[i] + (self.timeStep / self.speed) for i in range(len(self.spikes))]
            while len(self.spikes) > 0 and max(self.spikes) > self.length:
                self.spikes.remove(max(self.spikes))
                self.boutonSpike()

                # TODO: Add Voltage dependent stdp receptors?
                # TODO: Add current-based NMDA, GABA, and AMPA receptors