from __future__ import division
from numpy import *

'''
Created on Jan 8, 2011

An axon class that provides external input in the form of a poisson distribution of spikes.  It still maxes out at one spike per time point, however, so we need a number of them to get a serious distribution.

@author: stefan
'''

class NeuronNeuronAxon:
    def __init__(self, timeStep, lengths, source, targets, width=0.1, myelin=0.1):
        self.lengths = [length/1000 for length in lengths] # One length for each target too.
        self.spikes = [[] for t in range(len(targets))] # One list of spikes for each target
        
        # Each target passes itself when getting input
        # Check given the length appropriate to that target 
        
        self.source = source
        self.source.addOutput(self)
        self.targets = targets
#        for target in self.targets:
#            target.addInput(self)
        self.time = 0
        self.timeStep = timeStep
        self.width = width
        self.myelin = myelin
        self.speed = (self.width**(2-self.myelin)) * 1000
        
    def addTarget(self, target):
        self.targets.append(target)
        self.spikes.append([])
        self.lengths.append(1/1000)
        
    def getInput(self, target):
        targetIndex = self.targets.index(target)
        if len(self.spikes[targetIndex]) == 0:
            return 0
        elif max(self.spikes[targetIndex]) > self.lengths[targetIndex]:
            self.spikes[targetIndex].remove(max(self.spikes[targetIndex]))
            return 1
        else:
            return 0
    
    def enqueue(self):
        for s in range(len(self.spikes)):
            self.spikes[s].append(0) 

    def step(self):
        self.time += self.timeStep
        self.speed = (self.width**(2-self.myelin)) * 1000
        self.spikes = [[self.spikes[a][i]+(self.timeStep/self.speed) for i in range(len(self.spikes[a]))] for a in range(len(self.spikes))]    
    