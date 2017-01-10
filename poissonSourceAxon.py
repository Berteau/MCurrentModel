from __future__ import division
from numpy import *
from scipy.stats import poisson

'''
Created on Jan 8, 2011

An axon class that provides external input in the form of a poisson distribution of spikes.  It still maxes out at one spike per time point, however, so we need a number of them to get a serious distribution.

@author: stefan
'''

class PoissonAxon:
    def __init__(self, timeStep, lambdaConst = 0.5, onsetTime = 0, offsetTime = 10000000):
        self.lambdaConst = lambdaConst
        self.timeStep = timeStep
        self.time = 0
        self.onsetTime = onsetTime
        self.offsetTime = offsetTime
        self.distribution = poisson(self.lambdaConst)

    def getInput(self, target):
        if self.time < self.onsetTime or self.time > self.offsetTime:
            return 0
        else:
            return self.distribution.rvs((1,))[0]

    def addTarget(self, target):
        return True

    def step(self): 
        self.time += self.timeStep