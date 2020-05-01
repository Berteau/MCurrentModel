from DiffuseReceptor import *

class AxonalSerotoninReceptor(DiffuseReceptor):
    def __init__(self, target, typeString, weight, initialLevel):
        self.typeString = typeString
        self.target = target
        self.weight = weight
        self.level = initialLevel
        if self.target is not None:
            self.unmodifiedFailureRate = self.target.failureRate

    def setTarget(self, target):
        if self.target is not None:
            self.target.v_peak = self.unmodifiedFailureRate
        self.target = target
        self.unmodifiedFailureRate = self.target.failureRate

    def setLevel(self, level):
        if isinstance(level, (int,float)):
            self.level = level
            # Update threshold
            self.target.failureRate = self.unmodifiedFailureRate - (self.weight * self.level)
        else:
            raise ValueError("Error: setSerotoninLevel requires an integer or floating point input")

    def doActivity(self):
        return None
