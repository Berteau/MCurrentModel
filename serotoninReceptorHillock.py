from DiffuseReceptor import *

class SerotoninReceptorHillock(DiffuseReceptor):
    def __init__(self, target, receptorType, weight, initialLevel):
        self.receptorType = receptorType
        self.target = target
        self.weight = weight
        self.level = initialLevel
        self.unmodifiedThreshold = target.v_peak

    def setTarget(self, target):
        self.target.v_peak = self.unmodifiedThreshold
        self.target = target
        self.unmodifiedThreshold = self.target.v_peak

    def setLevel(self, level):
        if isinstance(level, (int,float)):
            self.level = level
            # Update threshold
            self.target.v_peak = self.unmodifiedThreshold + (self.weight * self.level)
        else:
            raise ValueError("Error: setSerotoninLevel requires an integer or floating point input")

    def doActivity(self):
        return None


