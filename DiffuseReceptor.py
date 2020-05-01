from abc import ABC

class DiffuseReceptor(ABC):
    def __init__(self, target, typeString, weight, initialLevel):
        self.typeString = typeString
        self.target = target
        self.weight = weight
        self.level = initialLevel

    def setLevel(self, level):
        if isinstance(level, (int, float)):
            self.level = level
        else:
            raise ValueError("Error: setLevel requires an integer or floating point input")

    def setTarget(self, target):
        self.target = target

    def setTypeString(self, typeString):
        self.typeString = typeString

    def getTypeString(self):
        return self.typeString

    def doActivity(self):
        pass

