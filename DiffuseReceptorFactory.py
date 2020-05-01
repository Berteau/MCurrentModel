from abc import ABC

class DiffuseReceptorFactory(ABC):
    def __init__(self, typeString, weightFunction, weightFunctionParamters):
        self.typeString = typeString
        self.weightFunction = weightFunction
        self.weightFunctionParamters

    def getTypeString(self):
        return self.typeString

    def constructReceptor(self, initialLevel=0):
        return None
