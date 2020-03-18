import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.QuantumNumbers import QuantumNumbers
from pyarts.classes.SpeciesTag import SpeciesTag


class QuantumIdentifier:
    """ ARTS QuantumIdentifier data

    Properties:
        type:
            Type of identifier (get: Index; set: Index or str)

        spec:
            Species of the identifier (Index)

        isot:
            Isotopologue of the identifier (Index)

        lowerqn:
            Lower state quantum numbers of the identifier (QuantumNumbers)

        upperqn:
            Upper state quantum numbers of the identifier (QuantumNumbers)

        levelqn:
            State quantum numbers of the identifier (QuantumNumbers)
    """
    def __init__(self, type="NONE", spec=0, isot=0, qns1=QuantumNumbers(), qns2=QuantumNumbers()):
        if isinstance(type, c.c_void_p):
            self.__delete__ = False
            self.__data__ = type
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumIdentifier())
            self.type = type
            self.spec = spec
            self.isot = isot
            self.upperqn = qns1
            self.lowerqn = qns2

    @property
    def type(self):
        """ Type of identifier (get: Index; set: Index or str) """
        return lib.getQuantumIdentifierType(self.__data__)

    @type.setter
    def type(self, type):
        if isinstance(type, str):
            type = type.encode("ascii")
            if lib.setQuantumIdentifierTypeFromString(self.__data__, type):
                raise ValueError("Invalid type")
        else:
            type = int(type)
            if lib.setQuantumIdentifierTypeFromIndex(self.__data__, type):
                raise ValueError("Invalid type")

    @property
    def spec(self):
        """ Species of the identifier (Index) """
        return lib.getQuantumIdentifierSpecies(self.__data__)

    @spec.setter
    def spec(self, val):
        if not SpeciesTag.validSpecies(val):
            raise ValueError("Invalid species")
        lib.setQuantumIdentifierSpecies(self.__data__, c.c_long(val))

    @property
    def isot(self):
        """ Isotopologue of the identifier (Index) """
        return lib.getQuantumIdentifierIsotopologue(self.__data__)

    @isot.setter
    def isot(self, val):
        if not SpeciesTag.validIsotopologue(self.spec, val):
            raise ValueError("Invalid isotopologue")
        lib.setQuantumIdentifierIsotopologue(self.__data__, int(val))

    @property
    def lowerqn(self):
        """ Lower state quantum numbers of the identifier (QuantumNumbers) """
        return QuantumNumbers(c.c_void_p(lib.getQuantumIdentifierLowerQuantumNumbers(self.__data__)))

    @lowerqn.setter
    def lowerqn(self, val):
        self.lowerqn.set(val)

    @property
    def upperqn(self):
        """ Upper state quantum numbers of the identifier (QuantumNumbers) """
        return QuantumNumbers(c.c_void_p(lib.getQuantumIdentifierUpperQuantumNumbers(self.__data__)))

    @upperqn.setter
    def upperqn(self, val):
        self.upperqn.set(val)

    @property
    def levelqn(self):
        return QuantumNumbers(c.c_void_p(lib.getQuantumIdentifierLevelQuantumNumbers(self.__data__)))

    @levelqn.setter
    def levelqn(self, val):
        """ State quantum numbers of the identifier (QuantumNumbers) """
        self.levelqn.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumIdentifier(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumIdentifier(self.__data__)

    def __repr__(self):
        return "ARTS QuantumIdentifier"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, QuantumIdentifier):
            self.type = other.type
            self.spec = other.spec
            self.isot = other.isot
            self.lowerqn = other.lowerqn
            self.upperqn = other.upperqn
            self.levelqn = other.levelqn  # repeat but for now keep...
        else:
            raise TypeError("Expects QuantumIdentifier")


lib.createQuantumIdentifier.restype = c.c_void_p
lib.createQuantumIdentifier.argtypes = []

lib.deleteQuantumIdentifier.restype = None
lib.deleteQuantumIdentifier.argtypes = [c.c_void_p]

lib.printQuantumIdentifier.restype = None
lib.printQuantumIdentifier.argtypes = [c.c_void_p]

lib.getQuantumIdentifierType.restype = c.c_long
lib.getQuantumIdentifierType.argtypes = [c.c_void_p]

lib.setQuantumIdentifierTypeFromString.restype = c.c_long
lib.setQuantumIdentifierTypeFromString.argtypes = [c.c_void_p, c.c_char_p]

lib.setQuantumIdentifierTypeFromIndex.restype = c.c_long
lib.setQuantumIdentifierTypeFromIndex.argtypes = [c.c_void_p, c.c_long]

lib.getQuantumIdentifierSpecies.restype = c.c_long
lib.getQuantumIdentifierSpecies.argtypes = [c.c_void_p]

lib.setQuantumIdentifierSpecies.restype = None
lib.setQuantumIdentifierSpecies.argtypes = [c.c_void_p, c.c_long]

lib.getQuantumIdentifierIsotopologue.restype = c.c_long
lib.getQuantumIdentifierIsotopologue.argtypes = [c.c_void_p]

lib.setQuantumIdentifierIsotopologue.restype = None
lib.setQuantumIdentifierIsotopologue.argtypes = [c.c_void_p, c.c_long]

lib.getQuantumIdentifierLevelQuantumNumbers.restype = c.c_void_p
lib.getQuantumIdentifierLevelQuantumNumbers.argtypes = [c.c_void_p]

lib.getQuantumIdentifierLowerQuantumNumbers.restype = c.c_void_p
lib.getQuantumIdentifierLowerQuantumNumbers.argtypes = [c.c_void_p]

lib.getQuantumIdentifierUpperQuantumNumbers.restype = c.c_void_p
lib.getQuantumIdentifierUpperQuantumNumbers.argtypes = [c.c_void_p]
