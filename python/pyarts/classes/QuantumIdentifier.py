import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.QuantumNumbers import QuantumNumbers
from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


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
        return lib.getTypeQuantumIdentifier(self.__data__)

    @type.setter
    def type(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexTypeQuantumIdentifier(self.__data__, type.encode("ascii")))
        else:
            type = int(type)
            if lib.setTypeQuantumIdentifier(self.__data__, type):
                raise ValueError("Invalid type")

    @property
    def spec(self):
        """ Species of the identifier (Index) """
        return lib.getSpeciesQuantumIdentifier(self.__data__)

    @spec.setter
    def spec(self, val):
        if not SpeciesTag.validSpecies(int(val)):
            raise ValueError("Invalid species")
        lib.setSpeciesQuantumIdentifier(self.__data__, int(val))

    @property
    def isot(self):
        """ Isotopologue of the identifier (Index) """
        return lib.getIsotopologueQuantumIdentifier(self.__data__)

    @isot.setter
    def isot(self, val):
        if not SpeciesTag.validIsotopologue(self.spec, int(val)):
            raise ValueError("Invalid isotopologue")
        lib.setIsotopologueQuantumIdentifier(self.__data__, int(val))

    @property
    def lowerqn(self):
        """ Lower state quantum numbers of the identifier (QuantumNumbers) """
        return QuantumNumbers(c.c_void_p(lib.getLowerQuantumNumbersQuantumIdentifier(self.__data__)))

    @lowerqn.setter
    def lowerqn(self, val):
        self.lowerqn.set(val)

    @property
    def upperqn(self):
        """ Upper state quantum numbers of the identifier (QuantumNumbers) """
        return QuantumNumbers(c.c_void_p(lib.getUpperQuantumNumbersQuantumIdentifier(self.__data__)))

    @upperqn.setter
    def upperqn(self, val):
        self.upperqn.set(val)

    @property
    def levelqn(self):
        return QuantumNumbers(c.c_void_p(lib.getEnergyLevelQuantumNumbersQuantumIdentifier(self.__data__)))

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
            lib.setSpeciesQuantumIdentifier(self.__data__, other.spec)
            lib.setIsotopologueQuantumIdentifier(self.__data__, other.isot)
            self.lowerqn = other.lowerqn
            self.upperqn = other.upperqn
            self.levelqn = other.levelqn  # repeat but for now keep...
        else:
            raise TypeError("Expects QuantumIdentifier")

    @staticmethod
    def name():
        return "QuantumIdentifier"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadQuantumIdentifier(self.__data__, correct_read_arguments(file)):
            raise OSError("Cannot read {}".format(file))

    def savexml(self, file, type="ascii", clobber=True):
        """ Saves the class to XML file

        Input:
            file:
                Filename to writable file (str)

            type:
                Filetype (str)

            clobber:
                Allow clobbering files? (any boolean)
        """
        if lib.xmlsaveQuantumIdentifier(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, QuantumIdentifier) and \
                self.type == other.type and \
                self.spec == other.spec and \
                self.isot == other.isot and \
                self.lowerqn == other.lowerqn and \
                self.upperqn == other.upperqn and \
                self.levelqn == other.levelqn:
            return True
        else:
            return False


# Generate ArrayOfQuantumIdentifier
exec(array_base(QuantumIdentifier))


lib.createQuantumIdentifier.restype = c.c_void_p
lib.createQuantumIdentifier.argtypes = []

lib.deleteQuantumIdentifier.restype = None
lib.deleteQuantumIdentifier.argtypes = [c.c_void_p]

lib.printQuantumIdentifier.restype = None
lib.printQuantumIdentifier.argtypes = [c.c_void_p]

lib.xmlreadQuantumIdentifier.restype = c.c_long
lib.xmlreadQuantumIdentifier.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveQuantumIdentifier.restype = c.c_long
lib.xmlsaveQuantumIdentifier.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getTypeQuantumIdentifier.restype = c.c_long
lib.getTypeQuantumIdentifier.argtypes = [c.c_void_p]

lib.setTypeQuantumIdentifier.restype = c.c_long
lib.setTypeQuantumIdentifier.argtypes = [c.c_void_p, c.c_long]

lib.string2indexTypeQuantumIdentifier.restype = c.c_long
lib.string2indexTypeQuantumIdentifier.argtypes = [c.c_void_p, c.c_char_p]

lib.getSpeciesQuantumIdentifier.restype = c.c_long
lib.getSpeciesQuantumIdentifier.argtypes = [c.c_void_p]

lib.setSpeciesQuantumIdentifier.restype = None
lib.setSpeciesQuantumIdentifier.argtypes = [c.c_void_p, c.c_long]

lib.getIsotopologueQuantumIdentifier.restype = c.c_long
lib.getIsotopologueQuantumIdentifier.argtypes = [c.c_void_p]

lib.setIsotopologueQuantumIdentifier.restype = None
lib.setIsotopologueQuantumIdentifier.argtypes = [c.c_void_p, c.c_long]

lib.getEnergyLevelQuantumNumbersQuantumIdentifier.restype = c.c_void_p
lib.getEnergyLevelQuantumNumbersQuantumIdentifier.argtypes = [c.c_void_p]

lib.getLowerQuantumNumbersQuantumIdentifier.restype = c.c_void_p
lib.getLowerQuantumNumbersQuantumIdentifier.argtypes = [c.c_void_p]

lib.getUpperQuantumNumbersQuantumIdentifier.restype = c.c_void_p
lib.getUpperQuantumNumbersQuantumIdentifier.argtypes = [c.c_void_p]
