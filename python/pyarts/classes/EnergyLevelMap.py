import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.QuantumIdentifier import ArrayOfQuantumIdentifier
from pyarts.classes.Vector import Vector
from pyarts.classes.Tensor4 import Tensor4
from pyarts.classes.io import correct_save_arguments, correct_read_arguments

class EnergyLevelMap:
    """ ARTS EnergyLevelMap data

    Properties:
        type:
            Type of the EnergyLevelMap (Index)

        levels:
            Identification (ArrayOfQuantumIdentifier)

        energies:
            Energy data (Vector)

        data:
            NLTE data (Tensor4)
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createEnergyLevelMap())
            if data is not None:
                raise ValueError("Cannot init with value")

    @property
    def OK(self):
        return bool(lib.getOKEnergyLevelMap(self.__data__))

    @property
    def type(self):
        """ Type of the EnergyLevelMap (Index) """
        return lib.getTypeEnergyLevelMap(self.__data__)

    @type.setter
    def type(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexTypeEnergyLevelMap(self.__data__, type.encode("ascii")))
        elif lib.setTypeEnergyLevelMap(self.__data__, int(type)):
            raise ValueError("Invalid type")

    @property
    def levels(self):
        """ Identification (ArrayOfQuantumIdentifier) """
        return ArrayOfQuantumIdentifier(c.c_void_p(lib.getLevelsEnergyLevelMap(self.__data__)))

    @levels.setter
    def levels(self, val):
        self.levels.set(val)

    @property
    def energies(self):
        """ Energy data (Vector) """
        return Vector(c.c_void_p(lib.getEnergiesEnergyLevelMap(self.__data__)))

    @energies.setter
    def energies(self, val):
        self.energies.set(val)

    @property
    def data(self):
        """ NLTE data (Tensor4) """
        return Tensor4(c.c_void_p(lib.getDataEnergyLevelMap(self.__data__)))

    @data.setter
    def data(self, val):
        self.data.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        if self.OK:
            lib.printEnergyLevelMap(self.__data__)
        else:
            raise RuntimeError("Class in bad state")

    def __del__(self):
        if self.__delete__:
            lib.deleteEnergyLevelMap(self.__data__)

    def __repr__(self):
        return "ARTS EnergyLevelMap"

    @staticmethod
    def name():
        """ Name of the class as a string.  Required for arrayification """
        return "EnergyLevelMap"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, EnergyLevelMap):
            self.type = other.type
            self.levels = other.levels
            self.energies = other.energies
            self.data = other.data
        else:
            raise TypeError("Expects EnergyLevelMap")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadEnergyLevelMap(self.__data__, correct_read_arguments(file)):
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
        if not self.OK:
            raise RuntimeError("Class in bad state")

        if lib.xmlsaveEnergyLevelMap(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, EnergyLevelMap) and \
                self.type == other.type and \
                self.levels == other.levels and \
                self.energies == other.energies and \
                self.data == other.data:
            return True
        else:
            return False

    def __bool__(self):
        return self.OK and bool(self.data)


lib.createEnergyLevelMap.restype = c.c_void_p
lib.createEnergyLevelMap.argtypes = []

lib.deleteEnergyLevelMap.restype = None
lib.deleteEnergyLevelMap.argtypes = [c.c_void_p]

lib.printEnergyLevelMap.restype = None
lib.printEnergyLevelMap.argtypes = [c.c_void_p]

lib.xmlreadEnergyLevelMap.restype = c.c_long
lib.xmlreadEnergyLevelMap.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveEnergyLevelMap.restype = c.c_long
lib.xmlsaveEnergyLevelMap.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getTypeEnergyLevelMap.restype = c.c_long
lib.getTypeEnergyLevelMap.argtypes = [c.c_void_p]

lib.setTypeEnergyLevelMap.restype = c.c_long
lib.setTypeEnergyLevelMap.argtypes = [c.c_void_p, c.c_long]

lib.string2indexTypeEnergyLevelMap.restype = c.c_long
lib.string2indexTypeEnergyLevelMap.argtypes = [c.c_void_p, c.c_char_p]

lib.getLevelsEnergyLevelMap.restype = c.c_void_p
lib.getLevelsEnergyLevelMap.argtypes = [c.c_void_p]

lib.getEnergiesEnergyLevelMap.restype = c.c_void_p
lib.getEnergiesEnergyLevelMap.argtypes = [c.c_void_p]

lib.getDataEnergyLevelMap.restype = c.c_void_p
lib.getDataEnergyLevelMap.argtypes = [c.c_void_p]

lib.getOKEnergyLevelMap.restype = c.c_bool
lib.getOKEnergyLevelMap.argtypes = [c.c_void_p]
