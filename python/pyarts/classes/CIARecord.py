import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes.GriddedField2 import ArrayOfGriddedField2
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class CIARecord:
    """ ARTS CIARecord data

    Properties:
        specs:
            The species (tuple of 2 Index)

        data:
            The data (ArrayOfGriddedField2)
    """
    def __init__(self, data=None, specs=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = False
            self.__data__ = c.c_void_p(lib.createCIARecord())
            if data is not None:
                self.data = data
                self.specs = specs

    @staticmethod
    def name():
        return "CIARecord"

    @property
    def data(self):
        """ The data (ArrayOfGriddedField2) """
        return ArrayOfGriddedField2(c.c_void_p(lib.getDataCIARecord(self.__data__)))

    @data.setter
    def data(self, x):
        self.data.set(x)

    @property
    def specs(self):
        """ The species (tuple of 2 Index) """
        return lib.getSpecies1CIARecord(self.__data__), lib.getSpecies2CIARecord(self.__data__)

    @specs.setter
    def specs(self, x):
        if isinstance(x, Sized) and len(x) == 2:
            if SpeciesTag.validSpecies(x[0]) and SpeciesTag.validSpecies(x[1]):
                lib.setSpeciesCIARecord(self.__data__, int(x[0]), int(x[1]))
            else:
                raise ValueError("Bad species")
        else:
            raise TypeError("Expected len-2 array-like")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printCIARecord(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteCIARecord(self.__data__)

    def __repr__(self):
        return "ARTS CIARecord"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, CIARecord):
            self.specs = other.specs
            self.data = other.data
        else:
            raise TypeError("Expects CIARecord")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadCIARecord(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveCIARecord(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, CIARecord) and self.specs == other.specs and self.data == other.data:
            return True
        else:
            return False



exec(array_base(CIARecord))


lib.createCIARecord.restype = c.c_void_p
lib.createCIARecord.argtypes = []

lib.deleteCIARecord.restype = None
lib.deleteCIARecord.argtypes = [c.c_void_p]

lib.printCIARecord.restype = None
lib.printCIARecord.argtypes = [c.c_void_p]

lib.xmlreadCIARecord.restype = c.c_long
lib.xmlreadCIARecord.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveCIARecord.restype = c.c_long
lib.xmlsaveCIARecord.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getSpecies1CIARecord.restype = c.c_long
lib.getSpecies1CIARecord.argtypes = [c.c_void_p]

lib.getSpecies2CIARecord.restype = c.c_long
lib.getSpecies2CIARecord.argtypes = [c.c_void_p]

lib.setSpeciesCIARecord.restype = None
lib.setSpeciesCIARecord.argtypes = [c.c_void_p, c.c_long, c.c_long]

lib.getDataCIARecord.restype = c.c_void_p
lib.getDataCIARecord.argtypes = [c.c_void_p]
