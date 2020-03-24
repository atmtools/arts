import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.GriddedField1 import ArrayOfGriddedField1
from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class SpeciesAuxData:
    """ ARTS SpeciesAuxData data

    Properties:
        data:
            The data (const list of lists of ArrayOfGriddedField1)

        types:
            The types (const list of lists of Index)
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSpeciesAuxData())
            lib.initSpeciesAuxData(self.__data__)

    @property
    def data(self):
        """ The data (const list of lists of ArrayOfGriddedField1) """
        x = []
        s = 0
        while len(x) == 0 or len(x[-1]) != 0:
            i = 0
            x.append([])
            while lib.validindexSpeciesAuxData(self.__data__, int(s), int(i)):
                x[-1].append(ArrayOfGriddedField1(c.c_void_p(lib.getDataSpeciesAuxData(self.__data__, int(s), int(i)))))
                i += 1
            s += 1
        x.pop()
        return x

    @property
    def types(self):
        """ The types (const list of lists of Index) """
        x = []
        s = 0
        while len(x) == 0 or len(x[-1]) != 0:
            i = 0
            x.append([])
            while lib.validindexSpeciesAuxData(self.__data__, int(s), int(i)):
                x[-1].append(lib.getTypeSpeciesAuxData(self.__data__, int(s), int(i)))
                i += 1
            s += 1
        return x

    def setType(self, s, i, t):
        """ Sets type of aux data

        Only for expert users

        Input:
            s:
                Species index (int)

            i:
                Isotopologue index (int)

            t:
                Type index (int)
        """
        if lib.validindexSpeciesAuxData(self.__data__, int(s), int(i)):
            if lib.setTypeFromIndexSpeciesAuxData(self.__data__, int(s), int(i), int(t)):
                raise TypeError("Invalid type")
        else:
            raise IndexError("Out of bounds")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSpeciesAuxData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesAuxData(self.__data__)

    def __repr__(self):
        return "ARTS SpeciesAuxData"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        raise TypeError("Cannot set full SpeciesAuxData")

    @staticmethod
    def name():
        return "SpeciesAuxData"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadSpeciesAuxData(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveSpeciesAuxData(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createSpeciesAuxData.restype = c.c_void_p
lib.createSpeciesAuxData.argtypes = []

lib.deleteSpeciesAuxData.restype = None
lib.deleteSpeciesAuxData.argtypes = [c.c_void_p]

lib.printSpeciesAuxData.restype = None
lib.printSpeciesAuxData.argtypes = [c.c_void_p]

lib.xmlreadSpeciesAuxData.restype = c.c_long
lib.xmlreadSpeciesAuxData.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveSpeciesAuxData.restype = c.c_long
lib.xmlsaveSpeciesAuxData.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.initSpeciesAuxData.restype = None
lib.initSpeciesAuxData.argtypes = [c.c_void_p]

lib.validindexSpeciesAuxData.restype = c.c_bool
lib.validindexSpeciesAuxData.argtypes = [c.c_void_p, c.c_long, c.c_long]

lib.getDataSpeciesAuxData.restype = c.c_void_p
lib.getDataSpeciesAuxData.argtypes = [c.c_void_p, c.c_long, c.c_long]

lib.setTypeFromIndexSpeciesAuxData.restype = c.c_long
lib.setTypeFromIndexSpeciesAuxData.argtypes = [c.c_void_p, c.c_long, c.c_long, c.c_long]

lib.getTypeSpeciesAuxData.restype = c.c_long
lib.getTypeSpeciesAuxData.argtypes = [c.c_void_p, c.c_long, c.c_long]
