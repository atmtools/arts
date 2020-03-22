import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class Index:
    """ ARTS Index data

    Properties:
        val:
            a value (int)
    """
    def __init__(self, value):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = False
            self.__data__ = lib.createIndex()
            self.val = value

    @staticmethod
    def name():
        return "Index"

    @property
    def val(self):
        """ a value (int) """
        return lib.getIndex(self.__data__)

    @val.setter
    def val(self, x):
        lib.setIndex(self.__data__, int(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printIndex(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteIndex(self.__data__)

    def __repr__(self):
        return "{}".format(self.val)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Index):
            self.val = other.val
        else:
            self.val = other

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadIndex(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveIndex(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


class Numeric:
    """ ARTS Numeric data

    Properties:
        val:
            a value (int)
    """
    def __init__(self, value):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = False
            self.__data__ = lib.createNumeric()
            self.val = value

    @staticmethod
    def name():
        return "Numeric"

    @property
    def val(self):
        """ a value (float) """
        return lib.getNumeric(self.__data__)

    @val.setter
    def val(self, x):
        lib.setNumeric(self.__data__, float(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printNumeric(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteNumeric(self.__data__)

    def __repr__(self):
        return "{}".format(self.val)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Numeric):
            self.val = other.val
        else:
            self.val = other

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadNumeric(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveNumeric(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(Index))

exec(array_base(ArrayOfIndex))


lib.createIndex.restype = c.c_void_p
lib.createIndex.argtypes = []

lib.deleteIndex.restype = None
lib.deleteIndex.argtypes = [c.c_void_p]

lib.printIndex.restype = None
lib.printIndex.argtypes = [c.c_void_p]

lib.getIndex.restype = c.c_long
lib.getIndex.argtypes = [c.c_void_p]

lib.setIndex.restype = None
lib.setIndex.argtypes = [c.c_void_p, c.c_long]

lib.xmlreadIndex.restype = c.c_long
lib.xmlreadIndex.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveIndex.restype = c.c_long
lib.xmlsaveIndex.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]



lib.createNumeric.restype = c.c_void_p
lib.createNumeric.argtypes = []

lib.deleteNumeric.restype = None
lib.deleteNumeric.argtypes = [c.c_void_p]

lib.printNumeric.restype = None
lib.printNumeric.argtypes = [c.c_void_p]

lib.getNumeric.restype = c.c_double
lib.getNumeric.argtypes = [c.c_void_p]

lib.setNumeric.restype = None
lib.setNumeric.argtypes = [c.c_void_p, c.c_double]

lib.xmlreadNumeric.restype = c.c_long
lib.xmlreadNumeric.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveNumeric.restype = c.c_long
lib.xmlsaveNumeric.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
