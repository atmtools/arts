import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class TelsemAtlas:
    """ ARTS TelsemAtlas data

    FIXME:  NO INTERFACE AVAILABLE

    Properties:
        None: Can only initialize from pointer and use XML
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTelsemAtlas())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "TelsemAtlas"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTelsemAtlas(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTelsemAtlas(self.__data__)

    def __repr__(self):
        return "ARTS TelsemAtlas"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, TelsemAtlas):
              raise RuntimeWarning("Cannot set TelsemAtlas, remains constant")
        else:
            raise TypeError("Expects TelsemAtlas")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTelsemAtlas(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTelsemAtlas(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(TelsemAtlas))


lib.createTelsemAtlas.restype = c.c_void_p
lib.createTelsemAtlas.argtypes = []

lib.deleteTelsemAtlas.restype = None
lib.deleteTelsemAtlas.argtypes = [c.c_void_p]

lib.printTelsemAtlas.restype = None
lib.printTelsemAtlas.argtypes = [c.c_void_p]

lib.xmlreadTelsemAtlas.restype = c.c_long
lib.xmlreadTelsemAtlas.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTelsemAtlas.restype = c.c_long
lib.xmlsaveTelsemAtlas.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
