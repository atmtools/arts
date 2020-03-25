import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class MCAntenna:
    """ ARTS MCAntenna data

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
            self.__data__ = c.c_void_p(lib.createMCAntenna())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "MCAntenna"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printMCAntenna(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteMCAntenna(self.__data__)

    def __repr__(self):
        return "ARTS MCAntenna"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, MCAntenna):
              raise RuntimeWarning("Cannot set MCAntenna, remains constant")
        else:
            raise TypeError("Expects MCAntenna")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadMCAntenna(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveMCAntenna(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createMCAntenna.restype = c.c_void_p
lib.createMCAntenna.argtypes = []

lib.deleteMCAntenna.restype = None
lib.deleteMCAntenna.argtypes = [c.c_void_p]

lib.printMCAntenna.restype = None
lib.printMCAntenna.argtypes = [c.c_void_p]

lib.xmlreadMCAntenna.restype = c.c_long
lib.xmlreadMCAntenna.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveMCAntenna.restype = c.c_long
lib.xmlsaveMCAntenna.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
