import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class Any:
    """ ARTS Any data

    Should not be used manually!

    Properties:
        None: Can only initialize from pointer and use XML
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createAny())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "Any"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAny(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAny(self.__data__)

    def __repr__(self):
        return "ARTS Any"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Any):
              raise RuntimeWarning("Cannot set Any, remains constant")
        else:
            raise TypeError("Expects Any")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadAny(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveAny(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createAny.restype = c.c_void_p
lib.createAny.argtypes = []

lib.deleteAny.restype = None
lib.deleteAny.argtypes = [c.c_void_p]

lib.printAny.restype = None
lib.printAny.argtypes = [c.c_void_p]

lib.xmlreadAny.restype = c.c_long
lib.xmlreadAny.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveAny.restype = c.c_long
lib.xmlsaveAny.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
