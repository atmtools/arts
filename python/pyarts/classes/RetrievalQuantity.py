import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class RetrievalQuantity:
    """ ARTS RetrievalQuantity data

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
            self.__data__ = c.c_void_p(lib.createRetrievalQuantity())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "RetrievalQuantity"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printRetrievalQuantity(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteRetrievalQuantity(self.__data__)

    def __repr__(self):
        return "ARTS RetrievalQuantity"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, RetrievalQuantity):
              raise RuntimeWarning("Cannot set RetrievalQuantity, remains constant")
        else:
            raise TypeError("Expects RetrievalQuantity")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadRetrievalQuantity(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveRetrievalQuantity(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(RetrievalQuantity))


lib.createRetrievalQuantity.restype = c.c_void_p
lib.createRetrievalQuantity.argtypes = []

lib.deleteRetrievalQuantity.restype = None
lib.deleteRetrievalQuantity.argtypes = [c.c_void_p]

lib.printRetrievalQuantity.restype = None
lib.printRetrievalQuantity.argtypes = [c.c_void_p]

lib.xmlreadRetrievalQuantity.restype = c.c_long
lib.xmlreadRetrievalQuantity.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveRetrievalQuantity.restype = c.c_long
lib.xmlsaveRetrievalQuantity.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
