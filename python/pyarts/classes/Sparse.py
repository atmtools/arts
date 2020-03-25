import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class Sparse:
    """ ARTS Sparse data

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
            self.__data__ = c.c_void_p(lib.createSparse())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "Sparse"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSparse(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSparse(self.__data__)

    def __repr__(self):
        return "ARTS Sparse"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Sparse):
              raise RuntimeWarning("Cannot set Sparse, remains constant")
        else:
            raise TypeError("Expects Sparse")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadSparse(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveSparse(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(Sparse))


lib.createSparse.restype = c.c_void_p
lib.createSparse.argtypes = []

lib.deleteSparse.restype = None
lib.deleteSparse.argtypes = [c.c_void_p]

lib.printSparse.restype = None
lib.printSparse.argtypes = [c.c_void_p]

lib.xmlreadSparse.restype = c.c_long
lib.xmlreadSparse.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveSparse.restype = c.c_long
lib.xmlsaveSparse.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
