import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class CovarianceMatrix:
    """ ARTS CovarianceMatrix data

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
            self.__data__ = c.c_void_p(lib.createCovarianceMatrix())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "CovarianceMatrix"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printCovarianceMatrix(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteCovarianceMatrix(self.__data__)

    def __repr__(self):
        return "ARTS CovarianceMatrix"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, CovarianceMatrix):
              raise RuntimeWarning("Cannot set CovarianceMatrix, remains constant")
        else:
            raise TypeError("Expects CovarianceMatrix")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadCovarianceMatrix(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveCovarianceMatrix(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createCovarianceMatrix.restype = c.c_void_p
lib.createCovarianceMatrix.argtypes = []

lib.deleteCovarianceMatrix.restype = None
lib.deleteCovarianceMatrix.argtypes = [c.c_void_p]

lib.printCovarianceMatrix.restype = None
lib.printCovarianceMatrix.argtypes = [c.c_void_p]

lib.xmlreadCovarianceMatrix.restype = c.c_long
lib.xmlreadCovarianceMatrix.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveCovarianceMatrix.restype = c.c_long
lib.xmlsaveCovarianceMatrix.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
