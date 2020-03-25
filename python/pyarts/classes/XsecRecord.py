import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class XsecRecord:
    """ ARTS XsecRecord data

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
            self.__data__ = c.c_void_p(lib.createXsecRecord())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "XsecRecord"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printXsecRecord(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteXsecRecord(self.__data__)

    def __repr__(self):
        return "ARTS XsecRecord"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, XsecRecord):
              raise RuntimeWarning("Cannot set XsecRecord, remains constant")
        else:
            raise TypeError("Expects XsecRecord")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadXsecRecord(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveXsecRecord(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(XsecRecord))


lib.createXsecRecord.restype = c.c_void_p
lib.createXsecRecord.argtypes = []

lib.deleteXsecRecord.restype = None
lib.deleteXsecRecord.argtypes = [c.c_void_p]

lib.printXsecRecord.restype = None
lib.printXsecRecord.argtypes = [c.c_void_p]

lib.xmlreadXsecRecord.restype = c.c_long
lib.xmlreadXsecRecord.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveXsecRecord.restype = c.c_long
lib.xmlsaveXsecRecord.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
