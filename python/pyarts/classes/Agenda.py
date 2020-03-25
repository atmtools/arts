import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class Agenda:
    """ ARTS Agenda data

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
            self.__data__ = c.c_void_p(lib.createAgenda())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "Agenda"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAgenda(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAgenda(self.__data__)

    def __repr__(self):
        return "ARTS Agenda"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Agenda):
              raise RuntimeWarning("Cannot set Agenda, remains constant")
        else:
            raise TypeError("Expects Agenda")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadAgenda(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveAgenda(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(Agenda))


lib.createAgenda.restype = c.c_void_p
lib.createAgenda.argtypes = []

lib.deleteAgenda.restype = None
lib.deleteAgenda.argtypes = [c.c_void_p]

lib.printAgenda.restype = None
lib.printAgenda.argtypes = [c.c_void_p]

lib.xmlreadAgenda.restype = c.c_long
lib.xmlreadAgenda.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveAgenda.restype = c.c_long
lib.xmlsaveAgenda.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
