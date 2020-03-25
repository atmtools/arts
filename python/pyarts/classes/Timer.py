import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class Timer:
    """ ARTS Timer data

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
            self.__data__ = c.c_void_p(lib.createTimer())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "Timer"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTimer(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTimer(self.__data__)

    def __repr__(self):
        return "ARTS Timer"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Timer):
              raise RuntimeWarning("Cannot set Timer, remains constant")
        else:
            raise TypeError("Expects Timer")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTimer(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTimer(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createTimer.restype = c.c_void_p
lib.createTimer.argtypes = []

lib.deleteTimer.restype = None
lib.deleteTimer.argtypes = [c.c_void_p]

lib.printTimer.restype = None
lib.printTimer.argtypes = [c.c_void_p]

lib.xmlreadTimer.restype = c.c_long
lib.xmlreadTimer.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTimer.restype = c.c_long
lib.xmlsaveTimer.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
