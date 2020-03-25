import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class Verbosity:
    """ ARTS Verbosity data

    Properties:
        data:
            The verbosity (agenda Index, screen Index, file Index, main bool)
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = False
            self.__data__ = c.c_void_p(lib.createVerbosity())
            if data is not None:
                self.data = data

    @staticmethod
    def name():
        return "Verbosity"

    @property
    def data(self, x):
        """ The verbosity (agenda Index, screen Index, file Index, main bool) """
        return lib.getAgendaVerbosity(self.__data__), lib.getScreenVerbosity(self.__data__), lib.getFileVerbosity(self.__data__), lib.getMainVerbosity(self.__data__)

    @data.setter
    def data(self, x):
        lib.setVerbosity(self.__data__, int(x[0]), int(x[1]), int(x[2]), bool(x[3]))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printVerbosity(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteVerbosity(self.__data__)

    def __repr__(self):
        return "ARTS Verbosity"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Verbosity):
            self.data = other.data
        else:
            raise TypeError("Expects Verbosity")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadVerbosity(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveVerbosity(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createVerbosity.restype = c.c_void_p
lib.createVerbosity.argtypes = []

lib.deleteVerbosity.restype = None
lib.deleteVerbosity.argtypes = [c.c_void_p]

lib.printVerbosity.restype = None
lib.printVerbosity.argtypes = [c.c_void_p]

lib.xmlreadVerbosity.restype = c.c_long
lib.xmlreadVerbosity.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveVerbosity.restype = c.c_long
lib.xmlsaveVerbosity.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getAgendaVerbosity.restype = c.c_long
lib.getAgendaVerbosity.argtypes = [c.c_void_p]

lib.getScreenVerbosity.restype = c.c_long
lib.getScreenVerbosity.argtypes = [c.c_void_p]

lib.getFileVerbosity.restype = c.c_long
lib.getFileVerbosity.argtypes = [c.c_void_p]

lib.getMainVerbosity.restype = c.c_long
lib.getMainVerbosity.argtypes = [c.c_void_p]

lib.setVerbosity.restype = None
lib.setVerbosity.argtypes = [c.c_void_p, c.c_long, c.c_long, c.c_long, c.c_bool]
