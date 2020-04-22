import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class Time:
    """ ARTS Time data

    Note that the local epoch is UNIX time in CPP20 onwards

    Properties:
        sec:
            Time in system time seconds since local epoch (Numeric)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTime())
            if data is not None:
                self.sec = data

    @staticmethod
    def name():
        return "Time"

    @property
    def sec(self):
        """ Time in system time seconds since local epoch (Numeric) """
        return lib.getSecondsTime(self.__data__)

    @sec.setter
    def sec(self, x):
        return lib.setSecondsTime(self.__data__, float(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTime(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTime(self.__data__)

    def __repr__(self):
        return "ARTS Time"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Time):
            self.sec = other.sec
        else:
            raise TypeError("Expects Time")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTime(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTime(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Time) and self.sec == other.sec:
            return True
        else:
            return False


# ArrayOfTime
exec(array_base(Time))


lib.createTime.restype = c.c_void_p
lib.createTime.argtypes = []

lib.deleteTime.restype = None
lib.deleteTime.argtypes = [c.c_void_p]

lib.printTime.restype = None
lib.printTime.argtypes = [c.c_void_p]

lib.xmlreadTime.restype = c.c_long
lib.xmlreadTime.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTime.restype = c.c_long
lib.xmlsaveTime.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getSecondsTime.restype = c.c_double
lib.getSecondsTime.argtypes = [c.c_void_p]

lib.setSecondsTime.restype = None
lib.setSecondsTime.argtypes = [c.c_void_p, c.c_double]