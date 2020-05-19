import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

from datetime import datetime
import numpy as np


class Time:
    """ ARTS Time data

    Note that the local epoch is UNIX time in CPP20 onwards, but that we assume it
    is the time on your system for this class to work.  Please confirm this before
    using any Time features

    Addition to class is treated like adding seconds.  Subtraction from class
    is also treated as seconds.  Subtracting Time from another time gives the
    number of seconds between two time stamps.  Note that since the internal
    clock is in units of scaled integers, even adding 0 to Time may round the
    clock enough that (t0 + 0 != t0) is True.

    Properties:
        sec:
            Time in system time seconds since local epoch (not-an-exact Numeric)
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

    @staticmethod
    def now():
        return Time()

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
        return self.strftime("%Y-%m-%d %H:%M:%S.%f")

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Time):
            lib.setTime(self.__data__, other.__data__)  # Workaround that "sec" is not exact
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

    def netcdfify(self):
        """ Create the NETCDF4 information required for writing this data
        
            Output: list that can be processed by netcdf.py, False arraytype
        """
        return [["val", str(self), str, {}]], False

    def denetcdf(self, group):
        """ Create the NETCDF4 information required for writing this data
        
            Output: list that can be processed by netcdf.py, False arraytype
        """
        self.sec = datetime.fromisoformat(group.val).timestamp()

    def __eq__(self, other):
        if isinstance(other, Time) and lib.equalTime(self.__data__, other.__data__):
            return True
        else:
            return False

    def __lt__(self, other):
        if isinstance(other, Time) and lib.lessTime(self.__data__, other.__data__):
            return True
        else:
            return False

    def __add__(self, other):
        if isinstance(other, Time):
            raise TypeError("Cannot add Time to itself")
        else:
            return Time(self.sec + other)

    def __sub__(self, other):
        if isinstance(other, Time):
            return self.sec - other.sec
        else:
            return Time(self.sec-other)

    def __radd__(self, other):
        if isinstance(other, Time):
            raise TypeError("Cannot add Time to itself")
        else:
            return Time(other + self.sec)

    def __rsub__(self, other):
        if isinstance(other, Time):
            return other.sec - self.sec
        else:
            return Time(self.sec-other)
    
    def strftime(self, fmt : str):
        """ Return self as string formatted by datetime """
        return datetime.fromtimestamp(self.sec).strftime(fmt)
    
    def fromisoformat(self, time: str):
        """ Sets self.sec from string by datetime fromisoformat """
        self.sec = datetime.fromisoformat(time).timestamp()


# ArrayOfTime
exec(array_base(Time))


def TimeGrid(ts):
    """ Return an array of datetimes from an ArrayOfTime or list of Time

    Helps to make calls to matplotlib look decent

    Example:
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> ts = [Time.now() + 3600*i for i in range(60*24)]  # Hourly time stamps for 60 days from now
        >>> y = [np.sin(4*np.pi * t.sec / (3600*60*24)) for t in ts]  # 2 Waves
        >>> plt.plot(TimeGrid(ts), y), plt.gcf().autofmt_xdate()  # Formatted plots

    Input:
        ts: An ArrayOfTime

    Output:
        List of datetime.datetime(s)
    """
    return np.array([datetime.fromtimestamp(t.sec) for t in ts])


# Makes example above work directly on ArrayOfTime without TimeGrid-conversion
ArrayOfTime.__array__ = TimeGrid


# ArrayOfArrayOfTime
exec(array_base(ArrayOfTime))


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

lib.setTime.restype = None
lib.setTime.argtypes = [c.c_void_p, c.c_void_p]

lib.equalTime.restype = c.c_bool
lib.equalTime.argtypes = [c.c_void_p, c.c_void_p]

lib.lessTime.restype = c.c_bool
lib.lessTime.argtypes = [c.c_void_p, c.c_void_p]

