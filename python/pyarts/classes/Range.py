import ctypes as c
from pyarts.workspace.api import arts_api as lib


class Range:
    """ ARTS Range data

    Properties:
        start:
            Returns the start index of the range (Index)

        extent:
            Returns the extent of the range (Index)

        stride:
            Returns the stride of the range (Index)
    """
    def __init__(self, start=0, extent=-1, stride=1):
        if isinstance(start, c.c_void_p):
            self.__delete__ = False
            self.__data__ = start
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createRange())

            assert start >= 0
            lib.setRange(self.__data__, int(start), int(extent), int(stride))

    @property
    def stop(self):
        """ Python sice-like stop """
        x = (self.start + self.stride * self.extent)
        if self.stride > 0:
            return x if self.extent >= 0 else None
        else:
            return x if self.extent >= 0 and x >= 0 else None

    @property
    def start(self):
        """ Returns the start index of the range (Index) """
        return lib.get_startRange(self.__data__)

    @property
    def extent(self):
        """ Returns the extent of the range (Index) """
        return lib.get_extentRange(self.__data__)

    @property
    def stride(self):
        """ Returns the stride of the range (Index) """
        return lib.get_strideRange(self.__data__)

    @start.setter
    def start(self, x):
        assert x >= 0
        lib.setRange(self.__data__, int(x), self.extent, self.stride)

    @extent.setter
    def extent(self, x):
        lib.setRange(self.__data__, self.start, int(x), self.stride)

    @stride.setter
    def stride(self, x):
        lib.setRange(self.__data__, self.start, self.extent, int(x))

    @staticmethod
    def name():
        return "Range"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printRange(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteRange(self.__data__)

    def __repr__(self):
        return "ARTS Range"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Range):
            lib.setRange(self.__data__, other.start, other.extent, other.stride)
        else:
            raise TypeError("Expects Range")

    def to_slice(self):
        """ Returns a python slice object from ARTS Range

        Warnings: May return bad extent for most python classes; may return
        bad start for most python classes
        """
        return slice(self.start, self.stop, self.stride)


lib.createRange.restype = c.c_void_p
lib.createRange.argtypes = []

lib.deleteRange.restype = None
lib.deleteRange.argtypes = [c.c_void_p]

lib.printRange.restype = None
lib.printRange.argtypes = [c.c_void_p]

lib.get_startRange.restype = c.c_long
lib.get_startRange.argtypes = [c.c_void_p]

lib.get_strideRange.restype = c.c_long
lib.get_strideRange.argtypes = [c.c_void_p]

lib.get_extentRange.restype = c.c_long
lib.get_extentRange.argtypes = [c.c_void_p]

lib.setRange.restype = None
lib.setRange.argtypes = [c.c_void_p, c.c_long, c.c_long, c.c_long]
