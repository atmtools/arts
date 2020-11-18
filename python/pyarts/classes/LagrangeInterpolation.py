import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Index, ArrayOfNumeric
from pyarts.classes.Vector import Vector
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class LagrangeInterpolation:
    """ ARTS LagrangeInterpolation data

    Properties:
        pos:
            Index of position below point (fIndex)

        lx:
            Weight of grids (Array<Numeric>)

        dlx:
            Weight of grids derivatives (Array<Numeric>)
        """
    def __init__(self, pos=0, lx=[1], dlx=[]):
        if isinstance(type, c.c_void_p):
            self.__data__ = pos
        else:
            self.__data__ = lib.createLagrangeInterpolation()
            self.pos = pos
            self.lx = lx
            self.dlx = dlx

    @staticmethod
    def name():
        return "LagrangeInterpolation"

    @property
    def pos(self):
        """ Index of position below point (ArrayOfIndex) """
        return Index(c.c_void_p(lib.getposLagrangeInterpolation(self.__data__)))

    @pos.setter
    def pos(self, val):
        self.pos.set(val)

    @property
    def lx(self):
        """ Weight of grids (Vector) """
        return ArrayOfNumeric(c.c_void_p(lib.getlxLagrangeInterpolation(self.__data__)))

    @lx.setter
    def lx(self, val):
        self.lx.set(val)

    @property
    def dlx(self):
        """ Weight of grids (Vector) """
        return ArrayOfNumeric(c.c_void_p(lib.getdlxLagrangeInterpolation(self.__data__)))

    @dlx.setter
    def dlx(self, val):
        self.dlx.set(val)

    def __del__(self):
        if self.__delete__:
            lib.deleteLagrangeInterpolation(self.__data__)

    def __repr__(self):
        return "ARTS LagrangeInterpolation"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printLagrangeInterpolation(self.__data__)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, LagrangeInterpolation):
            self.pos = other.pos
            self.lx = other.lx
            self.dlx = other.dlx
        else:
            raise TypeError("Expects LagrangeInterpolation")

    def __eq__(self, other):
        if isinstance(other, LagrangeInterpolation) and self.pos == other.pos and self.lx == other.lx and self.dlx == other.dlx:
            return True
        else:
            return False


exec(array_base(LagrangeInterpolation))


lib.createLagrangeInterpolation.restype = c.c_void_p
lib.createLagrangeInterpolation.argtypes = []

lib.deleteLagrangeInterpolation.restype = None
lib.deleteLagrangeInterpolation.argtypes = [c.c_void_p]

lib.printLagrangeInterpolation.restype = None
lib.printLagrangeInterpolation.argtypes = [c.c_void_p]

lib.getposLagrangeInterpolation.restype = c.c_void_p
lib.getposLagrangeInterpolation.argtypes = [c.c_void_p]

lib.getlxLagrangeInterpolation.restype = c.c_void_p
lib.getlxLagrangeInterpolation.argtypes = [c.c_void_p]

lib.getdlxLagrangeInterpolation.restype = c.c_void_p
lib.getdlxLagrangeInterpolation.argtypes = [c.c_void_p]
