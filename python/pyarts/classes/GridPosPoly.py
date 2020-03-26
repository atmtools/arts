import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import ArrayOfIndex
from pyarts.classes.Vector import Vector
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class GridPosPoly:
    """ ARTS GridPosPoly data

    Properties:
        idx:
            Index of position below point (ArrayOfIndex)

        w:
            Weight of grids (Vector)
        """
    def __init__(self, idx=0, fd=[0, 0]):
        if isinstance(type, c.c_void_p):
            self.__data__ = idx
        else:
            self.__data__ = lib.createGridPosPoly()
            self.idx = idx
            self.fd = fd

    @staticmethod
    def name():
        return "GridPosPoly"

    @property
    def idx(self):
        """ Index of position below point (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getidxGridPosPoly(self.__data__)))

    @idx.setter
    def idx(self, val):
        self.idx.set(val)

    @property
    def w(self):
        """ Weight of grids (Vector) """
        return Vector(c.c_void_p(lib.getwGridPosPoly(self.__data__)))

    @w.setter
    def w(self, val):
        self.w.set(val)

    def __del__(self):
        if self.__delete__:
            lib.deleteGridPosPoly(self.__data__)

    def __repr__(self):
        return "ARTS GridPosPoly"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, GridPosPoly):
            self.idx = other.idx
            self.w = other.w
        else:
            raise TypeError("Expects GridPosPoly")


exec(array_base(GridPosPoly))


lib.createGridPosPoly.restype = c.c_void_p
lib.createGridPosPoly.argtypes = []

lib.deleteGridPosPoly.restype = None
lib.deleteGridPosPoly.argtypes = [c.c_void_p]

lib.printGridPosPoly.restype = None
lib.printGridPosPoly.argtypes = [c.c_void_p]

lib.getidxGridPosPoly.restype = c.c_void_p
lib.getidxGridPosPoly.argtypes = [c.c_void_p]

lib.getwGridPosPoly.restype = c.c_void_p
lib.getwGridPosPoly.argtypes = [c.c_void_p]
