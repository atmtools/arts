import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.LineShapeSingleSpeciesModel import LineShapeSingleSpeciesModel


class LineShapeModel:
    """ ARTS LineShape::Model data

    Getter and setter for data-property is implemented

    Properties:
        size:
            Number of LineShapeSingleSpeciesModel(s) (Index)

        data:
            The data (list of LineShapeSingleSpeciesModel(s))
        """
    def __init__(self, data=[]):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createLineShapeModel())
            self.data = data

    @property
    def size(self):
        """ Number of LineShapeSingleSpeciesModel(s) (Index) """
        return lib.sizeLineShapeModel(self.__data__)
    @size.setter
    def size(self, size):
        size = int(size)
        if size < 0:
            raise ValueError("size must be 0 or larger")
        lib.resizeLineShapeModel(size, self.__data__)

    @property
    def data(self):
        """ The data (list of LineShapeSingleSpeciesModel(s)) """
        x = []
        n = self.size
        for i in range(n):
            x.append(self[i])
        return x

    @data.setter
    def data(self, val):
        if isinstance(val, Sized):
            self.size = len(val)
            n = self.size
            for i in range(n):
                self[i] = val[i]
        else:
            raise TypeError("Only accepts array-like input")

    def __getitem__(self, ind):
        if ind >= 0 and ind < self.size:
            return LineShapeSingleSpeciesModel(c.c_void_p(lib.getelemLineShapeModel(ind, self.__data__)))
        else:
            raise IndexError("Out of bounds")

    def __setitem__(self, ind, val):
        self[ind].set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printLineShapeModel(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteLineShapeModel(self.__data__)

    def __repr__(self):
        return "ARTS LineShape::Model with {} LineShape::SingleSpeciesModel(s)".format(self.size)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, LineShapeModel):
            self.data = other.data
        else:
            raise TypeError("Expects LineShapeModel")


lib.createLineShapeModel.restype = c.c_void_p
lib.createLineShapeModel.argtypes = []

lib.deleteLineShapeModel.restype = None
lib.deleteLineShapeModel.argtypes = [c.c_void_p]

lib.printLineShapeModel.restype = None
lib.printLineShapeModel.argtypes = [c.c_void_p]

lib.getelemLineShapeModel.restype = c.c_void_p
lib.getelemLineShapeModel.argtypes = [c.c_long, c.c_void_p]

lib.sizeLineShapeModel.restype = c.c_long
lib.sizeLineShapeModel.argtypes = []

lib.resizeLineShapeModel.restype = None
lib.resizeLineShapeModel.argtypes = [c.c_long]
