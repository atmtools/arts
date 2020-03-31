import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np


class Vector:
    """ ARTS Vector data

    Properties:
        dims:
            Dimensionality of Vector (constexpr Index)

        data:
            The data (numpy-array)

        size:
            Shape of the Vector (tuple)

    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createVector())
            if data is not None:
                self.data = data

    @property
    def dims(self):
        return 1

    @property
    def data(self):
        """ The data (numpy-array) """
        if not np.any(self.shape == np.zeros(shape=(self.dims))):
            return np.ctypeslib.as_array(lib.getDataVector(self.__data__), self.shape)
        else:
            return np.ndarray(self.shape, buffer=np.array([]))

    @data.setter
    def data(self, val):
        if not hasattr(val, 'shape'):
            val = np.array(val, dtype=float)

        if val.dtype != float:
            raise ValueError("Expects Numeric-like type. Got: {}".format(val.dtype))

        self.shape = val.shape
        self.data.flat[:] = val.flat[:]

    @property
    def shape(self):
        """ Shape of the Vector (tuple) """
        return lib.nelemVector(self.__data__),

    @shape.setter
    def shape(self, shape):
        if isinstance(shape, int):
            self.shape = shape,
        elif len(shape) == self.dims:
            lib.resizeVector(*shape, self.__data__)
        elif len(shape) < self.dims:
            shape = list(shape)
            while len(shape) < self.dims:
                shape.insert(0, 1)
            self.shape = shape
        else:
            raise RuntimeError("Cannot resize to larger dimensions")


    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printVector(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteVector(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Vector):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects Vector")

    @staticmethod
    def name():
        return "Vector"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadVector(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveVector(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Vector):
            return (self.data == other.data).all()
        else:
            return self.data == other

    def __lt__(self, other):
        if isinstance(other, Vector):
            return (self.data < other.data).all()
        else:
            return self.data < other

    def __bool__(self):
        return bool((self.data != 0).any())


# ArrayOfVector
exec(array_base(Vector))


# ArrayOfArrayOfVector
exec(array_base(ArrayOfVector))


lib.createVector.restype = c.c_void_p
lib.createVector.argtypes = []

lib.deleteVector.restype = None
lib.deleteVector.argtypes = [c.c_void_p]

lib.printVector.restype = None
lib.printVector.argtypes = [c.c_void_p]

lib.xmlreadVector.restype = c.c_long
lib.xmlreadVector.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveVector.restype = c.c_long
lib.xmlsaveVector.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.nelemVector.restype = c.c_long
lib.nelemVector.argtypes = [c.c_void_p]

lib.resizeVector.restype = None
lib.resizeVector.argtypes = [c.c_long, c.c_void_p]

lib.getDataVector.restype = c.POINTER(c.c_double)
lib.getDataVector.argtypes = [c.c_void_p]
