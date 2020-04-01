import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np
from copy import deepcopy as copy


class Matrix:
    """ ARTS Matrix data

    Properties:
        dims:
            Dimensionality of Matrix (constexpr Index)

        data:
            The data (numpy-array)

        size:
            Shape of the Matrix (tuple)

    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createMatrix())
            if data is not None:
                self.data = data

    @property
    def dims(self):
        return 2

    @property
    def flat(self):
        return self.data.flat

    @property
    def data(self):
        """ The data (numpy-array) """
        if not np.any(self.shape == np.zeros(shape=(self.dims))):
            return np.ctypeslib.as_array(lib.getDataMatrix(self.__data__), self.shape)
        else:
            return np.ndarray(self.shape, buffer=np.array([]))

    @data.setter
    def data(self, val):
        if not hasattr(val, 'shape'):
            val = np.array(val, dtype=float)

        if val.dtype != float:
            raise ValueError("Expects Numeric-like type. Got: {}".format(val.dtype))

        x = copy(val)
        self.shape = val.shape
        self.data.flat[:] = x.flat[:]

    @property
    def shape(self):
        """ Shape of the Matrix (tuple) """
        return lib.rowsMatrix(self.__data__), lib.colsMatrix(self.__data__),

    @shape.setter
    def shape(self, shape):
        if len(shape) == self.dims:
            lib.resizeMatrix(*shape, self.__data__)
        elif len(shape) < self.dims:
            shape = list(shape)
            while len(shape) < self.dims:
                shape.insert(0, 1)
            self.shape = shape
        else:
            raise RuntimeError("Cannot resize to larger dimensions")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printMatrix(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteMatrix(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Matrix):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects Matrix")

    @staticmethod
    def name():
        return "Matrix"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadMatrix(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveMatrix(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Matrix):
            return (self.data == other.data).all()
        else:
            return self.data == other

    def __lt__(self, other):
        if isinstance(other, Matrix):
            return (self.data < other.data).all()
        else:
            return self.data < other

    def __bool__(self):
        return bool((self.data != 0).any())

    def __array__(self):
        return self.data

    def __add__(self, val):
        return Matrix(self.data + np.array(val))

    def __sub__(self, val):
        return Matrix(self.data - np.array(val))

    def __mul__(self, val):
        return Matrix(self.data * np.array(val))

    def __matmul__(self, val):
        return self.data @ np.array(val)

    def __truediv__(self, val):
        return Matrix(self.data / np.array(val))

    def __pow__(self, val):
        return Matrix(self.data ** np.array(val))


# ArrayOfMatrix
exec(array_base(Matrix))


# ArrayOfArrayOfMatrix
exec(array_base(ArrayOfMatrix))


lib.createMatrix.restype = c.c_void_p
lib.createMatrix.argtypes = []

lib.deleteMatrix.restype = None
lib.deleteMatrix.argtypes = [c.c_void_p]

lib.printMatrix.restype = None
lib.printMatrix.argtypes = [c.c_void_p]

lib.xmlreadMatrix.restype = c.c_long
lib.xmlreadMatrix.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveMatrix.restype = c.c_long
lib.xmlsaveMatrix.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.colsMatrix.restype = c.c_long
lib.colsMatrix.argtypes = [c.c_void_p]

lib.rowsMatrix.restype = c.c_long
lib.rowsMatrix.argtypes = [c.c_void_p]

lib.resizeMatrix.restype = None
lib.resizeMatrix.argtypes = [c.c_long, c.c_long, c.c_void_p]

lib.getDataMatrix.restype = c.POINTER(c.c_double)
lib.getDataMatrix.argtypes = [c.c_void_p]
