import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np


class Tensor5:
    """ ARTS Tensor5 data

    Properties:
        dims:
            Dimensionality of Tensor5 (constexpr Index)

        data:
            The data (numpy-array)

        size:
            Shape of the Tensor5 (tuple)

    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTensor5())
            if data is not None:
                self.data = data

    @property
    def dims(self):
        return 5

    @property
    def data(self):
        """ The data (numpy-array) """
        if not np.any(self.shape == np.zeros(shape=(self.dims))):
            return np.ctypeslib.as_array(lib.getDataTensor5(self.__data__), self.shape)
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
        """ Shape of the Tensor5 (tuple) """
        return lib.shelvesTensor5(self.__data__), lib.booksTensor5(self.__data__), lib.pagesTensor5(self.__data__), lib.rowsTensor5(self.__data__), lib.colsTensor5(self.__data__),

    @shape.setter
    def shape(self, shape):
        if len(shape) == self.dims:
            lib.resizeTensor5(*shape, self.__data__)
        elif len(shape) < self.dims:
            shape = list(shape)
            while len(shape) < self.dims:
                shape.insert(0, 1)
            self.shape = shape
        else:
            raise RuntimeError("Cannot resize to larger dimensions")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTensor5(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTensor5(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Tensor5):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects Tensor5")

    @staticmethod
    def name():
        return "Tensor5"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTensor5(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTensor5(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


# ArrayOfTensor5
exec(array_base(Tensor5))


# ArrayOfArrayOfTensor5
# exec(array_base(ArrayOfTensor5))


lib.createTensor5.restype = c.c_void_p
lib.createTensor5.argtypes = []

lib.deleteTensor5.restype = None
lib.deleteTensor5.argtypes = [c.c_void_p]

lib.printTensor5.restype = None
lib.printTensor5.argtypes = [c.c_void_p]

lib.xmlreadTensor5.restype = c.c_long
lib.xmlreadTensor5.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTensor5.restype = c.c_long
lib.xmlsaveTensor5.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.colsTensor5.restype = c.c_long
lib.colsTensor5.argtypes = [c.c_void_p]

lib.rowsTensor5.restype = c.c_long
lib.rowsTensor5.argtypes = [c.c_void_p]

lib.pagesTensor5.restype = c.c_long
lib.pagesTensor5.argtypes = [c.c_void_p]

lib.booksTensor5.restype = c.c_long
lib.booksTensor5.argtypes = [c.c_void_p]

lib.shelvesTensor5.restype = c.c_long
lib.shelvesTensor5.argtypes = [c.c_void_p]

lib.resizeTensor5.restype = None
lib.resizeTensor5.argtypes = [c.c_long, c.c_long, c.c_long, c.c_long, c.c_long, c.c_void_p]

lib.getDataTensor5.restype = c.POINTER(c.c_double)
lib.getDataTensor5.argtypes = [c.c_void_p]
