import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np


class Tensor3:
    """ ARTS Tensor3 data

    Properties:
        dims:
            Dimensionality of Tensor3 (constexpr Index)

        data:
            The data (numpy-array)

        size:
            Shape of the Tensor3 (tuple)

    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTensor3())
            if data is not None:
                self.data = data

    @property
    def dims(self):
        return 3

    @property
    def data(self):
        """ The data (numpy-array) """
        if not np.any(self.shape == np.zeros(shape=(self.dims))):
            return np.ctypeslib.as_array(lib.getDataTensor3(self.__data__), self.shape)
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
        """ Shape of the Tensor3 (tuple) """
        return lib.pagesTensor3(self.__data__), lib.rowsTensor3(self.__data__), lib.colsTensor3(self.__data__),

    @shape.setter
    def shape(self, shape):
        if len(shape) == self.dims:
            lib.resizeTensor3(*shape, self.__data__)
        elif len(shape) < self.dims:
            shape = list(shape)
            while len(shape) < self.dims:
                shape.insert(0, 1)
            self.shape = shape
        else:
            raise RuntimeError("Cannot resize to larger dimensions")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTensor3(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTensor3(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Tensor3):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects Tensor3")

    @staticmethod
    def name():
        return "Tensor3"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTensor3(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTensor3(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


# ArrayOfTensor3
exec(array_base(Tensor3))


# ArrayOfArrayOfTensor3
exec(array_base(ArrayOfTensor3))


lib.createTensor3.restype = c.c_void_p
lib.createTensor3.argtypes = []

lib.deleteTensor3.restype = None
lib.deleteTensor3.argtypes = [c.c_void_p]

lib.printTensor3.restype = None
lib.printTensor3.argtypes = [c.c_void_p]

lib.xmlreadTensor3.restype = c.c_long
lib.xmlreadTensor3.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTensor3.restype = c.c_long
lib.xmlsaveTensor3.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.colsTensor3.restype = c.c_long
lib.colsTensor3.argtypes = [c.c_void_p]

lib.rowsTensor3.restype = c.c_long
lib.rowsTensor3.argtypes = [c.c_void_p]

lib.pagesTensor3.restype = c.c_long
lib.pagesTensor3.argtypes = [c.c_void_p]

lib.resizeTensor3.restype = None
lib.resizeTensor3.argtypes = [c.c_long, c.c_long, c.c_long, c.c_void_p]

lib.getDataTensor3.restype = c.POINTER(c.c_double)
lib.getDataTensor3.argtypes = [c.c_void_p]
