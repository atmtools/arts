import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np
from copy import deepcopy as copy


class Tensor6:
    """ ARTS Tensor6 data

    Properties:
        dims:
            Dimensionality of Tensor6 (constexpr Index)

        data:
            The data (numpy-array)

        size:
            Shape of the Tensor6 (tuple)

    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTensor6())
            if data is not None:
                self.data = data

    @property
    def dims(self):
        return 6

    @property
    def flat(self):
        return self.data.flat

    @property
    def data(self):
        """ The data (numpy-array) """
        if not np.any(self.shape == np.zeros(shape=(self.dims))):
            return np.ctypeslib.as_array(lib.getDataTensor6(self.__data__), self.shape)
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
        """ Shape of the Tensor6 (tuple) """
        return lib.vitrinesTensor6(self.__data__), lib.shelvesTensor6(self.__data__), lib.booksTensor6(self.__data__), lib.pagesTensor6(self.__data__), lib.rowsTensor6(self.__data__), lib.colsTensor6(self.__data__),

    @shape.setter
    def shape(self, shape):
        if len(shape) == self.dims:
            lib.resizeTensor6(*shape, self.__data__)
        elif len(shape) < self.dims:
            shape = list(shape)
            while len(shape) < self.dims:
                shape.insert(0, 1)
            self.shape = shape
        else:
            raise RuntimeError("Cannot resize to larger dimensions")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTensor6(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTensor6(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Tensor6):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects Tensor6")

    @staticmethod
    def name():
        return "Tensor6"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTensor6(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTensor6(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Tensor6):
            return (self.data == other.data).all()
        else:
            return self.data == other

    def __lt__(self, other):
        if isinstance(other, Tensor6):
            return (self.data < other.data).all()
        else:
            return self.data < other

    def __bool__(self):
        return bool((self.data != 0).any())

    def __array__(self):
        return self.data

    def __iadd__(self, val):
        self.data += np.array(val)
        return self

    def __isub__(self, val):
        self.data -= np.array(val)
        return self

    def __imul__(self, val):
        self.data *= np.array(val)
        return self

    def __imatmul__(self, val):
        self.data @= np.array(val)
        return self

    def __itruediv__(self, val):
        self.data /= np.array(val)
        return self

    def __ipow__(self, val):
        self.data **= np.array(val)
        return self

    def __add__(self, val):
        return Tensor6(self.data + np.array(val))

    def __sub__(self, val):
        return Tensor6(self.data - np.array(val))

    def __mul__(self, val):
        return Tensor6(self.data * np.array(val))

    def __matmul__(self, val):
        return self.data @ np.array(val)

    def __truediv__(self, val):
        return Tensor6(self.data / np.array(val))

    def __pow__(self, val):
        return Tensor6(self.data ** np.array(val))

    def __radd__(self, val):
        return Tensor6(np.array(val) + self.data)

    def __rsub__(self, val):
        return Tensor6(np.array(val) - self.data)

    def __rmul__(self, val):
        return Tensor6(np.array(val) * self.data)

    def __rmatmul__(self, val):
        return np.array(val) @ self.data

    def __rtruediv__(self, val):
        return Tensor6(np.array(val) / self.data)

    def __rpow__(self, val):
        return Tensor6(np.array(val) ** self.data)


# ArrayOfTensor6
exec(array_base(Tensor6))


# ArrayOfArrayOfTensor6
exec(array_base(ArrayOfTensor6))


lib.createTensor6.restype = c.c_void_p
lib.createTensor6.argtypes = []

lib.deleteTensor6.restype = None
lib.deleteTensor6.argtypes = [c.c_void_p]

lib.printTensor6.restype = None
lib.printTensor6.argtypes = [c.c_void_p]

lib.xmlreadTensor6.restype = c.c_long
lib.xmlreadTensor6.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTensor6.restype = c.c_long
lib.xmlsaveTensor6.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.colsTensor6.restype = c.c_long
lib.colsTensor6.argtypes = [c.c_void_p]

lib.rowsTensor6.restype = c.c_long
lib.rowsTensor6.argtypes = [c.c_void_p]

lib.pagesTensor6.restype = c.c_long
lib.pagesTensor6.argtypes = [c.c_void_p]

lib.booksTensor6.restype = c.c_long
lib.booksTensor6.argtypes = [c.c_void_p]

lib.shelvesTensor6.restype = c.c_long
lib.shelvesTensor6.argtypes = [c.c_void_p]

lib.vitrinesTensor6.restype = c.c_long
lib.vitrinesTensor6.argtypes = [c.c_void_p]

lib.resizeTensor6.restype = None
lib.resizeTensor6.argtypes = [c.c_long, c.c_long, c.c_long, c.c_long, c.c_long, c.c_long, c.c_void_p]

lib.getDataTensor6.restype = c.POINTER(c.c_double)
lib.getDataTensor6.argtypes = [c.c_void_p]
