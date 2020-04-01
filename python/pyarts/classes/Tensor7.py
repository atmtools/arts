import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np
from copy import deepcopy as copy

class Tensor7:
    """ ARTS Tensor7 data

    Properties:
        dims:
            Dimensionality of Tensor7 (constexpr Index)

        data:
            The data (numpy-array)

        size:
            Shape of the Tensor7 (tuple)

    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTensor7())
            if data is not None:
                self.data = data

    @property
    def dims(self):
        return 7

    @property
    def flat(self):
        return self.data.flat

    @property
    def data(self):
        """ The data (numpy-array) """
        if not np.any(self.shape == np.zeros(shape=(self.dims))):
            return np.ctypeslib.as_array(lib.getDataTensor7(self.__data__), self.shape)
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
        """ Shape of the Tensor7 (tuple) """
        return lib.librariesTensor7(self.__data__), lib.vitrinesTensor7(self.__data__), lib.shelvesTensor7(self.__data__), lib.booksTensor7(self.__data__), lib.pagesTensor7(self.__data__), lib.rowsTensor7(self.__data__), lib.colsTensor7(self.__data__),

    @shape.setter
    def shape(self, shape):
        if len(shape) == self.dims:
            lib.resizeTensor7(*shape, self.__data__)
        elif len(shape) < self.dims:
            shape = list(shape)
            while len(shape) < self.dims:
                shape.insert(0, 1)
            self.shape = shape
        else:
            raise RuntimeError("Cannot resize to larger dimensions")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTensor7(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTensor7(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Tensor7):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects Tensor7")

    @staticmethod
    def name():
        return "Tensor7"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTensor7(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTensor7(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Tensor7):
            return (self.data == other.data).all()
        else:
            return self.data == other

    def __lt__(self, other):
        if isinstance(other, Tensor7):
            return (self.data < other.data).all()
        else:
            return self.data < other

    def __bool__(self):
        return bool((self.data != 0).any())


# ArrayOfTensor7
exec(array_base(Tensor7))


# ArrayOfArrayOfTensor7
# exec(array_base(ArrayOfTensor7))


lib.createTensor7.restype = c.c_void_p
lib.createTensor7.argtypes = []

lib.deleteTensor7.restype = None
lib.deleteTensor7.argtypes = [c.c_void_p]

lib.printTensor7.restype = None
lib.printTensor7.argtypes = [c.c_void_p]

lib.xmlreadTensor7.restype = c.c_long
lib.xmlreadTensor7.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTensor7.restype = c.c_long
lib.xmlsaveTensor7.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.colsTensor7.restype = c.c_long
lib.colsTensor7.argtypes = [c.c_void_p]

lib.rowsTensor7.restype = c.c_long
lib.rowsTensor7.argtypes = [c.c_void_p]

lib.pagesTensor7.restype = c.c_long
lib.pagesTensor7.argtypes = [c.c_void_p]

lib.booksTensor7.restype = c.c_long
lib.booksTensor7.argtypes = [c.c_void_p]

lib.shelvesTensor7.restype = c.c_long
lib.shelvesTensor7.argtypes = [c.c_void_p]

lib.vitrinesTensor7.restype = c.c_long
lib.vitrinesTensor7.argtypes = [c.c_void_p]

lib.librariesTensor7.restype = c.c_long
lib.librariesTensor7.argtypes = [c.c_void_p]

lib.resizeTensor7.restype = None
lib.resizeTensor7.argtypes = [c.c_long, c.c_long, c.c_long, c.c_long, c.c_long, c.c_long, c.c_long, c.c_void_p]

lib.getDataTensor7.restype = c.POINTER(c.c_double)
lib.getDataTensor7.argtypes = [c.c_void_p]
