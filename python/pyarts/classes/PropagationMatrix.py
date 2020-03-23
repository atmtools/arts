import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Tensor4 import Tensor4
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np


class PropagationMatrix:
    """ ARTS PropagationMatrix data

    Properties:
        data:
            The data (numpy-array)

        stokes:
            Stokes dimension (const Index)

        freqs:
            Frequency dimension (const Index)

        zeniths:
            Zenith dimension (const Index)

        azimuths:
            Azimuth dimension (const Index)

    """
    def __init__(self, data=0.0, freqs=0, stokes=1, za=1, aa=1):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = freqs
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createPropagationMatrix())
            self.setData(freqs, stokes, za, aa, data)

    @property
    def data(self):
        """ The data (numpy-array) """
        return Tensor4(c.c_void_p(lib.getDataPropagationMatrix(self.__data__))).data

    @data.setter
    def data(self, value):
        if value.shape == self.data.shape:
            self.data.flat[:] = value[:]
        else:
            raise TypeError("Expects same-size data")

    @property
    def stokes(self):
        """ Stokes dimension (const Index) """
        return lib.stokesPropagationMatrix(self.__data__)

    @property
    def frequencies(self):
        """ Frequency dimension (const Index) """
        return lib.frequenciesPropagationMatrix(self.__data__)

    @property
    def zeniths(self):
        """ Zenith dimension (const Index) """
        return lib.zenithsPropagationMatrix(self.__data__)

    @property
    def azimuths(self):
        """ Azimuth dimension (const Index) """
        return lib.azimuthsPropagationMatrix(self.__data__)

    def setData(self, freqs, stokes, za, aa, data):
        """ Sets the data by reinitialization to new size """
        if lib.setPropagationMatrix(self.__data__, int(freqs), int(stokes), int(za), int(aa), float(data)):
            raise ValueError("Bad input")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printPropagationMatrix(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deletePropagationMatrix(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, PropagationMatrix):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects PropagationMatrix")

    @staticmethod
    def name():
        return "PropagationMatrix"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadPropagationMatrix(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsavePropagationMatrix(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


# ArrayOfPropagationMatrix
exec(array_base(PropagationMatrix))


# ArrayOfArrayOfPropagationMatrix
exec(array_base(ArrayOfPropagationMatrix))


lib.createPropagationMatrix.restype = c.c_void_p
lib.createPropagationMatrix.argtypes = []

lib.deletePropagationMatrix.restype = None
lib.deletePropagationMatrix.argtypes = [c.c_void_p]

lib.printPropagationMatrix.restype = None
lib.printPropagationMatrix.argtypes = [c.c_void_p]

lib.xmlreadPropagationMatrix.restype = c.c_long
lib.xmlreadPropagationMatrix.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsavePropagationMatrix.restype = c.c_long
lib.xmlsavePropagationMatrix.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.stokesPropagationMatrix.restype = c.c_long
lib.stokesPropagationMatrix.argtypes = [c.c_void_p]

lib.frequenciesPropagationMatrix.restype = c.c_long
lib.frequenciesPropagationMatrix.argtypes = [c.c_void_p]

lib.zenithsPropagationMatrix.restype = c.c_long
lib.zenithsPropagationMatrix.argtypes = [c.c_void_p]

lib.azimuthsPropagationMatrix.restype = c.c_long
lib.azimuthsPropagationMatrix.argtypes = [c.c_void_p]

lib.getDataPropagationMatrix.restype = c.c_void_p
lib.getDataPropagationMatrix.argtypes = [c.c_void_p]

lib.setPropagationMatrix.restype = c.c_long
lib.setPropagationMatrix.argtypes = [c.c_void_p, c.c_long, c.c_long, c.c_long, c.c_long, c.c_double]
