import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Tensor4 import Tensor4
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np


class StokesVector:
    """ ARTS StokesVector data

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
            self.__data__ = c.c_void_p(lib.createStokesVector())
            self.setData(freqs, stokes, za, aa, data)

    @property
    def data(self):
        """ The data (numpy-array) """
        return Tensor4(c.c_void_p(lib.getDataStokesVector(self.__data__))).data

    @data.setter
    def data(self, value):
        if value.shape == self.data.shape:
            self.data.flat[:] = value[:]
        else:
            raise TypeError("Expects same-size data")

    @property
    def stokes(self):
        """ Stokes dimension (const Index) """
        return lib.stokesStokesVector(self.__data__)

    @property
    def frequencies(self):
        """ Frequency dimension (const Index) """
        return lib.frequenciesStokesVector(self.__data__)

    @property
    def zeniths(self):
        """ Zenith dimension (const Index) """
        return lib.zenithsStokesVector(self.__data__)

    @property
    def azimuths(self):
        """ Azimuth dimension (const Index) """
        return lib.azimuthsStokesVector(self.__data__)

    def setData(self, freqs, stokes, za, aa, data):
        """ Sets the data by reinitialization to new size """
        if lib.setStokesVector(self.__data__, int(freqs), int(stokes), int(za), int(aa), float(data)):
            raise ValueError("Bad input")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printStokesVector(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteStokesVector(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, StokesVector):
            self.shape = other.shape
            self.data = other.data
        else:
            raise TypeError("Expects StokesVector")

    @staticmethod
    def name():
        return "StokesVector"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadStokesVector(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveStokesVector(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


# ArrayOfStokesVector
exec(array_base(StokesVector))


# ArrayOfArrayOfStokesVector
exec(array_base(ArrayOfStokesVector))


lib.createStokesVector.restype = c.c_void_p
lib.createStokesVector.argtypes = []

lib.deleteStokesVector.restype = None
lib.deleteStokesVector.argtypes = [c.c_void_p]

lib.printStokesVector.restype = None
lib.printStokesVector.argtypes = [c.c_void_p]

lib.xmlreadStokesVector.restype = c.c_long
lib.xmlreadStokesVector.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveStokesVector.restype = c.c_long
lib.xmlsaveStokesVector.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.stokesStokesVector.restype = c.c_long
lib.stokesStokesVector.argtypes = [c.c_void_p]

lib.frequenciesStokesVector.restype = c.c_long
lib.frequenciesStokesVector.argtypes = [c.c_void_p]

lib.zenithsStokesVector.restype = c.c_long
lib.zenithsStokesVector.argtypes = [c.c_void_p]

lib.azimuthsStokesVector.restype = c.c_long
lib.azimuthsStokesVector.argtypes = [c.c_void_p]

lib.getDataStokesVector.restype = c.c_void_p
lib.getDataStokesVector.argtypes = [c.c_void_p]

lib.setStokesVector.restype = c.c_long
lib.setStokesVector.argtypes = [c.c_void_p, c.c_long, c.c_long, c.c_long, c.c_long, c.c_double]
