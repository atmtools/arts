import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np


class TransmissionMatrix:
    """ ARTS TransmissionMatrix data

    Properties:
        stokes:
            Stokes dimensionality (const Index)

        freqs:
            Number of frequency points (const Index)

        data:
            The data (list of numpy-arrays)
    """
    def __init__(self, stokes=1, freqs=0):
        if isinstance(stokes, c.c_void_p):
            self.__delete__ = False
            self.__data__ = stokes
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTransmissionMatrix())
            self.setTransmissionMatrix(stokes, freqs)

    @staticmethod
    def name():
        return "TransmissionMatrix"

    def setTransmissionMatrix(self, stokes, freqs):
        """ Reinitialize Transmission matrix to new size """
        if stokes > 4 or stokes < 1 or freqs < 0:
            raise ValueError("Bad input")
        else:
            lib.setTransmissionMatrix(self.__data__, int(stokes), int(freqs))

    @property
    def stokes(self):
        """ Stokes dimensionality (const Index) """
        return lib.getStokesDimTransmissionMatrix(self.__data__)

    @property
    def freqs(self):
        """ Number of frequency points (const Index) """
        return lib.getFrequenciesTransmissionMatrix(self.__data__)

    @property
    def data(self):
        """ The data (list of numpy-arrays) """
        x = []
        s = self.stokes
        f = self.freqs
        if s == 1:
            for i in range(f):
                x.append(np.ctypeslib.as_array(lib.getMat1TransmissionMatrix(int(i), self.__data__), (s, s)))
        elif s == 2:
            for i in range(f):
                x.append(np.ctypeslib.as_array(lib.getMat2TransmissionMatrix(int(i), self.__data__), (s, s)))
        elif s == 3:
            for i in range(f):
                x.append(np.ctypeslib.as_array(lib.getMat3TransmissionMatrix(int(i), self.__data__), (s, s)))
        elif s == 4:
            for i in range(f):
                x.append(np.ctypeslib.as_array(lib.getMat4TransmissionMatrix(int(i), self.__data__), (s, s)))
        return x

    @data.setter
    def data(self, val):
        if not isinstance(val, Sized):
            raise TypeError("Only accepts array-like input")
        else:
            f = len(val)

        if f and not hasattr(val[0], "shape"):
            raise TypeError("Only accepts array-like input")
        else:
            s = val[0].shape[0] if f else 1

        if s > 4 or s < 1:
            raise ValueError("Bad input")

        for v in val:
            if not hasattr(v, "shape") or not v.shape == (s, s):
                raise TypeError("Not same types in list")

        self.setTransmissionMatrix(s, f)
        x = self.data

        for i in range(f):
            x[i].flat[:] = val[i].flat[:]

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTransmissionMatrix(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTransmissionMatrix(self.__data__)

    def __repr__(self):
        return "ARTS TransmissionMatrix"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, TransmissionMatrix):
            self.data = other.data
        else:
            raise TypeError("Expects TransmissionMatrix")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTransmissionMatrix(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTransmissionMatrix(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(TransmissionMatrix))


exec(array_base(ArrayOfTransmissionMatrix))


lib.createTransmissionMatrix.restype = c.c_void_p
lib.createTransmissionMatrix.argtypes = []

lib.deleteTransmissionMatrix.restype = None
lib.deleteTransmissionMatrix.argtypes = [c.c_void_p]

lib.printTransmissionMatrix.restype = None
lib.printTransmissionMatrix.argtypes = [c.c_void_p]

lib.xmlreadTransmissionMatrix.restype = c.c_long
lib.xmlreadTransmissionMatrix.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTransmissionMatrix.restype = c.c_long
lib.xmlsaveTransmissionMatrix.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getMat1TransmissionMatrix.restype = c.POINTER(c.c_double)
lib.getMat1TransmissionMatrix.argtypes = [c.c_long, c.c_void_p]

lib.getMat2TransmissionMatrix.restype = c.POINTER(c.c_double)
lib.getMat2TransmissionMatrix.argtypes = [c.c_long, c.c_void_p]

lib.getMat3TransmissionMatrix.restype = c.POINTER(c.c_double)
lib.getMat3TransmissionMatrix.argtypes = [c.c_long, c.c_void_p]

lib.getMat4TransmissionMatrix.restype = c.POINTER(c.c_double)
lib.getMat4TransmissionMatrix.argtypes = [c.c_long, c.c_void_p]

lib.setTransmissionMatrix.restype = None
lib.setTransmissionMatrix.argtypes = [c.c_void_p, c.c_long, c.c_long]

lib.getStokesDimTransmissionMatrix.restype = c.c_long
lib.getStokesDimTransmissionMatrix.argtypes = [c.c_void_p]

lib.getFrequenciesTransmissionMatrix.restype = c.c_long
lib.getFrequenciesTransmissionMatrix.argtypes = [c.c_void_p]
