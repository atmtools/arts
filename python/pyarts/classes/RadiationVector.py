import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np


class RadiationVector:
    """ ARTS RadiationVector data

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
            self.__data__ = c.c_void_p(lib.createRadiationVector())
            self.setRadiationVector(stokes, freqs)

    @staticmethod
    def name():
        return "RadiationVector"

    def setRadiationVector(self, stokes, freqs):
        """ Reinitialize Transmission matrix to new size """
        if stokes > 4 or stokes < 1 or freqs < 0:
            raise ValueError("Bad input")
        else:
            lib.setRadiationVector(self.__data__, int(stokes), int(freqs))

    @property
    def stokes(self):
        """ Stokes dimensionality (const Index) """
        return lib.getStokesDimRadiationVector(self.__data__)

    @property
    def freqs(self):
        """ Number of frequency points (const Index) """
        return lib.getFrequenciesRadiationVector(self.__data__)

    @property
    def data(self):
        """ The data (list of numpy-arrays) """
        x = []
        for i in range(self.freqs):
            x.append(self[i])
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
            s = val[0].shape[0] if f else self.stokes

        if s > 4 or s < 1:
            raise ValueError("Bad input")

        self.setRadiationVector(s, f)
        for i in range(f):
            self[i] = val[i]

    def __getitem__(self, i):
        if self.stokes == 1 and i < self.freqs:
            return np.ctypeslib.as_array(lib.getVec1RadiationVector(int(i) % self.freqs, self.__data__), (1,))
        elif self.stokes == 2:
            return np.ctypeslib.as_array(lib.getVec2RadiationVector(int(i) % self.freqs, self.__data__), (2,))
        elif self.stokes == 3:
            return np.ctypeslib.as_array(lib.getVec3RadiationVector(int(i) % self.freqs, self.__data__), (3,))
        elif self.stokes == 4:
            return np.ctypeslib.as_array(lib.getVec4RadiationVector(int(i) % self.freqs, self.__data__), (4,))

    def __setitem__(self, i, x):
        self[i].flat[:] = x.flat[:]

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printRadiationVector(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteRadiationVector(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def __len__(self):
        return self.freqs

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, RadiationVector):
            lib.setRadiationVector(self.__data__, self.stokes, self.freqs)
            self.data = other.data
        else:
            raise TypeError("Expects RadiationVector")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadRadiationVector(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveRadiationVector(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Sized) and len(other) == len(self):
            for i in range(len(self)):
                if (self[i] != other[i]).any():
                    return False
            return True
        else:
            return False


exec(array_base(RadiationVector))


exec(array_base(ArrayOfRadiationVector))


lib.createRadiationVector.restype = c.c_void_p
lib.createRadiationVector.argtypes = []

lib.deleteRadiationVector.restype = None
lib.deleteRadiationVector.argtypes = [c.c_void_p]

lib.printRadiationVector.restype = None
lib.printRadiationVector.argtypes = [c.c_void_p]

lib.xmlreadRadiationVector.restype = c.c_long
lib.xmlreadRadiationVector.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveRadiationVector.restype = c.c_long
lib.xmlsaveRadiationVector.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getVec1RadiationVector.restype = c.POINTER(c.c_double)
lib.getVec1RadiationVector.argtypes = [c.c_long, c.c_void_p]

lib.getVec2RadiationVector.restype = c.POINTER(c.c_double)
lib.getVec2RadiationVector.argtypes = [c.c_long, c.c_void_p]

lib.getVec3RadiationVector.restype = c.POINTER(c.c_double)
lib.getVec3RadiationVector.argtypes = [c.c_long, c.c_void_p]

lib.getVec4RadiationVector.restype = c.POINTER(c.c_double)
lib.getVec4RadiationVector.argtypes = [c.c_long, c.c_void_p]

lib.setRadiationVector.restype = None
lib.setRadiationVector.argtypes = [c.c_void_p, c.c_long, c.c_long]

lib.getStokesDimRadiationVector.restype = c.c_long
lib.getStokesDimRadiationVector.argtypes = [c.c_void_p]

lib.getFrequenciesRadiationVector.restype = c.c_long
lib.getFrequenciesRadiationVector.argtypes = [c.c_void_p]
