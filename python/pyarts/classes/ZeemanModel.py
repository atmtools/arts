import ctypes as c
from pyarts.workspace.api import arts_api as lib

from math import isnan, nan

class ZeemanModel:
    """ ARTS Zeeman::Model data

    Properties:
        gu:
            Upper level Zeeman splitting (Numeric)

        gl:
            Lower level Zeeman splitting (Numeric)
        """

    def __init__(self, gu=nan, gl=nan):
        if isinstance(gu, c.c_void_p):
            self.__delete__ = False
            self.__data__ = gu
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createZeemanModel())
            self.gu = gu
            self.gl = gl

    @property
    def gu(self):
        """ Upper level Zeeman splitting (Numeric) """
        return lib.getguZeemanModel(self.__data__)

    @gu.setter
    def gu(self, x):
        lib.setguZeemanModel(self.__data__, float(x))

    @property
    def gl(self):
        """ Lower level Zeeman splitting (Numeric) """
        return lib.getglZeemanModel(self.__data__)

    @gl.setter
    def gl(self, x):
        lib.setglZeemanModel(self.__data__, float(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printZeemanModel(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteZeemanModel(self.__data__)

    def __repr__(self):
        return "ARTS Zeeman::Model {} {}".format(self.gu, self.gl)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, ZeemanModel):
            self.gl = other.gl
            self.gu = other.gu
        else:
            raise TypeError("Expects ZeemanModel")

    def __eq__(self, other):
        if isinstance(other, ZeemanModel) and \
            (self.gl == other.gl or (isnan(self.gl) and isnan(other.gl)))  and \
            (self.gu == other.gu or (isnan(self.gu) and isnan(other.gu))):
            return True
        else:
            return False


lib.createZeemanModel.restype = c.c_void_p
lib.createZeemanModel.argtypes = []

lib.deleteZeemanModel.restype = None
lib.deleteZeemanModel.argtypes = [c.c_void_p]

lib.printZeemanModel.restype = None
lib.printZeemanModel.argtypes = [c.c_void_p]

lib.getguZeemanModel.restype = c.c_double
lib.getguZeemanModel.argtypes = [c.c_void_p]

lib.getglZeemanModel.restype = c.c_double
lib.getglZeemanModel.argtypes = [c.c_void_p]

lib.setguZeemanModel.restype = None
lib.setguZeemanModel.argtypes = [c.c_void_p, c.c_double]

lib.setglZeemanModel.restype = None
lib.setglZeemanModel.argtypes = [c.c_void_p, c.c_double]
