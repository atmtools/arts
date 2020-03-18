import ctypes as c
from pyarts.workspace.api import arts_api as lib


class ZeemanModel:
    """ ARTS Zeeman::Model data

    Properties:
        gu:
            Upper level Zeeman splitting (Numeric)

        gl:
            Lower level Zeeman splitting (Numeric)
        """

    def __init__(self, gu=float('nan'), gl=float('nan')):
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
        return lib.getZeemanModelGU(self.__data__)

    @gu.setter
    def gu(self, x):
        x = float(x)
        lib.setZeemanModelGU(self.__data__, x)

    @property
    def gl(self):
        """ Lower level Zeeman splitting (Numeric) """
        return lib.getZeemanModelGL(self.__data__)

    @gl.setter
    def gl(self, x):
        x = float(x)
        lib.setZeemanModelGL(self.__data__, x)

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


lib.createZeemanModel.restype = c.c_void_p
lib.createZeemanModel.argtypes = []

lib.deleteZeemanModel.restype = None
lib.deleteZeemanModel.argtypes = [c.c_void_p]

lib.printZeemanModel.restype = None
lib.printZeemanModel.argtypes = [c.c_void_p]

lib.getZeemanModelGU.restype = c.c_double
lib.getZeemanModelGU.argtypes = [c.c_void_p]

lib.getZeemanModelGL.restype = c.c_double
lib.getZeemanModelGL.argtypes = [c.c_void_p]

lib.setZeemanModelGU.restype = None
lib.setZeemanModelGU.argtypes = [c.c_void_p, c.c_double]

lib.setZeemanModelGL.restype = None
lib.setZeemanModelGL.argtypes = [c.c_void_p, c.c_double]
