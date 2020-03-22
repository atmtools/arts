import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.LineShapeModelParameters import LineShapeModelParameters
from pyarts.classes.LineShapeModelParameters import InternalLineShapeModelParameters


class LineShapeSingleSpeciesModel:
    """ ARTS LineShape::SingleSpeciesModel data

    Properties:
        g0:
            Speed independent broadening model (LineShapeModelParameters)

        d0:
            Speed independent shifting model (LineShapeModelParameters)

        g2:
            Speed dependent broadening model (LineShapeModelParameters)

        d2:
            Speed dependent shifting model (LineShapeModelParameters)

        fvc:
            Frequency velocity changing model (LineShapeModelParameters)

        eta:
            Correlation model (LineShapeModelParameters)

        y:
            1st order line mixing model (LineShapeModelParameters)

        g:
            2nd order strength line mixing model (LineShapeModelParameters)

        dv:
            2nd order shift line mixing model (LineShapeModelParameters)
        """
    def __init__(self, g0=LineShapeModelParameters(), d0=LineShapeModelParameters(),
                       g2=LineShapeModelParameters(), d2=LineShapeModelParameters(),
                       fvc=LineShapeModelParameters(), eta=LineShapeModelParameters(),
                       y=LineShapeModelParameters(), g=LineShapeModelParameters(),
                       dv=LineShapeModelParameters()):
        if isinstance(g0, c.c_void_p):
            self.__delete__ = False
            self.__data__ = g0
        else:
            self.__delete__ = True
            self.__data__ = c.cast(lib.createLineShapeSingleSpeciesModel(), c.c_void_p)
            self.g0 = g0
            self.d0 = d0
            self.g2 = g2
            self.d2 = d2
            self.fvc = fvc
            self.eta = eta
            self.y = y
            self.g = g
            self.dv = dv

    @property
    def g0(self):
        """ Speed independent broadening model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getG0LineShapeSingleSpeciesModel(self.__data__).contents)

    @g0.setter
    def g0(self, x):
        self.g0.set(x)

    @property
    def d0(self):
        """ Speed independent shifting model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getD0LineShapeSingleSpeciesModel(self.__data__).contents)

    @d0.setter
    def d0(self, x):
        self.d0.set(x)

    @property
    def g2(self):
        """ Speed dependent broadening model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getG2LineShapeSingleSpeciesModel(self.__data__).contents)

    @g2.setter
    def g2(self, x):
        self.g2.set(x)

    @property
    def d2(self):
        """ Speed dependent shifting model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getD2LineShapeSingleSpeciesModel(self.__data__).contents)

    @d2.setter
    def d2(self, x):
        self.d2.set(x)

    @property
    def fvc(self):
        """ Frequency velocity changing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getFVCLineShapeSingleSpeciesModel(self.__data__).contents)

    @fvc.setter
    def fvc(self, x):
        self.fvc.set(x)

    @property
    def eta(self):
        """ Correlation model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getETALineShapeSingleSpeciesModel(self.__data__).contents)

    @eta.setter
    def eta(self, x):
        self.eta.set(x)

    @property
    def y(self):
        """ 1st order line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getYLineShapeSingleSpeciesModel(self.__data__).contents)

    @y.setter
    def y(self, x):
        self.y.set(x)

    @property
    def g(self):
        """ 2nd order strength line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getGLineShapeSingleSpeciesModel(self.__data__).contents)

    @g.setter
    def g(self, x):
        self.g.set(x)

    @property
    def dv(self):
        """ 2nd order shift line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getDVLineShapeSingleSpeciesModel(self.__data__).contents)

    @dv.setter
    def dv(self, x):
        self.dv.set(x)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printLineShapeSingleSpeciesModel(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteLineShapeSingleSpeciesModel(self.__data__)

    def __repr__(self):
        return "ARTS LineShape::SingleSpeciesModel"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, LineShapeSingleSpeciesModel):
            self.g0 = other.g0
            self.d0 = other.d0
            self.g2 = other.g2
            self.d2 = other.d2
            self.fvc = other.fvc
            self.eta = other.eta
            self.y = other.y
            self.g = other.g
            self.dv = other.dv
        else:
            raise TypeError("Expects LineShapeSingleSpeciesModel")


lib.createLineShapeSingleSpeciesModel.restype = c.c_void_p
lib.createLineShapeSingleSpeciesModel.argtypes = []

lib.deleteLineShapeSingleSpeciesModel.restype = None
lib.deleteLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.printLineShapeSingleSpeciesModel.restype = None
lib.printLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getG0LineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getG0LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getD0LineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getD0LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getG2LineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getG2LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getD2LineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getD2LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getFVCLineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getFVCLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getETALineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getETALineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getYLineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getYLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getGLineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getGLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getDVLineShapeSingleSpeciesModel.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getDVLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]
