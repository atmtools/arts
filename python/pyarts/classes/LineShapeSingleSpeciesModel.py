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
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelG0(self.__data__).contents)

    @g0.setter
    def g0(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelG0(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelG0(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def d0(self):
        """ Speed independent shifting model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelD0(self.__data__).contents)

    @d0.setter
    def d0(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelD0(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelD0(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def g2(self):
        """ Speed dependent broadening model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelG2(self.__data__).contents)

    @g2.setter
    def g2(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelG2(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelG2(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def d2(self):
        """ Speed dependent shifting model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelD2(self.__data__).contents)

    @d2.setter
    def d2(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelD2(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelD2(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def fvc(self):
        """ Frequency velocity changing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelFVC(self.__data__).contents)

    @fvc.setter
    def fvc(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelFVC(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelFVC(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def eta(self):
        """ Correlation model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelETA(self.__data__).contents)

    @eta.setter
    def eta(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelETA(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelETA(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def y(self):
        """ 1st order line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelY(self.__data__).contents)

    @y.setter
    def y(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelY(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelY(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def g(self):
        """ 2nd order strength line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelG(self.__data__).contents)

    @g.setter
    def g(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelG(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelG(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

    @property
    def dv(self):
        """ 2nd order shift line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(lib.getLineShapeSingleSpeciesModelDV(self.__data__).contents)

    @dv.setter
    def dv(self, x):
        if isinstance(x, LineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelDV(self.__data__, c.pointer(x.data))
        elif isinstance(x, InternalLineShapeModelParameters):
            lib.setLineShapeSingleSpeciesModelDV(self.__data__, c.pointer(x))
        else:
            raise TypeError("Only accepts LineShapeModelParameters input")

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

lib.getLineShapeSingleSpeciesModelG0.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelG0.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelD0.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelD0.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelG2.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelG2.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelD2.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelD2.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelFVC.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelFVC.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelETA.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelETA.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelY.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelY.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelG.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelG.argtypes = [c.c_void_p]

lib.getLineShapeSingleSpeciesModelDV.restype = c.POINTER(InternalLineShapeModelParameters)
lib.getLineShapeSingleSpeciesModelDV.argtypes = [c.c_void_p]

lib.setLineShapeSingleSpeciesModelG0.restype = None
lib.setLineShapeSingleSpeciesModelG0.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelD0.restype = None
lib.setLineShapeSingleSpeciesModelD0.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelG2.restype = None
lib.setLineShapeSingleSpeciesModelG2.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelD2.restype = None
lib.setLineShapeSingleSpeciesModelD2.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelFVC.restype = None
lib.setLineShapeSingleSpeciesModelFVC.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelETA.restype = None
lib.setLineShapeSingleSpeciesModelETA.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelY.restype = None
lib.setLineShapeSingleSpeciesModelY.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelG.restype = None
lib.setLineShapeSingleSpeciesModelG.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]

lib.setLineShapeSingleSpeciesModelDV.restype = None
lib.setLineShapeSingleSpeciesModelDV.argtypes = [c.c_void_p, c.POINTER(InternalLineShapeModelParameters)]
