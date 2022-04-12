import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.LineShapeModelParameters import LineShapeModelParameters


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
        return LineShapeModelParameters(c.c_void_p(lib.getG0LineShapeSingleSpeciesModel(self.__data__)))

    @g0.setter
    def g0(self, x):
        self.g0.set(x)

    @property
    def d0(self):
        """ Speed independent shifting model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getD0LineShapeSingleSpeciesModel(self.__data__)))

    @d0.setter
    def d0(self, x):
        self.d0.set(x)

    @property
    def g2(self):
        """ Speed dependent broadening model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getG2LineShapeSingleSpeciesModel(self.__data__)))

    @g2.setter
    def g2(self, x):
        self.g2.set(x)

    @property
    def d2(self):
        """ Speed dependent shifting model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getD2LineShapeSingleSpeciesModel(self.__data__)))

    @d2.setter
    def d2(self, x):
        self.d2.set(x)

    @property
    def fvc(self):
        """ Frequency velocity changing model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getFVCLineShapeSingleSpeciesModel(self.__data__)))

    @fvc.setter
    def fvc(self, x):
        self.fvc.set(x)

    @property
    def eta(self):
        """ Correlation model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getETALineShapeSingleSpeciesModel(self.__data__)))

    @eta.setter
    def eta(self, x):
        self.eta.set(x)

    @property
    def y(self):
        """ 1st order line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getYLineShapeSingleSpeciesModel(self.__data__)))

    @y.setter
    def y(self, x):
        self.y.set(x)

    @property
    def g(self):
        """ 2nd order strength line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getGLineShapeSingleSpeciesModel(self.__data__)))

    @g.setter
    def g(self, x):
        self.g.set(x)

    @property
    def dv(self):
        """ 2nd order shift line mixing model (LineShapeModelParameters) """
        return LineShapeModelParameters(c.c_void_p(lib.getDVLineShapeSingleSpeciesModel(self.__data__)))

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
        out = ""
        if self.g0.type != "None":
            if len(out): out += ' '
            out += f"G0 {self.g0}"
        if self.d0.type != "None":
            if len(out): out += ' '
            out += f"D0 {self.d0}"
        if self.g2.type != "None":
            if len(out): out += ' '
            out += f"G2 {self.g2}"
        if self.d2.type != "None":
            if len(out): out += ' '
            out += f"D2 {self.d2}"
        if self.fvc.type != "None":
            if len(out): out += ' '
            out += f"FVC {self.fvc}"
        if self.eta.type != "None":
            if len(out): out += ' '
            out += f"ETA {self.eta}"
        if self.y.type != "None":
            if len(out): out += ' '
            out += f"Y {self.y}"
        if self.g.type != "None":
            if len(out): out += ' '
            out += f"G {self.g}"
        if self.dv.type != "None":
            if len(out): out += ' '
            out += f"DV {self.dv}"
        return out

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

    def __eq__(self, other):
        if isinstance(other, LineShapeSingleSpeciesModel) and \
                self.g0 == other.g0 and \
                self.d0 == other.d0 and \
                self.g2 == other.g2 and \
                self.d2 == other.d2 and \
                self.fvc == other.fvc and \
                self.eta == other.eta and \
                self.y == other.y and \
                self.g == other.g and \
                self.dv == other.dv:
            return True
        else:
            return False


lib.createLineShapeSingleSpeciesModel.restype = c.c_void_p
lib.createLineShapeSingleSpeciesModel.argtypes = []

lib.deleteLineShapeSingleSpeciesModel.restype = None
lib.deleteLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.printLineShapeSingleSpeciesModel.restype = None
lib.printLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getG0LineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getG0LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getD0LineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getD0LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getG2LineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getG2LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getD2LineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getD2LineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getFVCLineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getFVCLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getETALineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getETALineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getYLineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getYLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getGLineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getGLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]

lib.getDVLineShapeSingleSpeciesModel.restype = c.c_void_p
lib.getDVLineShapeSingleSpeciesModel.argtypes = [c.c_void_p]
