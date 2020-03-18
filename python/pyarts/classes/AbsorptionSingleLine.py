import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.LineShapeModel import LineShapeModel
from pyarts.classes.Rational import Rational
from pyarts.classes.ZeemanModel import ZeemanModel


class AbsorptionSingleLine:
    """ ARTS Absorption::SingleLine data

    Properties:
        f0:
            Line frequency (Numeric)

        i0:
            Line strength (Numeric)

        e0:
            Line lower state energy (Numeric)

        gl:
            Lower state statistical weight (Numeric)

        gu:
            Upper state statistical weight (Numeric)

        a:
            Einstein coefficient (Numeric)

        zeeman:
            Zeeman model data (ZeemanModel)

        lsm:
            Line shape model (LineShapeModel)

        sizequpp:
            Number of upper state quantum numbers (Index)

        sizeqlow:
            Number of lower state quantum numbers (Index)

        qupp:
           Upper state quantum numbers (list of Rational)

        qlow:
            Lower state quantum numbers (list of Rational)
        """
    def __init__(self, f0=0, i0=0, e0=0, gl=0, gu=0, a=0, zm=ZeemanModel(),
                 lsm=LineShapeModel(), qupp=[], qlow=[]):
        if isinstance(f0, c.c_void_p):
            self.__delete__ = False
            self.__data__ = f0
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createAbsorptionSingleLine())
            self.f0 = f0
            self.i0 = i0
            self.e0 = e0
            self.gl = gl
            self.gu = gu
            self.a = a
            self.zeeman = zm
            self.lsm = lsm
            self.qupp = qupp
            self.qlow = qlow

    @property
    def f0(self):
        """ Line frequency (Numeric) """
        return lib.getAbsorptionSingleLineF0(self.__data__)

    @f0.setter
    def f0(self, x):
        lib.setAbsorptionSingleLineF0(self.__data__, c.c_double(x))

    @property
    def i0(self):
        """ Line strength (Numeric) """
        return lib.getAbsorptionSingleLineI0(self.__data__)

    @i0.setter
    def i0(self, x):
        lib.setAbsorptionSingleLineI0(self.__data__, c.c_double(x))

    @property
    def e0(self):
        """ Line lower state energy (Numeric) """
        return lib.getAbsorptionSingleLineE0(self.__data__)

    @e0.setter
    def e0(self, x):
        lib.setAbsorptionSingleLineE0(self.__data__, c.c_double(x))

    @property
    def gl(self):
        """ Lower state statistical weight (Numeric) """
        return lib.getAbsorptionSingleLineGL(self.__data__)

    @gl.setter
    def gl(self, x):
        lib.setAbsorptionSingleLineGL(self.__data__, c.c_double(x))

    @property
    def gu(self):
        """ Upper state statistical weight (Numeric) """
        return lib.getAbsorptionSingleLineGU(self.__data__)

    @gu.setter
    def gu(self, x):
        lib.setAbsorptionSingleLineGU(self.__data__, c.c_double(x))

    @property
    def a(self):
        return lib.getAbsorptionSingleLineA(self.__data__)

    @a.setter
    def a(self, x):
        """ Einstein coefficient (Numeric) """
        lib.setAbsorptionSingleLineA(self.__data__, c.c_double(x))

    @property
    def zeeman(self):
        """ Zeeman model data (ZeemanModel) """
        return ZeemanModel(c.c_void_p(lib.getAbsorptionSingleLineZeemanModel(self.__data__)))

    @zeeman.setter
    def zeeman(self, x):
        self.zeeman.set(x)

    @property
    def lsm(self):
        """ Line shape model (LineShapeModel) """
        return LineShapeModel(c.c_void_p(lib.getAbsorptionSingleLineLineShapeModel(self.__data__)))

    @lsm.setter
    def lsm(self, x):
        self.lsm.set(x)

    @property
    def sizequpp(self):
        """ Number of upper state quantum numbers (Index)

        Note that the upper and lower state quantum numbers must be equally
        many for the class to be considered OK
        """
        return lib.sizeAbsorptionSingleLineUpperQuantas(self.__data__)

    @sizequpp.setter
    def sizequpp(self, n):
        n = int(n)
        lib.resizeAbsorptionSingleLineUpperQuantas(c.c_long(n), self.__data__)

    @property
    def sizeqlow(self):
        """ Number of lower state quantum numbers (Index)

        Note that the upper and lower state quantum numbers must be equally
        many for the class to be considered OK
        """
        return lib.sizeAbsorptionSingleLineLowerQuantas(self.__data__)

    @sizeqlow.setter
    def sizeqlow(self, n):
        n = int(n)
        lib.resizeAbsorptionSingleLineLowerQuantas(c.c_long(n), self.__data__)

    @property
    def qupp(self):
        """ Upper state quantum numbers (list of Rational) """
        x = []
        n = self.sizequpp
        for i in range(n):
            x.append(Rational(c.c_void_p(lib.getAbsorptionSingleLineUpperQuanta(i, self.__data__))))
        return x

    @qupp.setter
    def qupp(self, val):
        if isinstance(val, Sized):
            self.sizequpp = len(val)
            n = self.sizequpp
            x = self.qupp
            for i in range(n):
                x[i].set(val[i])
        else:
            raise TypeError("Only accepts array-like input")

    @property
    def qlow(self):
        """ Lower state quantum numbers (list of Rational) """
        x = []
        n = self.sizeqlow
        for i in range(n):
            x.append(Rational(c.c_void_p(lib.getAbsorptionSingleLineLowerQuanta(i, self.__data__))))
        return x

    @qlow.setter
    def qlow(self, val):
        if isinstance(val, Sized):
            self.sizeqlow = len(val)
            n = self.sizeqlow
            x = self.qlow
            for i in range(n):
                x[i].set(val[i])
        else:
            raise TypeError("Only accepts array-like input")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAbsorptionSingleLine(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionSingleLine(self.__data__)

    def __repr__(self):
        return "ARTS Absorption::SingleLine"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, AbsorptionSingleLine):
            self.f0 = other.f0
            self.i0 = other.i0
            self.e0 = other.e0
            self.gl = other.gl
            self.gu = other.gu
            self.a = other.a
            self.zeeman = other.zeeman
            self.lsm = other.lsm
            self.qupp = other.qupp
            self.qlow = other.qlow
        else:
            raise TypeError("Expects AbsorptionSingleLine")


lib.createAbsorptionSingleLine.restype = c.c_void_p
lib.createAbsorptionSingleLine.argtypes = []

lib.deleteAbsorptionSingleLine.restype = None
lib.deleteAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.printAbsorptionSingleLine.restype = None
lib.printAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getAbsorptionSingleLineF0.restype = c.c_double
lib.getAbsorptionSingleLineF0.argtypes = [c.c_void_p]

lib.getAbsorptionSingleLineI0.restype = c.c_double
lib.getAbsorptionSingleLineI0.argtypes = [c.c_void_p]

lib.getAbsorptionSingleLineE0.restype = c.c_double
lib.getAbsorptionSingleLineE0.argtypes = [c.c_void_p]

lib.getAbsorptionSingleLineGL.restype = c.c_double
lib.getAbsorptionSingleLineGL.argtypes = [c.c_void_p]

lib.getAbsorptionSingleLineGU.restype = c.c_double
lib.getAbsorptionSingleLineGU.argtypes = [c.c_void_p]

lib.getAbsorptionSingleLineA.restype = c.c_double
lib.getAbsorptionSingleLineA.argtypes = [c.c_void_p]

lib.setAbsorptionSingleLineF0.restype = None
lib.setAbsorptionSingleLineF0.argtypes = [c.c_void_p, c.c_double]

lib.setAbsorptionSingleLineI0.restype = None
lib.setAbsorptionSingleLineI0.argtypes = [c.c_void_p, c.c_double]

lib.setAbsorptionSingleLineE0.restype = None
lib.setAbsorptionSingleLineE0.argtypes = [c.c_void_p, c.c_double]

lib.setAbsorptionSingleLineGL.restype = None
lib.setAbsorptionSingleLineGL.argtypes = [c.c_void_p, c.c_double]

lib.setAbsorptionSingleLineGU.restype = None
lib.setAbsorptionSingleLineGU.argtypes = [c.c_void_p, c.c_double]

lib.setAbsorptionSingleLineA.restype = None
lib.setAbsorptionSingleLineA.argtypes = [c.c_void_p, c.c_double]

lib.getAbsorptionSingleLineZeemanModel.restype = c.c_void_p
lib.getAbsorptionSingleLineZeemanModel.argtypes = [c.c_void_p]

lib.getAbsorptionSingleLineLineShapeModel.restype = c.c_void_p
lib.getAbsorptionSingleLineLineShapeModel.argtypes = [c.c_void_p]

lib.sizeAbsorptionSingleLineLowerQuantas.restype = c.c_long
lib.sizeAbsorptionSingleLineLowerQuantas.argtypes = [c.c_void_p]

lib.sizeAbsorptionSingleLineUpperQuantas.restype = c.c_long
lib.sizeAbsorptionSingleLineUpperQuantas.argtypes = [c.c_void_p]

lib.resizeAbsorptionSingleLineLowerQuantas.restype = None
lib.resizeAbsorptionSingleLineLowerQuantas.argtypes = [c.c_long, c.c_void_p]

lib.resizeAbsorptionSingleLineUpperQuantas.restype = None
lib.resizeAbsorptionSingleLineUpperQuantas.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionSingleLineLowerQuanta.restype = c.c_void_p
lib.getAbsorptionSingleLineLowerQuanta.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionSingleLineUpperQuanta.restype = c.c_void_p
lib.getAbsorptionSingleLineUpperQuanta.argtypes = [c.c_long, c.c_void_p]
