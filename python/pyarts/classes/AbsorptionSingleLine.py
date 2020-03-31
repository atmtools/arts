import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.LineShapeModel import LineShapeModel
from pyarts.classes.Rational import Rational
from pyarts.classes.ZeemanModel import ZeemanModel

from math import isnan


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
        return lib.getF0AbsorptionSingleLine(self.__data__)

    @f0.setter
    def f0(self, x):
        x = float(x)
        lib.setF0AbsorptionSingleLine(self.__data__, x)

    @property
    def i0(self):
        """ Line strength (Numeric) """
        return lib.getI0AbsorptionSingleLine(self.__data__)

    @i0.setter
    def i0(self, x):
        x = float(x)
        lib.setI0AbsorptionSingleLine(self.__data__, x)

    @property
    def e0(self):
        """ Line lower state energy (Numeric) """
        return lib.getE0AbsorptionSingleLine(self.__data__)

    @e0.setter
    def e0(self, x):
        x = float(x)
        lib.setE0AbsorptionSingleLine(self.__data__, x)

    @property
    def gl(self):
        """ Lower state statistical weight (Numeric) """
        return lib.getg_lowAbsorptionSingleLine(self.__data__)

    @gl.setter
    def gl(self, x):
        x = float(x)
        lib.setg_lowAbsorptionSingleLine(self.__data__, x)

    @property
    def gu(self):
        """ Upper state statistical weight (Numeric) """
        return lib.getg_uppAbsorptionSingleLine(self.__data__)

    @gu.setter
    def gu(self, x):
        x = float(x)
        lib.setg_uppAbsorptionSingleLine(self.__data__, x)

    @property
    def a(self):
        """ Einstein coefficient (Numeric) """
        return lib.getAAbsorptionSingleLine(self.__data__)

    @a.setter
    def a(self, x):
        x = float(x)
        lib.setAAbsorptionSingleLine(self.__data__, x)

    @property
    def zeeman(self):
        """ Zeeman model data (ZeemanModel) """
        return ZeemanModel(c.c_void_p(lib.getZeemanAbsorptionSingleLine(self.__data__)))

    @zeeman.setter
    def zeeman(self, x):
        self.zeeman.set(x)

    @property
    def lsm(self):
        """ Line shape model (LineShapeModel) """
        return LineShapeModel(c.c_void_p(lib.getLineShapeAbsorptionSingleLine(self.__data__)))

    @lsm.setter
    def lsm(self, x):
        self.lsm.set(x)

    @property
    def sizequpp(self):
        """ Number of upper state quantum numbers (Index)

        Note that the upper and lower state quantum numbers must be equally
        many for the class to be considered OK
        """
        return lib.sizeUpperQuantumNumbersAbsorptionSingleLine(self.__data__)

    @sizequpp.setter
    def sizequpp(self, n):
        n = int(n)
        lib.resizeUpperQuantumNumbersAbsorptionSingleLine(c.c_long(n), self.__data__)

    @property
    def sizeqlow(self):
        """ Number of lower state quantum numbers (Index)

        Note that the upper and lower state quantum numbers must be equally
        many for the class to be considered OK
        """
        return lib.sizeLowerQuantumNumbersAbsorptionSingleLine(self.__data__)

    @sizeqlow.setter
    def sizeqlow(self, n):
        n = int(n)
        lib.resizeLowerQuantumNumbersAbsorptionSingleLine(c.c_long(n), self.__data__)

    @property
    def qupp(self):
        """ Upper state quantum numbers (list of Rational) """
        x = []
        n = self.sizequpp
        for i in range(n):
            x.append(Rational(c.c_void_p(lib.getelemUpperQuantumNumbersAbsorptionSingleLine(i, self.__data__))))
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
            x.append(Rational(c.c_void_p(lib.getelemLowerQuantumNumbersAbsorptionSingleLine(i, self.__data__))))
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

    def __eq__(self, other):
        if isinstance(other, AbsorptionSingleLine) and \
                (self.f0 == other.f0 or (isnan(self.f0 and isnan(other.f0)))) and \
                (self.i0 == other.i0 or (isnan(self.i0 and isnan(other.i0)))) and \
                (self.e0 == other.e0 or (isnan(self.e0 and isnan(other.e0)))) and \
                (self.gl == other.gl or (isnan(self.gl and isnan(other.gl)))) and \
                (self.gu == other.gu or (isnan(self.gu and isnan(other.gu)))) and \
                (self.a == other.a or (isnan(self.a and isnan(other.a)))) and \
                self.zeeman == other.zeeman and \
                self.lsm == other.lsm and \
                self.qupp == other.qupp and \
                self.qlow == other.qlow:
            return True
        else:
            return False


lib.createAbsorptionSingleLine.restype = c.c_void_p
lib.createAbsorptionSingleLine.argtypes = []

lib.deleteAbsorptionSingleLine.restype = None
lib.deleteAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.printAbsorptionSingleLine.restype = None
lib.printAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getF0AbsorptionSingleLine.restype = c.c_double
lib.getF0AbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getI0AbsorptionSingleLine.restype = c.c_double
lib.getI0AbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getE0AbsorptionSingleLine.restype = c.c_double
lib.getE0AbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getg_lowAbsorptionSingleLine.restype = c.c_double
lib.getg_lowAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getg_uppAbsorptionSingleLine.restype = c.c_double
lib.getg_uppAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getAAbsorptionSingleLine.restype = c.c_double
lib.getAAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.setF0AbsorptionSingleLine.restype = None
lib.setF0AbsorptionSingleLine.argtypes = [c.c_void_p, c.c_double]

lib.setI0AbsorptionSingleLine.restype = None
lib.setI0AbsorptionSingleLine.argtypes = [c.c_void_p, c.c_double]

lib.setE0AbsorptionSingleLine.restype = None
lib.setE0AbsorptionSingleLine.argtypes = [c.c_void_p, c.c_double]

lib.setg_lowAbsorptionSingleLine.restype = None
lib.setg_lowAbsorptionSingleLine.argtypes = [c.c_void_p, c.c_double]

lib.setg_uppAbsorptionSingleLine.restype = None
lib.setg_uppAbsorptionSingleLine.argtypes = [c.c_void_p, c.c_double]

lib.setAAbsorptionSingleLine.restype = None
lib.setAAbsorptionSingleLine.argtypes = [c.c_void_p, c.c_double]

lib.getZeemanAbsorptionSingleLine.restype = c.c_void_p
lib.getZeemanAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.getLineShapeAbsorptionSingleLine.restype = c.c_void_p
lib.getLineShapeAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.sizeLowerQuantumNumbersAbsorptionSingleLine.restype = c.c_long
lib.sizeLowerQuantumNumbersAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.sizeUpperQuantumNumbersAbsorptionSingleLine.restype = c.c_long
lib.sizeUpperQuantumNumbersAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.resizeLowerQuantumNumbersAbsorptionSingleLine.restype = None
lib.resizeLowerQuantumNumbersAbsorptionSingleLine.argtypes = [c.c_long, c.c_void_p]

lib.resizeUpperQuantumNumbersAbsorptionSingleLine.restype = None
lib.resizeUpperQuantumNumbersAbsorptionSingleLine.argtypes = [c.c_long, c.c_void_p]

lib.getelemLowerQuantumNumbersAbsorptionSingleLine.restype = c.c_void_p
lib.getelemLowerQuantumNumbersAbsorptionSingleLine.argtypes = [c.c_long, c.c_void_p]

lib.getelemUpperQuantumNumbersAbsorptionSingleLine.restype = c.c_void_p
lib.getelemUpperQuantumNumbersAbsorptionSingleLine.argtypes = [c.c_long, c.c_void_p]
