import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib
from pyarts.classes.ArrayBase import array_base

from pyarts.classes.LineShapeModel import LineShapeModel
from pyarts.classes.Rational import Rational
from pyarts.classes.ZeemanModel import ZeemanModel
from pyarts.classes.BasicTypes import Numeric
from pyarts.classes.quantum import QuantumNumberLocalState

from pyarts.classes.macros import BasicInterfaceCAPI, EnumMacroInterfaceCAP, VoidStructGetterCAPI

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
        return Numeric(c.c_void_p(lib.getF0AbsorptionSingleLine(self.__data__)))

    @f0.setter
    def f0(self, x):
        self.f0.set(x)

    @property
    def i0(self):
        """ Line strength (Numeric) """
        return Numeric(c.c_void_p(lib.getI0AbsorptionSingleLine(self.__data__)))

    @i0.setter
    def i0(self, x):
        self.i0.set(x)

    @property
    def e0(self):
        """ Line lower state energy (Numeric) """
        return Numeric(c.c_void_p(lib.getE0AbsorptionSingleLine(self.__data__)))

    @e0.setter
    def e0(self, x):
        self.e0.set(x)

    @property
    def gl(self):
        """ Lower state statistical weight (Numeric) """
        return Numeric(c.c_void_p(lib.getglowAbsorptionSingleLine(self.__data__)))

    @gl.setter
    def gl(self, x):
        self.gl.set(x)

    @property
    def gu(self):
        """ Upper state statistical weight (Numeric) """
        return Numeric(c.c_void_p(lib.getguppAbsorptionSingleLine(self.__data__)))

    @gu.setter
    def gu(self, x):
        self.gu.set(x)

    @property
    def a(self):
        """ Einstein coefficient (Numeric) """
        return Numeric(c.c_void_p(lib.getAAbsorptionSingleLine(self.__data__)))

    @a.setter
    def a(self, x):
        self.a.set(x)

    @property
    def zeeman(self):
        """ Zeeman model data (ZeemanModel) """
        return ZeemanModel(c.c_void_p(lib.getzeemanAbsorptionSingleLine(self.__data__)))

    @zeeman.setter
    def zeeman(self, x):
        self.zeeman.set(x)

    @property
    def lsm(self):
        """ Line shape model (LineShapeModel) """
        return LineShapeModel(c.c_void_p(lib.getlineshapeAbsorptionSingleLine(self.__data__)))

    @lsm.setter
    def lsm(self, x):
        self.lsm.set(x)

    @property
    def localquanta(self):
        """ Upper state quantum numbers (QuantumNumberLocalState) """
        return QuantumNumberLocalState(c.c_void_p(lib.getlocalquantaAbsorptionSingleLine(self.__data__)))

    @localquanta.setter
    def localquanta(self, x):
        self.localquanta.set(x)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAbsorptionSingleLine(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionSingleLine(self.__data__)
    
    @staticmethod
    def name():
        return "AbsorptionSingleLine"

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
            self.localquanta = other.localquanta
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
                self.localquanta == other.localquanta:
            return True
        else:
            return False


exec(array_base(AbsorptionSingleLine))

lib.createAbsorptionSingleLine.restype = c.c_void_p
lib.createAbsorptionSingleLine.argtypes = []

lib.deleteAbsorptionSingleLine.restype = None
lib.deleteAbsorptionSingleLine.argtypes = [c.c_void_p]

lib.printAbsorptionSingleLine.restype = None
lib.printAbsorptionSingleLine.argtypes = [c.c_void_p]

VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "F0")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "I0")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "E0")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "glow")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "gupp")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "A")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "zeeman")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "lineshape")
VoidStructGetterCAPI(lib, "AbsorptionSingleLine", "localquanta")
