import ctypes as c
from pyarts.workspace.api import arts_api as lib

from math import isnan, nan


class InternalLineShapeModelParameters(c.Structure):
    """ ARTS Internal data layout... do not use """
    _fields_ = [("type", c.c_long),
                ("x0", c.c_double),
                ("x1", c.c_double),
                ("x2", c.c_double),
                ("x3", c.c_double),]

    def __init__(self, type=0, x0=nan, x1=nan, x2=nan, x3=nan):
        self.type = int(type)
        self.x0 = float(x0)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)

    def print(self):
        """ Print to cout the ARTS representation of the struct """
        lib.printLineShapeModelParameters(c.pointer(self))


class LineShapeModelParameters:
    """ ARTS LineShape::ModelParameters data

    Properties:
        type:
            Type of model (get: Index; set: str)

        x0:
            Model parameter #1 (Numeric)

        x1:
            Model parameter #2 (Numeric)

        x2:
            Model parameter #3 (Numeric)

        x3:
            Model parameter #4 (Numeric)
        """
    def __init__(self, type="#", x0=nan, x1=nan, x2=nan, x3=nan):
        if isinstance(type, InternalLineShapeModelParameters):
            self.__data__ = type
        else:
            self.__data__ = InternalLineShapeModelParameters()
            self.type = type
            self.x0 = x0
            self.x1 = x1
            self.x2 = x2
            self.x3 = x3

    @property
    def type(self):
        """ Type of model (get: Index; set: str) """
        return self.__data__.type

    @type.setter
    def type(self, val):
        val = str(val).encode("ascii")
        n = lib.getLineShapeModelParametersType(val)
        if n < 0:
            raise ValueError("Invalid type")
        self.__data__.type = n

    @property
    def x0(self):
        """ Model parameter #1 (Numeric) """
        return self.__data__.x0

    @x0.setter
    def x0(self, val):
        self.__data__.x0 = float(val)

    @property
    def x1(self):
        """ Model parameter #2 (Numeric) """
        return self.__data__.x1

    @x1.setter
    def x1(self, val):
        self.__data__.x1 = float(val)

    @property
    def x2(self):
        """ Model parameter #3 (Numeric) """
        return self.__data__.x2

    @x2.setter
    def x2(self, val):
        self.__data__.x2 = float(val)

    @property
    def x3(self):
        """ Model parameter #4 (Numeric) """
        return self.__data__.x3

    @x3.setter
    def x3(self, val):
        self.__data__.x3 = float(val)

    @property
    def data(self):
        """ Full data representation (InternalLineShapeModelParameters) """
        return self.__data__

    def print(self):
        """ Print to cout the ARTS representation of the class """
        self.__data__.print()

    def __repr__(self):
        return "ARTS LineShape::ModelParameters"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, LineShapeModelParameters):
            self.type = other.type
            self.x0 = other.x0
            self.x1 = other.x1
            self.x2 = other.x2
            self.x3 = other.x3
        else:
            raise TypeError("Expects LineShapeModelParameters")

    def __eq__(self, other):
        if isinstance(other, LineShapeModelParameters) and \
                self.type == other.type and \
                (self.x0 == other.x0 or (isnan(self.x0) and isnan(other.x0))) and \
                (self.x1 == other.x1 or (isnan(self.x1) and isnan(other.x1))) and \
                (self.x2 == other.x2 or (isnan(self.x2) and isnan(other.x2))) and \
                (self.x3 == other.x3 or (isnan(self.x3) and isnan(other.x3))):
            return True
        else:
            return False


lib.printLineShapeModelParameters.restype = None
lib.printLineShapeModelParameters.argtypes = [c.c_void_p]

lib.getLineShapeModelParametersType.restype = c.c_long
lib.getLineShapeModelParametersType.argtypes = [c.c_char_p]
