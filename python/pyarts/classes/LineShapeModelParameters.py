import ctypes as c
from pyarts.workspace.api import arts_api as lib


class InternalLineShapeModelParameters(c.Structure):
    """ ARTS Internal data layout... do not use """
    _fields_ = [("type", c.c_long),
                ("x0", c.c_double),
                ("x1", c.c_double),
                ("x2", c.c_double),
                ("x3", c.c_double),]

    def __init__(self, type=0, x0=float('nan'), x1=float('nan'), x2=float('nan'), x3=float('nan')):
        self.type = type
        self.x0 = c.c_double(x0)
        self.x1 = c.c_double(x1)
        self.x2 = c.c_double(x2)
        self.x3 = c.c_double(x3)

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
    def __init__(self, type="#", x0=float('nan'), x1=float('nan'), x2=float('nan'), x3=float('nan')):
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
        self.__data__.x0 = c.c_double(val)

    @property
    def x1(self):
        """ Model parameter #2 (Numeric) """
        return self.__data__.x1

    @x1.setter
    def x1(self, val):
        self.__data__.x1 = c.c_double(val)

    @property
    def x2(self):
        """ Model parameter #3 (Numeric) """
        return self.__data__.x2

    @x2.setter
    def x2(self, val):
        self.__data__.x2 = c.c_double(val)

    @property
    def x3(self):
        """ Model parameter #4 (Numeric) """
        return self.__data__.x3

    @x3.setter
    def x3(self, val):
        self.__data__.x3 = c.c_double(val)

    @property
    def data(self):
        """ Full data representation (InternalLineShapeModelParameters) """
        return self.__data__

    def print(self):
        """ Print to cout the ARTS representation of the class """
        self.__data__.print()

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


lib.printLineShapeModelParameters.restype = None
lib.printLineShapeModelParameters.argtypes = [c.c_void_p]

lib.getLineShapeModelParametersType.restype = c.c_long
lib.getLineShapeModelParametersType.argtypes = [c.c_char_p]
