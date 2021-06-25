import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import String, Numeric

from math import isnan, nan


class LineShapeTemperatureModel:
    """ ARTS LineShapeTemperatureModel data
    
    Properties:
        name:
            Name of type (String)
    """
    def __init__(self, data):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createLineShapeTemperatureModel())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getLineShapeTemperatureModelString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setLineShapeTemperatureModelString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad LineShapeTemperatureModel: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printLineShapeTemperatureModel(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteLineShapeTemperatureModel(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, LineShapeTemperatureModel) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, LineShapeTemperatureModel) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createLineShapeTemperatureModel.restype = c.c_void_p
lib.createLineShapeTemperatureModel.argtypes = []

lib.deleteLineShapeTemperatureModel.restype = None
lib.deleteLineShapeTemperatureModel.argtypes = [c.c_void_p]

lib.printLineShapeTemperatureModel.restype = None
lib.printLineShapeTemperatureModel.argtypes = [c.c_void_p]

lib.getLineShapeTemperatureModelString.restype = c.c_void_p
lib.getLineShapeTemperatureModelString.argtypes = [c.c_void_p]

lib.setLineShapeTemperatureModelString.restype = c.c_int
lib.setLineShapeTemperatureModelString.argtypes = [c.c_void_p, c.c_char_p]


class LineShapeModelParameters:
    """ ARTS LineShape::ModelParameters data

    Properties:
        type:
            Type of model (LineShapeTemperatureModel)

        x0:
            Model parameter #1 (Numeric)

        x1:
            Model parameter #2 (Numeric)

        x2:
            Model parameter #3 (Numeric)

        x3:
            Model parameter #4 (Numeric)
        """
    def __init__(self, type="None", x0=nan, x1=nan, x2=nan, x3=nan):
        if isinstance(type, c.c_void_p):
            self.__delete__ = False
            self.__data__ = type
        else:
            self.__data__ = c.c_void_p(lib.createLineShapeModelParameters())
            self.type = type
            self.x0 = x0
            self.x1 = x1
            self.x2 = x2
            self.x3 = x3

    @property
    def type(self):
        """ Type of model (LineShapeTemperatureModel) """
        return LineShapeTemperatureModel(c.c_void_p(lib.gettypeLineShapeModelParameters(self.__data__)))

    @type.setter
    def type(self, val):
        self.type.set(val)

    @property
    def x0(self):
        """ Model parameter #1 (Numeric) """
        return Numeric(c.c_void_p(lib.getX0LineShapeModelParameters(self.__data__)))

    @x0.setter
    def x0(self, val):
        self.x0.set(val)

    @property
    def x1(self):
        """ Model parameter #2 (Numeric) """
        return Numeric(c.c_void_p(lib.getX1LineShapeModelParameters(self.__data__)))

    @x1.setter
    def x1(self, val):
        self.x1.set(val)

    @property
    def x2(self):
        """ Model parameter #3 (Numeric) """
        return Numeric(c.c_void_p(lib.getX2LineShapeModelParameters(self.__data__)))

    @x2.setter
    def x2(self, val):
        self.x2.set(val)

    @property
    def x3(self):
        """ Model parameter #4 (Numeric) """
        return Numeric(c.c_void_p(lib.getX3LineShapeModelParameters(self.__data__)))

    @x3.setter
    def x3(self, val):
        self.x3.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        self.__data__.print()

    def __repr__(self):
        return f"{self.type} {self.x0} {self.x1} {self.x2} {self.x3}"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, LineShapeModelParameters):
            self.__data__.type = other.type
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


lib.createLineShapeModelParameters.restype = c.c_void_p
lib.createLineShapeModelParameters.argtypes = []

lib.deleteLineShapeModelParameters.restype = None
lib.deleteLineShapeModelParameters.argtypes = [c.c_void_p]

lib.printLineShapeModelParameters.restype = None
lib.printLineShapeModelParameters.argtypes = [c.c_void_p]

lib.gettypeLineShapeModelParameters.restype = c.c_void_p
lib.gettypeLineShapeModelParameters.argtypes = [c.c_void_p]

lib.getX0LineShapeModelParameters.restype = c.c_void_p
lib.getX0LineShapeModelParameters.argtypes = [c.c_void_p]

lib.getX1LineShapeModelParameters.restype = c.c_void_p
lib.getX1LineShapeModelParameters.argtypes = [c.c_void_p]

lib.getX2LineShapeModelParameters.restype = c.c_void_p
lib.getX2LineShapeModelParameters.argtypes = [c.c_void_p]

lib.getX3LineShapeModelParameters.restype = c.c_void_p
lib.getX3LineShapeModelParameters.argtypes = [c.c_void_p]
