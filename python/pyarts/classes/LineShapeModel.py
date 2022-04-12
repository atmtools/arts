import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.LineShapeSingleSpeciesModel import LineShapeSingleSpeciesModel
from pyarts.classes.BasicTypes import String


class LineShapeType:
    """ ARTS LineShapeType data
    
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
            self.__data__ = c.c_void_p(lib.createLineShapeType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getLineShapeTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setLineShapeTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad LineShapeType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printLineShapeType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteLineShapeType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, LineShapeType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, LineShapeType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createLineShapeType.restype = c.c_void_p
lib.createLineShapeType.argtypes = []

lib.deleteLineShapeType.restype = None
lib.deleteLineShapeType.argtypes = [c.c_void_p]

lib.printLineShapeType.restype = None
lib.printLineShapeType.argtypes = [c.c_void_p]

lib.getLineShapeTypeString.restype = c.c_void_p
lib.getLineShapeTypeString.argtypes = [c.c_void_p]

lib.setLineShapeTypeString.restype = c.c_int
lib.setLineShapeTypeString.argtypes = [c.c_void_p, c.c_char_p]


class LineShapeModel:
    """ ARTS LineShape::Model data

    Getter and setter for data-property is implemented

    Properties:
        size:
            Number of LineShapeSingleSpeciesModel(s) (Index)

        data:
            The data (list of LineShapeSingleSpeciesModel(s))
        """
    def __init__(self, data=[]):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createLineShapeModel())
            self.data = data

    @property
    def size(self):
        """ Number of LineShapeSingleSpeciesModel(s) (Index) """
        return lib.sizeLineShapeModel(self.__data__)
    @size.setter
    def size(self, size):
        size = int(size)
        if size < 0:
            raise ValueError("size must be 0 or larger")
        lib.resizeLineShapeModel(size, self.__data__)

    @property
    def data(self):
        """ The data (list of LineShapeSingleSpeciesModel(s)) """
        x = []
        n = self.size
        for i in range(n):
            x.append(self[i])
        return x

    @data.setter
    def data(self, val):
        if isinstance(val, Sized):
            self.size = len(val)
            n = self.size
            for i in range(n):
                self[i] = val[i]
        else:
            raise TypeError("Only accepts array-like input")

    def __getitem__(self, ind):
        if ind >= 0 and ind < self.size:
            return LineShapeSingleSpeciesModel(c.c_void_p(lib.getelemLineShapeModel(ind, self.__data__)))
        else:
            raise IndexError("Out of bounds")

    def __setitem__(self, ind, val):
        self[ind].set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printLineShapeModel(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteLineShapeModel(self.__data__)

    def __repr__(self):
        return "ARTS LineShape::Model with {} LineShape::SingleSpeciesModel(s)".format(self.size)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, LineShapeModel):
            self.data = other.data
        else:
            raise TypeError("Expects LineShapeModel")

    def __eq__(self, other):
        if isinstance(other, LineShapeModel) and self.data == other.data:
            return True
        else:
            return False


lib.createLineShapeModel.restype = c.c_void_p
lib.createLineShapeModel.argtypes = []

lib.deleteLineShapeModel.restype = None
lib.deleteLineShapeModel.argtypes = [c.c_void_p]

lib.printLineShapeModel.restype = None
lib.printLineShapeModel.argtypes = [c.c_void_p]

lib.getelemLineShapeModel.restype = c.c_void_p
lib.getelemLineShapeModel.argtypes = [c.c_long, c.c_void_p]

lib.sizeLineShapeModel.restype = c.c_long
lib.sizeLineShapeModel.argtypes = [c.c_void_p]

lib.resizeLineShapeModel.restype = None
lib.resizeLineShapeModel.argtypes = [c.c_long, c.c_void_p]
