import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Index
from pyarts.classes.Matrix import Matrix
from pyarts.classes.Vector import Vector
from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class TessemNN:
    """ ARTS TessemNN data

    Properties:
        nb_inputs:
            (Index)

        nb_outputs:
            (Index)

        nb_cache:
            (Numeric)

        b1:
            (Vector)

        b2:
            (Vector)

        w1:
            (Matrix)

        w2:
            (Matrix)

        x_min:
            (Vector)

        x_max:
            (Vector)

        y_min:
            (Vector)

        y_max:
            (Vector)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTessemNN())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "TessemNN"

    @property
    def nb_inputs(self):
        """ (Index) """
        return Index(c.c_void_p(lib.getnb_inputsTessemNN(self.__data__)))

    @nb_inputs.setter
    def nb_inputs(self, val):
        self.nb_inputs.set(val)

    @property
    def nb_outputs(self):
        """ (Index) """
        return Index(c.c_void_p(lib.getnb_outputsTessemNN(self.__data__)))

    @nb_outputs.setter
    def nb_outputs(self, val):
        self.nb_outputs.set(val)

    @property
    def nb_cache(self):
        """ (Index) """
        return Index(c.c_void_p(lib.getnb_cacheTessemNN(self.__data__)))

    @nb_cache.setter
    def nb_cache(self, val):
        self.nb_cache.set(val)

    @property
    def b1(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getb1TessemNN(self.__data__)))

    @b1.setter
    def b1(self, val):
        self.b1.set(val)

    @property
    def b2(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getb2TessemNN(self.__data__)))

    @b2.setter
    def b2(self, val):
        self.b2.set(val)

    @property
    def x_min(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getx_minTessemNN(self.__data__)))

    @x_min.setter
    def x_min(self, val):
        self.x_min.set(val)

    @property
    def x_max(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getx_maxTessemNN(self.__data__)))

    @x_max.setter
    def x_max(self, val):
        self.x_max.set(val)

    @property
    def y_min(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.gety_minTessemNN(self.__data__)))

    @y_min.setter
    def y_min(self, val):
        self.y_min.set(val)

    @property
    def y_max(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.gety_maxTessemNN(self.__data__)))

    @y_max.setter
    def y_max(self, val):
        self.y_max.set(val)

    @property
    def w1(self):
        """ (Matrix) """
        return Matrix(c.c_void_p(lib.getw1TessemNN(self.__data__)))

    @w1.setter
    def w1(self, val):
        self.w1.set(val)

    @property
    def w2(self):
        """ (Matrix) """
        return Matrix(c.c_void_p(lib.getw2TessemNN(self.__data__)))

    @w2.setter
    def w2(self, val):
        self.w2.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTessemNN(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTessemNN(self.__data__)

    def __repr__(self):
        return "ARTS TessemNN"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, TessemNN):
              self.nb_inputs = other.nb_inputs
              self.nb_outputs = other.nb_outputs
              self.nb_cache = other.nb_cache
              self.b1 = other.b1
              self.b2 = other.b2
              self.w1 = other.w1
              self.w2 = other.w2
              self.x_min = other.x_min
              self.x_max = other.x_max
              self.y_min = other.y_min
              self.y_max = other.y_max
        else:
            raise TypeError("Expects TessemNN")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTessemNN(self.__data__, correct_read_arguments(file)):
            raise OSError("Cannot read {}".format(file))

    def savexml(self, file, type="ascii", clobber=True):
        """ Saves the class to XML file

        Input:
            file:
                Filename to writable file (str)

            type:
                Filetype (str)

            clobber:
                Allow clobbering files? (any boolean)
        """
        if lib.xmlsaveTessemNN(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createTessemNN.restype = c.c_void_p
lib.createTessemNN.argtypes = []

lib.deleteTessemNN.restype = None
lib.deleteTessemNN.argtypes = [c.c_void_p]

lib.printTessemNN.restype = None
lib.printTessemNN.argtypes = [c.c_void_p]

lib.xmlreadTessemNN.restype = c.c_long
lib.xmlreadTessemNN.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTessemNN.restype = c.c_long
lib.xmlsaveTessemNN.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getnb_inputsTessemNN.restype = c.c_void_p
lib.getnb_inputsTessemNN.argtypes = [c.c_void_p]

lib.getnb_outputsTessemNN.restype = c.c_void_p
lib.getnb_outputsTessemNN.argtypes = [c.c_void_p]

lib.getnb_cacheTessemNN.restype = c.c_void_p
lib.getnb_cacheTessemNN.argtypes = [c.c_void_p]

lib.getb1TessemNN.restype = c.c_void_p
lib.getb1TessemNN.argtypes = [c.c_void_p]

lib.getb2TessemNN.restype = c.c_void_p
lib.getb2TessemNN.argtypes = [c.c_void_p]

lib.getw1TessemNN.restype = c.c_void_p
lib.getw1TessemNN.argtypes = [c.c_void_p]

lib.getw2TessemNN.restype = c.c_void_p
lib.getw2TessemNN.argtypes = [c.c_void_p]

lib.getx_minTessemNN.restype = c.c_void_p
lib.getx_minTessemNN.argtypes = [c.c_void_p]

lib.getx_maxTessemNN.restype = c.c_void_p
lib.getx_maxTessemNN.argtypes = [c.c_void_p]

lib.gety_minTessemNN.restype = c.c_void_p
lib.gety_minTessemNN.argtypes = [c.c_void_p]

lib.gety_maxTessemNN.restype = c.c_void_p
lib.gety_maxTessemNN.argtypes = [c.c_void_p]
