import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import String, Index
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes.JacobianTarget import JacobianTarget
from pyarts.classes.Matrix import Matrix
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class RetrievalQuantity:
    """ ARTS RetrievalQuantity data

    Properties:
        subtag:
            (String)

        subsubtag:
            (String)

        mode:
            (String)

        target:
            (JacobianTarget)

        grids:
            (ArrayOfVector)

        transformation_func:
            (String)

        t_func_parameters:
            (Vector)

        transformation:
            (Matrix)

        offset:
            (Vector)

        integration:
            (bool)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createRetrievalQuantity())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def offset(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getOffsetRetrievalQuantity(self.__data__)))

    @offset.setter
    def offset(self, val):
        self.offset.set(val)

    @property
    def transformation(self):
        """ (Matrix) """
        return Matrix(c.c_void_p(lib.getTransformationRetrievalQuantity(self.__data__)))

    @transformation.setter
    def transformation(self, val):
        self.transformation.set(val)

    @property
    def t_func_parameters(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getTFuncParametersRetrievalQuantity(self.__data__)))

    @t_func_parameters.setter
    def t_func_parameters(self, val):
        self.t_func_parameters.set(val)

    @property
    def transformation_func(self):
        """ (String) """
        return String(c.c_void_p(lib.getTransformationFuncRetrievalQuantity(self.__data__)))

    @transformation_func.setter
    def transformation_func(self, val):
        self.transformation_func.set(val)
        
    @property
    def grids(self):
        """ (ArrayOfVector) """
        return ArrayOfVector(c.c_void_p(lib.getGridsRetrievalQuantity(self.__data__)))

    @grids.setter
    def grids(self, val):
        self.grids.set(val)

    @property
    def mode(self):
        """ (String) """
        return String(c.c_void_p(lib.getModeRetrievalQuantity(self.__data__)))

    @mode.setter
    def mode(self, val):
        self.mode.set(val)

    @property
    def subsubtag(self):
        """ (String) """
        return String(c.c_void_p(lib.getSubSubTagRetrievalQuantity(self.__data__)))

    @subsubtag.setter
    def subsubtag(self, val):
        self.subsubtag.set(val)

    @property
    def subtag(self):
        """ (String) """
        return String(c.c_void_p(lib.getSubTagRetrievalQuantity(self.__data__)))

    @subtag.setter
    def subtag(self, val):
        self.subtag.set(val)

    @property
    def target(self):
        """ (JacobianTarget) """
        return JacobianTarget(c.c_void_p(lib.getTargetRetrievalQuantity(self.__data__)))

    @target.setter
    def target(self, val):
        self.target.set(val)

    @staticmethod
    def name():
        return "RetrievalQuantity"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printRetrievalQuantity(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteRetrievalQuantity(self.__data__)

    def __repr__(self):
        return "ARTS RetrievalQuantity"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, RetrievalQuantity):
            self.subtag = other.subtag
            self.subsubtag = other.subsubtag
            self.mode = other.mode
            self.target = other.target
            self.grids = other.grids
            self.transformation_func = other.transformation_func
            self.t_func_parameters = other.t_func_parameters
            self.transformation = other.transformation
            self.offset = other.offset
        else:
            raise TypeError("Expects RetrievalQuantity")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadRetrievalQuantity(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveRetrievalQuantity(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, RetrievalQuantity) and \
                self.subtag == other.subtag and \
                self.subsubtag == other.subsubtag and \
                self.mode == other.mode and \
                self.target == other.target and \
                self.grids == other.grids and \
                self.transformation_func == other.transformation_func and \
                self.t_func_parameters == other.t_func_parameters and \
                self.transformation == other.transformation and \
                self.offset == other.offset:
            return True
        else:
            return False


exec(array_base(RetrievalQuantity))


lib.createRetrievalQuantity.restype = c.c_void_p
lib.createRetrievalQuantity.argtypes = []

lib.deleteRetrievalQuantity.restype = None
lib.deleteRetrievalQuantity.argtypes = [c.c_void_p]

lib.printRetrievalQuantity.restype = None
lib.printRetrievalQuantity.argtypes = [c.c_void_p]

lib.xmlreadRetrievalQuantity.restype = c.c_long
lib.xmlreadRetrievalQuantity.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveRetrievalQuantity.restype = c.c_long
lib.xmlsaveRetrievalQuantity.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getSubTagRetrievalQuantity.restype = c.c_void_p
lib.getSubTagRetrievalQuantity.argtypes = [c.c_void_p]

lib.getSubSubTagRetrievalQuantity.restype = c.c_void_p
lib.getSubSubTagRetrievalQuantity.argtypes = [c.c_void_p]

lib.getModeRetrievalQuantity.restype = c.c_void_p
lib.getModeRetrievalQuantity.argtypes = [c.c_void_p]

lib.getGridsRetrievalQuantity.restype = c.c_void_p
lib.getGridsRetrievalQuantity.argtypes = [c.c_void_p]

lib.getTargetRetrievalQuantity.restype = c.c_void_p
lib.getTargetRetrievalQuantity.argtypes = [c.c_void_p]

lib.getTransformationFuncRetrievalQuantity.restype = c.c_void_p
lib.getTransformationFuncRetrievalQuantity.argtypes = [c.c_void_p]

lib.getTFuncParametersRetrievalQuantity.restype = c.c_void_p
lib.getTFuncParametersRetrievalQuantity.argtypes = [c.c_void_p]

lib.getTransformationRetrievalQuantity.restype = c.c_void_p
lib.getTransformationRetrievalQuantity.argtypes = [c.c_void_p]

lib.getOffsetRetrievalQuantity.restype = c.c_void_p
lib.getOffsetRetrievalQuantity.argtypes = [c.c_void_p]
