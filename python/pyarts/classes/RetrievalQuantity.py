import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import String, Numeric, Index
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes.Matrix import Matrix
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class RetrievalQuantity:
    """ ARTS RetrievalQuantity data

    Properties:
        maintag:
            (String)

        subtag:
            (String)

        subsubtag:
            (String)

        mode:
            (String)

        analytical:
            (Index)

        perturbation:
            (Numeric)

        grids:
            (ArrayOfVector)

        quantumidentity:
            (QuantumIdentifier)

        transformation_func:
            (String)

        t_func_parameters:
            (Vector)

        transformation:
            (Matrix)

        offset:
            (Vector)

        type:
            (Index)

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
    def type(self):
        """ (Index) """
        return lib.getTypeRetrievalQuantity(self.__data__)

    @type.setter
    def type(self, val):
        if lib.setTypeRetrievalQuantity(self.__data__, int(val)):
            raise ValueError("Bad input")

    @property
    def integration(self):
        """ (bool) """
        return lib.getIntegrationRetrievalQuantity(self.__data__)

    @integration.setter
    def integration(self, val):
        lib.setIntegrationRetrievalQuantity(self.__data__, bool(val))

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
    def quantumidentity(self):
        """ (QuantumIdentifier) """
        return QuantumIdentifier(c.c_void_p(lib.getQuantumIdentityRetrievalQuantity(self.__data__)))

    @quantumidentity.setter
    def quantumidentity(self, val):
        self.quantumidentity.set(val)

    @property
    def grids(self):
        """ (ArrayOfVector) """
        return ArrayOfVector(c.c_void_p(lib.getGridsRetrievalQuantity(self.__data__)))

    @grids.setter
    def grids(self, val):
        self.grids.set(val)

    @property
    def perturbation(self):
        """ (Numeric) """
        return Numeric(c.c_void_p(lib.getPerturbationRetrievalQuantity(self.__data__)))

    @perturbation.setter
    def perturbation(self, val):
        self.perturbation.set(val)

    @property
    def analytical(self):
        """ (Index) """
        return Index(c.c_void_p(lib.getAnalyticalRetrievalQuantity(self.__data__)))

    @analytical.setter
    def analytical(self, val):
        self.analytical.set(val)

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
    def maintag(self):
        """ (String) """
        return String(c.c_void_p(lib.getMainTagRetrievalQuantity(self.__data__)))

    @maintag.setter
    def maintag(self, val):
        self.maintag.set(val)

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
              raise RuntimeWarning("Cannot set RetrievalQuantity, remains constant")
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

lib.getTypeRetrievalQuantity.restype = c.c_long
lib.getTypeRetrievalQuantity.argtypes = [c.c_void_p]

lib.setTypeRetrievalQuantity.restype = c.c_long
lib.setTypeRetrievalQuantity.argtypes = [c.c_void_p, c.c_long]

lib.getIntegrationRetrievalQuantity.restype = c.c_long
lib.getIntegrationRetrievalQuantity.argtypes = [c.c_void_p]

lib.setIntegrationRetrievalQuantity.restype = c.c_long
lib.setIntegrationRetrievalQuantity.argtypes = [c.c_void_p, c.c_long]

lib.getMainTagRetrievalQuantity.restype = c.c_void_p
lib.getMainTagRetrievalQuantity.argtypes = [c.c_void_p]

lib.getSubTagRetrievalQuantity.restype = c.c_void_p
lib.getSubTagRetrievalQuantity.argtypes = [c.c_void_p]

lib.getSubSubTagRetrievalQuantity.restype = c.c_void_p
lib.getSubSubTagRetrievalQuantity.argtypes = [c.c_void_p]

lib.getModeRetrievalQuantity.restype = c.c_void_p
lib.getModeRetrievalQuantity.argtypes = [c.c_void_p]

lib.getAnalyticalRetrievalQuantity.restype = c.c_void_p
lib.getAnalyticalRetrievalQuantity.argtypes = [c.c_void_p]

lib.getPerturbationRetrievalQuantity.restype = c.c_void_p
lib.getPerturbationRetrievalQuantity.argtypes = [c.c_void_p]

lib.getGridsRetrievalQuantity.restype = c.c_void_p
lib.getGridsRetrievalQuantity.argtypes = [c.c_void_p]

lib.getQuantumIdentityRetrievalQuantity.restype = c.c_void_p
lib.getQuantumIdentityRetrievalQuantity.argtypes = [c.c_void_p]

lib.getTransformationFuncRetrievalQuantity.restype = c.c_void_p
lib.getTransformationFuncRetrievalQuantity.argtypes = [c.c_void_p]

lib.getTFuncParametersRetrievalQuantity.restype = c.c_void_p
lib.getTFuncParametersRetrievalQuantity.argtypes = [c.c_void_p]

lib.getTransformationRetrievalQuantity.restype = c.c_void_p
lib.getTransformationRetrievalQuantity.argtypes = [c.c_void_p]

lib.getOffsetRetrievalQuantity.restype = c.c_void_p
lib.getOffsetRetrievalQuantity.argtypes = [c.c_void_p]
