import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Numeric, String
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


from math import isnan

class JacobianTarget:
    """ ARTS JacobianTarget data

    Properties:
        type:
            Enumeration type of Jacobian (String)

        subtype:
            Enumerated subtype of Jacobian (String)

        perturbation:
            Perturbation for some calculations (Numeric)

        quantumidentity:
            Identifier of the Jacobian line parameter (QuantumIdentifier)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createJacobianTarget())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def type(self):
        """ Enumeration type of Jacobian (String) """
        return String(c.c_void_p(lib.enumgetTargetTypeJacobianTarget(self.__data__)), delete=True)

    @type.setter
    def type(self, val): 
        if lib.enumsetTargetTypeJacobianTarget(self.__data__, str(val).encode("utf-8")):
            raise ValueError("Bad input, class in undefined mode")

    @property
    def subtype(self):
        """ Enumeration subtype of Jacobian (String) """
        return String(c.c_void_p(lib.enumgetTargetSubTypeJacobianTarget(self.__data__)), delete=True)

    @subtype.setter
    def subtype(self, val): 
        if lib.enumsetTargetSubTypeJacobianTarget(self.__data__, str(val).encode("utf-8")):
            raise ValueError("Bad input, class in undefined mode")

    @property
    def quantumidentity(self):
        """ Identifier of the Jacobian line parameter (QuantumIdentifier) """
        return QuantumIdentifier(c.c_void_p(lib.getQuantumIdentityJacobianTarget(self.__data__)))

    @quantumidentity.setter
    def quantumidentity(self, val):
        self.quantumidentity.set(val)

    @property
    def perturbation(self):
        """ Perturbation for some calculations (Numeric) """
        return Numeric(c.c_void_p(lib.getPerturbationJacobianTarget(self.__data__)))

    @perturbation.setter
    def perturbation(self, val):
        self.perturbation.set(val)

    @staticmethod
    def name():
        return "JacobianTarget"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printJacobianTarget(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteJacobianTarget(self.__data__)

    def __repr__(self):
        return "ARTS JacobianTarget"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, JacobianTarget):
            self.type = other.type
            self.subtype = other.subtype
            self.perturbation = other.perturbation
            self.quantumidentity = other.quantumidentity
        else:
            raise TypeError("Expects JacobianTarget")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadJacobianTarget(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveJacobianTarget(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, JacobianTarget) and \
                self.type == other.type and\
                self.subtype == other.subtype and \
                (self.perturbation == other.perturbation or (isnan(self.perturbation) and isnan(other.perturbation))) and \
                self.quantumidentity == other.quantumidentity:
            return True
        else:
            return False


exec(array_base(JacobianTarget))


lib.createJacobianTarget.restype = c.c_void_p
lib.createJacobianTarget.argtypes = []

lib.deleteJacobianTarget.restype = None
lib.deleteJacobianTarget.argtypes = [c.c_void_p]

lib.printJacobianTarget.restype = None
lib.printJacobianTarget.argtypes = [c.c_void_p]

lib.xmlreadJacobianTarget.restype = c.c_long
lib.xmlreadJacobianTarget.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveJacobianTarget.restype = c.c_long
lib.xmlsaveJacobianTarget.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.enumgetTargetTypeJacobianTarget.restype = c.c_void_p
lib.enumgetTargetTypeJacobianTarget.argtypes = [c.c_void_p]

lib.enumsetTargetTypeJacobianTarget.restype = c.c_bool
lib.enumsetTargetTypeJacobianTarget.argtypes = [c.c_void_p, c.c_char_p]

lib.enumgetTargetSubTypeJacobianTarget.restype = c.c_void_p
lib.enumgetTargetSubTypeJacobianTarget.argtypes = [c.c_void_p]

lib.enumsetTargetSubTypeJacobianTarget.restype = c.c_bool
lib.enumsetTargetSubTypeJacobianTarget.argtypes = [c.c_void_p, c.c_char_p]

lib.getPerturbationJacobianTarget.restype = c.c_void_p
lib.getPerturbationJacobianTarget.argtypes = [c.c_void_p]

lib.getQuantumIdentityJacobianTarget.restype = c.c_void_p
lib.getQuantumIdentityJacobianTarget.argtypes = [c.c_void_p]
