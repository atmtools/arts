import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Numeric, String
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class ScatteringMetaData:
    """ ARTS ScatteringMetaData data

    Properties:
        description:
            (String)

        source:
            (String)

        refr_index:
            (String)

        mass:
            (Numeric)

        diameter_max:
            (Numeric)

        diameter_volume_equ:
            (Numeric)

        diameter_area_equ_aerodynamical:
            (Numeric)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createScatteringMetaData())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "ScatteringMetaData"

    @property
    def description(self):
        """ (String) """
        return String(c.c_void_p(lib.getdescriptionScatteringMetaData(self.__data__)))

    @description.setter
    def description(self, val):
        self.description.set(val)

    @property
    def source(self):
        """ (String) """
        return String(c.c_void_p(lib.getsourceScatteringMetaData(self.__data__)))

    @source.setter
    def source(self, val):
        self.source.set(val)

    @property
    def refr_index(self):
        """ (String) """
        return String(c.c_void_p(lib.getrefr_indexScatteringMetaData(self.__data__)))

    @refr_index.setter
    def refr_index(self, val):
        self.refr_index.set(val)

    @property
    def mass(self):
        """ (Numeric) """
        return Numeric(c.c_void_p(lib.getmassScatteringMetaData(self.__data__)))

    @mass.setter
    def mass(self, val):
        self.mass.set(val)

    @property
    def diameter_max(self):
        """ (Numeric) """
        return Numeric(c.c_void_p(lib.getdiameter_maxScatteringMetaData(self.__data__)))

    @diameter_max.setter
    def diameter_max(self, val):
        self.diameter_max.set(val)

    @property
    def diameter_volume_equ(self):
        """ (Numeric) """
        return Numeric(c.c_void_p(lib.getdiameter_volume_equScatteringMetaData(self.__data__)))

    @diameter_volume_equ.setter
    def diameter_volume_equ(self, val):
        self.diameter_volume_equ.set(val)

    @property
    def diameter_area_equ_aerodynamical(self):
        """ (Numeric) """
        return Numeric(c.c_void_p(lib.getdiameter_area_equ_aerodynamicalScatteringMetaData(self.__data__)))

    @diameter_area_equ_aerodynamical.setter
    def diameter_area_equ_aerodynamical(self, val):
        self.diameter_area_equ_aerodynamical.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printScatteringMetaData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteScatteringMetaData(self.__data__)

    def __repr__(self):
        return "ARTS ScatteringMetaData"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, ScatteringMetaData):
            self.description = other.description
            self.source = other.source
            self.refr_index = other.refr_index
            self.mass = other.mass
            self.diameter_max = other.diameter_max
            self.diameter_volume_equ = other.diameter_volume_equ
            self.diameter_area_equ_aerodynamical = other.diameter_area_equ_aerodynamical
        else:
            raise TypeError("Expects ScatteringMetaData")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadScatteringMetaData(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveScatteringMetaData(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


exec(array_base(ScatteringMetaData))


exec(array_base(ArrayOfScatteringMetaData))


lib.createScatteringMetaData.restype = c.c_void_p
lib.createScatteringMetaData.argtypes = []

lib.deleteScatteringMetaData.restype = None
lib.deleteScatteringMetaData.argtypes = [c.c_void_p]

lib.printScatteringMetaData.restype = None
lib.printScatteringMetaData.argtypes = [c.c_void_p]

lib.xmlreadScatteringMetaData.restype = c.c_long
lib.xmlreadScatteringMetaData.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveScatteringMetaData.restype = c.c_long
lib.xmlsaveScatteringMetaData.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getdescriptionScatteringMetaData.restype = c.c_void_p
lib.getdescriptionScatteringMetaData.argtypes = [c.c_void_p]

lib.getsourceScatteringMetaData.restype = c.c_void_p
lib.getsourceScatteringMetaData.argtypes = [c.c_void_p]

lib.getrefr_indexScatteringMetaData.restype = c.c_void_p
lib.getrefr_indexScatteringMetaData.argtypes = [c.c_void_p]

lib.getmassScatteringMetaData.restype = c.c_void_p
lib.getmassScatteringMetaData.argtypes = [c.c_void_p]

lib.getdiameter_maxScatteringMetaData.restype = c.c_void_p
lib.getdiameter_maxScatteringMetaData.argtypes = [c.c_void_p]

lib.getdiameter_volume_equScatteringMetaData.restype = c.c_void_p
lib.getdiameter_volume_equScatteringMetaData.argtypes = [c.c_void_p]

lib.getdiameter_area_equ_aerodynamicalScatteringMetaData.restype = c.c_void_p
lib.getdiameter_area_equ_aerodynamicalScatteringMetaData.argtypes = [c.c_void_p]
