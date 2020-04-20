import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import String
from pyarts.classes.Tensor5 import Tensor5
from pyarts.classes.Tensor7 import Tensor7
from pyarts.classes.Vector import Vector
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class SingleScatteringData:
    """ ARTS SingleScatteringData data

    Properties:
        ptype:
            (Index)  FIXME: NOT INDEX BUT ENUM, CAN I CHANGE C++ CODE?

        description:
            (String)

        f_grid:
            (Vector)

        T_grid:
            (Vector)

        za_grid:
            (Vector)

        aa_grid:
            (Vector)

        pha_mat_data:
            (Tensor7)

        ext_mat_data:
            (Tensor5)

        abs_vec_data:
            (Tensor5)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSingleScatteringData())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "SingleScatteringData"

    @property
    def ptype(self):
        """ (Index) """
        return lib.getptypeSingleScatteringData(self.__data__).contents.value

    @ptype.setter
    def ptype(self, val):
        lib.getptypeSingleScatteringData(self.__data__)[0] = c.c_int(int(val))

    @property
    def description(self):
        """ (String) """
        return String(c.c_void_p(lib.getdescriptionSingleScatteringData(self.__data__)))

    @description.setter
    def description(self, val):
        self.description.set(val)

    @property
    def f_grid(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getf_gridSingleScatteringData(self.__data__)))

    @f_grid.setter
    def f_grid(self, val):
        self.f_grid.set(val)

    @property
    def T_grid(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getT_gridSingleScatteringData(self.__data__)))

    @T_grid.setter
    def T_grid(self, val):
        self.T_grid.set(val)

    @property
    def za_grid(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getza_gridSingleScatteringData(self.__data__)))

    @za_grid.setter
    def za_grid(self, val):
        self.za_grid.set(val)

    @property
    def aa_grid(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getaa_gridSingleScatteringData(self.__data__)))

    @aa_grid.setter
    def aa_grid(self, val):
        self.aa_grid.set(val)

    @property
    def pha_mat_data(self):
        """ (Tensor7) """
        return Tensor7(c.c_void_p(lib.getpha_mat_dataSingleScatteringData(self.__data__)))

    @pha_mat_data.setter
    def pha_mat_data(self, val):
        self.pha_mat_data.set(val)

    @property
    def ext_mat_data(self):
        """ (Tensor5) """
        return Tensor5(c.c_void_p(lib.getext_mat_dataSingleScatteringData(self.__data__)))

    @ext_mat_data.setter
    def ext_mat_data(self, val):
        self.ext_mat_data.set(val)

    @property
    def abs_vec_data(self):
        """ (Tensor5) """
        return Tensor5(c.c_void_p(lib.getabs_vec_dataSingleScatteringData(self.__data__)))

    @abs_vec_data.setter
    def abs_vec_data(self, val):
        self.abs_vec_data.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSingleScatteringData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSingleScatteringData(self.__data__)

    def __repr__(self):
        return "ARTS SingleScatteringData"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, SingleScatteringData):
            self.ptype = other.ptype
            self.description = other.description
            self.f_grid = other.f_grid
            self.T_grid = other.T_grid
            self.za_grid = other.za_grid
            self.aa_grid = other.aa_grid
            self.pha_mat_data = other.pha_mat_data
            self.ext_mat_data = other.ext_mat_data
            self.abs_vec_data = other.abs_vec_data
        else:
            raise TypeError("Expects SingleScatteringData")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadSingleScatteringData(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveSingleScatteringData(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, SingleScatteringData) and \
                self.ptype == other.ptype and \
                self.description == other.description and \
                self.f_grid == other.f_grid and \
                self.T_grid == other.T_grid and \
                self.za_grid == other.za_grid and \
                self.aa_grid == other.aa_grid and \
                self.pha_mat_data == other.pha_mat_data and \
                self.ext_mat_data == other.ext_mat_data and \
                self.abs_vec_data == other.abs_vec_data:
            return True
        else:
            return False


exec(array_base(SingleScatteringData))


exec(array_base(ArrayOfSingleScatteringData))


lib.createSingleScatteringData.restype = c.c_void_p
lib.createSingleScatteringData.argtypes = []

lib.deleteSingleScatteringData.restype = None
lib.deleteSingleScatteringData.argtypes = [c.c_void_p]

lib.printSingleScatteringData.restype = None
lib.printSingleScatteringData.argtypes = [c.c_void_p]

lib.xmlreadSingleScatteringData.restype = c.c_long
lib.xmlreadSingleScatteringData.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveSingleScatteringData.restype = c.c_long
lib.xmlsaveSingleScatteringData.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getptypeSingleScatteringData.restype = c.POINTER(c.c_int)
lib.getptypeSingleScatteringData.argtypes = [c.c_void_p]

lib.getdescriptionSingleScatteringData.restype = c.c_void_p
lib.getdescriptionSingleScatteringData.argtypes = [c.c_void_p]

lib.getf_gridSingleScatteringData.restype = c.c_void_p
lib.getf_gridSingleScatteringData.argtypes = [c.c_void_p]

lib.getT_gridSingleScatteringData.restype = c.c_void_p
lib.getT_gridSingleScatteringData.argtypes = [c.c_void_p]

lib.getza_gridSingleScatteringData.restype = c.c_void_p
lib.getza_gridSingleScatteringData.argtypes = [c.c_void_p]

lib.getaa_gridSingleScatteringData.restype = c.c_void_p
lib.getaa_gridSingleScatteringData.argtypes = [c.c_void_p]

lib.getpha_mat_dataSingleScatteringData.restype = c.c_void_p
lib.getpha_mat_dataSingleScatteringData.argtypes = [c.c_void_p]

lib.getext_mat_dataSingleScatteringData.restype = c.c_void_p
lib.getext_mat_dataSingleScatteringData.argtypes = [c.c_void_p]

lib.getabs_vec_dataSingleScatteringData.restype = c.c_void_p
lib.getabs_vec_dataSingleScatteringData.argtypes = [c.c_void_p]
