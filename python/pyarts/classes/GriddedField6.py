import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import ArrayOfString
from pyarts.classes.Vector import Vector
from pyarts.classes.Tensor6 import Tensor6
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class GriddedField6:
    """ ARTS GriddedField6 data

    Properties:
        dim:
            GriddedField dimension (const Index)

        fieldname:
            Name of field (String)

        gridtypes:
            List of gridtypes (const list of Index>; 0: Numeric grid, 1: String grid)

        shape:
            Size of the grids (tuple)

        gridnames:
            List of grid names (list of String)

        grids:
            List of grids (list of Vector and ArrayOfString)

        data:
            The data (Tensor6)
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createGriddedField6())

    @staticmethod
    def name():
        return "GriddedField6"

    @property
    def OK(self):
        return lib.checksizeGriddedField6(self.__data__)

    @property
    def dim(self):
        """ GriddedField dimension (const Index) """
        return lib.get_dimGriddedField6(self.__data__)

    @property
    def fieldname(self):
        """ Name of field (String) """
        return lib.get_nameGriddedField6(self.__data__).decode("utf-8")

    @fieldname.setter
    def fieldname(self, val):
        lib.set_nameGriddedField6(self.__data__, str(val).encode("ascii"))

    @property
    def gridtypes(self):
        """ List of gridtypes (list of Index>; 0: Numeric grid, 1: String grid) """
        x = []
        for i in range(self.dim):
            x.append(lib.get_grid_typeIndexGriddedField6(int(i), self.__data__))
        return x

    @property
    def gridnames(self):
        """ List of grid names (list of String) """
        x = []
        for i in range(self.dim):
            x.append(lib.get_grid_nameGriddedField6(int(i), self.__data__).decode("utf-8"))
        return x

    @gridnames.setter
    def gridnames(self, val):
        if isinstance(val, Sized) and len(val) == self.dim:
            for i in range(self.dim):
                if isinstance(val[i], str):
                    lib.set_grid_nameGriddedField6(int(i), self.__data__, val[i].encode("ascii"))
                else:
                    raise TypeError("Expects str input")
        else:
            raise TypeError("Only accepts array-like input of length {}".format(self.dim))

    @property
    def shape(self):
        """ Size of the grids (tuple) """
        x = []
        for i in range(self.dim):
            x.append(lib.get_grid_sizeGriddedField6(int(i), self.__data__))
        return tuple(x)

    @property
    def grids(self):
        types = self.gridtypes
        x = []
        for i in range(self.dim):
            if types[i] == 0:
                x.append(Vector(c.c_void_p(lib.get_numeric_gridGriddedField6(int(i), self.__data__))))
            else:
                x.append(ArrayOfString(c.c_void_p(lib.get_string_gridGriddedField6(int(i), self.__data__))))
        return x

    @grids.setter
    def grids(self, val):
        if isinstance(val, Sized) and len(val) == self.dim:
            for i in range(self.dim):
                if isinstance(val[i], Vector):
                    lib.set_gridGriddedField6(int(i), self.__data__, val[i].__data__, c.c_bool(1))  # Ugly use of pointer
                elif isinstance(val[i], ArrayOfString):
                    lib.set_gridGriddedField6(int(i), self.__data__, val[i].__data__, c.c_bool(0))  # Ugly use of pointer
                else:
                    raise TypeError("Expects Vector and/or ArrayOfString input")
        else:
            raise TypeError("Only accepts array-like input of length {}".format(self.dim))

    @property
    def data(self):
        """ The data (Tensor6) """
        return Tensor6(c.c_void_p(lib.dataGriddedField6(self.__data__)))

    @data.setter
    def data(self, val):
        self.data.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        if self.OK:
            lib.printGriddedField6(self.__data__)
        else:
            raise RuntimeError("Class is in bad state")

    def __del__(self):
        if self.__delete__:
            lib.deleteGriddedField6(self.__data__)

    def __repr__(self):
        return "ARTS GriddedField6"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, GriddedField6):
            self.fieldname = other.fieldname
            self.gridnames = other.gridnames
            self.grids = other.grids
            self.data = other.data
        else:
            raise TypeError("Expects GriddedField6")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadGriddedField6(self.__data__, correct_read_arguments(file)):
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
        if not self.OK:
            raise RuntimeError("Class is in bad state")

        if lib.xmlsaveGriddedField6(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, GriddedField6) and \
                self.fieldname == other.fieldname and \
                self.gridnames == other.gridnames and \
                self.grids == other.grids and \
                self.data == other.data:
            return True
        else:
            return False


lib.createGriddedField6.restype = c.c_void_p
lib.createGriddedField6.argtypes = []

lib.deleteGriddedField6.restype = None
lib.deleteGriddedField6.argtypes = [c.c_void_p]

lib.printGriddedField6.restype = None
lib.printGriddedField6.argtypes = [c.c_void_p]

lib.xmlreadGriddedField6.restype = c.c_long
lib.xmlreadGriddedField6.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveGriddedField6.restype = c.c_long
lib.xmlsaveGriddedField6.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.get_dimGriddedField6.restype = c.c_long
lib.get_dimGriddedField6.argtypes = [c.c_void_p]

lib.get_grid_typeIndexGriddedField6.restype = c.c_long
lib.get_grid_typeIndexGriddedField6.argtypes = [c.c_long, c.c_void_p]

lib.get_grid_sizeGriddedField6.restype = c.c_long
lib.get_grid_sizeGriddedField6.argtypes = [c.c_long, c.c_void_p]

lib.get_nameGriddedField6.restype = c.c_char_p
lib.get_nameGriddedField6.argtypes = [c.c_void_p]

lib.set_nameGriddedField6.restype = None
lib.set_nameGriddedField6.argtypes = [c.c_void_p, c.c_char_p]

lib.get_grid_nameGriddedField6.restype = c.c_char_p
lib.get_grid_nameGriddedField6.argtypes = [c.c_long, c.c_void_p]

lib.set_grid_nameGriddedField6.restype = None
lib.set_grid_nameGriddedField6.argtypes = [c.c_long, c.c_void_p, c.c_char_p]

lib.get_numeric_gridGriddedField6.restype = c.c_void_p
lib.get_numeric_gridGriddedField6.argtypes = [c.c_long, c.c_void_p]

lib.get_string_gridGriddedField6.restype = c.c_void_p
lib.get_string_gridGriddedField6.argtypes = [c.c_long, c.c_void_p]

lib.set_gridGriddedField6.restype = c.c_void_p
lib.set_gridGriddedField6.argtypes = [c.c_long, c.c_void_p, c.c_void_p, c.c_bool]

lib.dataGriddedField6.restype = c.c_void_p
lib.dataGriddedField6.argtypes = [c.c_void_p]

lib.checksizeGriddedField6.restype = c.c_bool
lib.checksizeGriddedField6.argtypes = [c.c_void_p]
