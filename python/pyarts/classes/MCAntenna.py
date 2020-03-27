import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Numeric
from pyarts.classes.Vector import Vector
from pyarts.classes.Matrix import Matrix
from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class MCAntenna:
    """ ARTS MCAntenna data

    Properties:
        type:
            (Index)

        sigma_aa:
            (Numeric)

        sigma_za:
            (Numeric)

        aa_grid:
            (Vector)

        za_grid:
            (Vector)

        g_lookup:
            (Matrix)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createMCAntenna())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def type(self):
        """ (Index) """
        return lib.getTypeMCAntenna(self.__data__)

    @type.setter
    def type(self, val):
        if lib.setTypeMCAntenna(self.__data__, int(val)):
            raise ValueError("Bad input")

    @property
    def sigma_aa(self):
        """ (Numeric) """
        return Numeric(c.c_void_p(lib.getsaaMCAntenna(self.__data__)))

    @sigma_aa.setter
    def sigma_aa(self, val):
        self.sigma_aa.set(val)

    @property
    def sigma_za(self):
        """ (Numeric) """
        return Numeric(c.c_void_p(lib.getszaMCAntenna(self.__data__)))

    @sigma_za.setter
    def sigma_za(self, val):
        self.sigma_za.set(val)

    @property
    def aa_grid(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getaagMCAntenna(self.__data__)))

    @aa_grid.setter
    def aa_grid(self, val):
        self.aa_grid.set(val)

    @property
    def za_grid(self):
        """ (Vector) """
        return Vector(c.c_void_p(lib.getzagMCAntenna(self.__data__)))

    @za_grid.setter
    def za_grid(self, val):
        self.za_grid.set(val)

    @property
    def g_lookup(self):
        """ (Matrix) """
        return Matrix(c.c_void_p(lib.getGMCAntenna(self.__data__)))

    @g_lookup.setter
    def g_lookup(self, val):
        self.g_lookup.set(val)

    @staticmethod
    def name():
        return "MCAntenna"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printMCAntenna(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteMCAntenna(self.__data__)

    def __repr__(self):
        return "ARTS MCAntenna"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, MCAntenna):
              raise RuntimeWarning("Cannot set MCAntenna, remains constant")
        else:
            raise TypeError("Expects MCAntenna")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadMCAntenna(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveMCAntenna(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createMCAntenna.restype = c.c_void_p
lib.createMCAntenna.argtypes = []

lib.deleteMCAntenna.restype = None
lib.deleteMCAntenna.argtypes = [c.c_void_p]

lib.printMCAntenna.restype = None
lib.printMCAntenna.argtypes = [c.c_void_p]

lib.xmlreadMCAntenna.restype = c.c_long
lib.xmlreadMCAntenna.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveMCAntenna.restype = c.c_long
lib.xmlsaveMCAntenna.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getTypeMCAntenna.restype = c.c_long
lib.getTypeMCAntenna.argtypes = [c.c_void_p]

lib.setTypeMCAntenna.restype = c.c_long
lib.setTypeMCAntenna.argtypes = [c.c_void_p, c.c_long]

lib.getsaaMCAntenna.restype = c.c_void_p
lib.getsaaMCAntenna.argtypes = [c.c_void_p]

lib.getszaMCAntenna.restype = c.c_void_p
lib.getszaMCAntenna.argtypes = [c.c_void_p]

lib.getaagMCAntenna.restype = c.c_void_p
lib.getaagMCAntenna.argtypes = [c.c_void_p]

lib.getzagMCAntenna.restype = c.c_void_p
lib.getzagMCAntenna.argtypes = [c.c_void_p]

lib.getGMCAntenna.restype = c.c_void_p
lib.getGMCAntenna.argtypes = [c.c_void_p]
