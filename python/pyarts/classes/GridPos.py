import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Index
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class GridPos:
    """ ARTS GridPos data

    Properties:
        idx:
            Index of position below point (Index)

        fd[2]:
            Weight of grids (Numeric, Numeric)
        """
    def __init__(self, idx=0, fd=[0, 0]):
        if isinstance(idx, c.c_void_p):
            self.__delete__ = False
            self.__data__ = idx
        else:
            self.__delete__ = True
            self.__data__ = lib.createGridPos()
            self.idx = idx
            self.fd = fd

    @staticmethod
    def name():
        return "GridPos"

    @property
    def idx(self):
        """ Index of position below point (Index) """
        return Index(c.c_void_p(lib.getidxGridPos(self.__data__)))

    @idx.setter
    def idx(self, val):
        self.idx.set(val)

    @property
    def fd(self):
        """ Weight of grids (Numeric, Numeric) """
        f = lib.getfdGridPos(self.__data__)
        return  c.cast(f, c.POINTER(c.c_double))[0], c.cast(f, c.POINTER(c.c_double))[1]

    @fd.setter
    def fd(self, val):
        f = lib.getfdGridPos(self.__data__)
        c.cast(f, c.POINTER(c.c_double))[0] = float(val[0])
        c.cast(f, c.POINTER(c.c_double))[1] = float(val[1])

    def __repr__(self):
        return "{}: {}".format(self.idx, self.fd)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        if self.OK:
            self.__data__.print()
        else:
            raise ValueError("Class is in bad state")

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, GridPos):
            self.idx = other.idx
            self.fd = other.fd
        else:
            raise TypeError("Expects GridPos")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadGridPos(self.__data__, correct_read_arguments(file)):
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
            raise ValueError("Class is in bad state")

        if lib.xmlsaveGridPos(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, GridPos) and self.idx == other.idx and self.fd == other.fd:
            return True
        else:
            return False

    @property
    def OK(self):
        return bool(self)

    def __bool__(self):
        return self.idx >= 0


exec(array_base(GridPos))


lib.createGridPos.restype = c.c_void_p
lib.createGridPos.argtypes = []

lib.deleteGridPos.restype = None
lib.deleteGridPos.argtypes = [c.c_void_p]

lib.printGridPos.restype = None
lib.printGridPos.argtypes = [c.c_void_p]

lib.xmlreadGridPos.restype = c.c_long
lib.xmlreadGridPos.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveGridPos.restype = c.c_long
lib.xmlsaveGridPos.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getfdGridPos.restype = c.c_void_p
lib.getfdGridPos.argtypes = [c.c_void_p]

lib.getidxGridPos.restype = c.c_void_p
lib.getidxGridPos.argtypes = [c.c_void_p]