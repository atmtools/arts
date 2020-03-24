import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class InternalGridPos(c.Structure):
    """ ARTS Internal data layout... do not use """
    _fields_ = [("idx", c.c_long),
                ("fd", c.c_double*2),]

    def __init__(self, idx=0, fd=[0, 0]):
        self.idx = int(idx)
        self.fd[0] = float(fd[0])
        self.fd[1] = float(fd[1])

    def print(self):
        """ Print to cout the ARTS representation of the struct """
        lib.printGridPos(c.pointer(self))


class GridPos:
    """ ARTS GridPos data

    Properties:
        idx:
            Index of position below point (Index)

        fd[2]:
            Weight of grids (Numeric, Numeric)
        """
    def __init__(self, idx=0, fd=[0, 0]):
        if isinstance(type, InternalGridPos):
            self.__data__ = idx
        else:
            self.__data__ = InternalGridPos()
            self.idx = idx
            self.fd = fd

    @staticmethod
    def name():
        return "GridPos"

    @property
    def idx(self):
        """ Index of position below point (Index) """
        return self.__data__.idx

    @idx.setter
    def idx(self, val):
        self.__data__.idx = int(val)

    @property
    def fd(self):
        """ Weight of grids (Numeric, Numeric) """
        return self.__data__.fd[0], self.__data__.fd[1]

    @fd.setter
    def fd(self, val):
        self.__data__.fd[0] = float(val[0])
        self.__data__.fd[1] = float(val[1])

    @property
    def data(self):
        """ Full data representation (InternalGridPos) """
        return self.__data__

    def print(self):
        """ Print to cout the ARTS representation of the class """
        self.__data__.print()

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
        if lib.xmlsaveGridPos(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))



exec(array_base(GridPos))


lib.printGridPos.restype = None
lib.printGridPos.argtypes = [c.c_void_p]
