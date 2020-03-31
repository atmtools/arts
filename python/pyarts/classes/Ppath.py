import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Index, Numeric, String
from pyarts.classes.Matrix import Matrix
from pyarts.classes.Vector import Vector
from pyarts.classes.GridPos import ArrayOfGridPos
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class Ppath:
    """ ARTS Ppath data

    Properties:
        dim:
            Atmospheric dimensionality (Index)

        np:
            Number of points describing the ppath (Index)

        constant:
            The propagation path constant (only used for 1D) (Numeric)

        background:
            Radiative background (String)

        start_pos:
            Start position (Vector)

        start_los:
            Start line-of-sight (Vector)

        start_lstep:
            Length between sensor and atmospheric boundary (Numeric)

        pos:
            The distance between start pos and the last position in pos (Matrix)

        los:
            Line-of-sight at each ppath point (Matrix)

        r:
            Radius of each ppath point (Vector)

        lstep:
            The length between ppath points (Vector)

        end_pos:
            End position (Vector)

        end_los:
            End line-of-sight (Vector)

        end_lstep:
            The distance between end pos and the first position in pos (Numeric)

        nreal:
            The real part of the refractive index at each path position (Vector)

        ngroup:
            The group index of refraction (Vector)

        gp_p:
            Index position with respect to the pressure grid (ArrayOfGridPos)

        gp_lat:
            Index position with respect to the latitude grid (ArrayOfGridPos)

        gp_lon:
            Index position with respect to the longitude grid (ArrayOfGridPos)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createPpath())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "Ppath"

    @property
    def dim(self):
        """ Atmospheric dimensionality (Index) """
        return Index(c.c_void_p(lib.getdimPpath(self.__data__)))

    @dim.setter
    def dim(self, val):
        self.dim.set(val)

    @property
    def np(self):
        """ Number of points describing the ppath (Index) """
        return Index(c.c_void_p(lib.getnpPpath(self.__data__)))

    @np.setter
    def np(self, val):
        self.np.set(val)

    @property
    def constant(self):
        """ The propagation path constant (only used for 1D) (Numeric) """
        return Numeric(c.c_void_p(lib.getconstantPpath(self.__data__)))

    @constant.setter
    def constant(self, val):
        self.constant.set(val)

    @property
    def background(self):
        """ Radiative background (String) """
        return String(c.c_void_p(lib.getbackgroundPpath(self.__data__)))

    @background.setter
    def background(self, val):
        self.background.set(val)

    @property
    def start_pos(self):
        """ Start position (Vector) """
        return Vector(c.c_void_p(lib.getstart_posPpath(self.__data__)))

    @start_pos.setter
    def start_pos(self, val):
        self.start_pos.set(val)

    @property
    def start_los(self):
        """ Start line-of-sight (Vector) """
        return Vector(c.c_void_p(lib.getstart_losPpath(self.__data__)))

    @start_los.setter
    def start_los(self, val):
        self.start_los.set(val)

    @property
    def start_lstep(self):
        """ Length between sensor and atmospheric boundary (Numeric) """
        return Numeric(c.c_void_p(lib.getstart_lstepPpath(self.__data__)))

    @start_lstep.setter
    def start_lstep(self, val):
        self.start_lstep.set(val)

    @property
    def pos(self):
        """ The distance between start pos and the last position in pos (Matrix) """
        return Matrix(c.c_void_p(lib.getposPpath(self.__data__)))

    @pos.setter
    def pos(self, val):
        self.pos.set(val)

    @property
    def los(self):
        """ Line-of-sight at each ppath point (Matrix) """
        return Matrix(c.c_void_p(lib.getlosPpath(self.__data__)))

    @los.setter
    def los(self, val):
        self.los.set(val)

    @property
    def r(self):
        """ Radius of each ppath point (Vector) """
        return Vector(c.c_void_p(lib.getrPpath(self.__data__)))

    @r.setter
    def r(self, val):
        self.r.set(val)

    @property
    def lstep(self):
        """ The length between ppath points (Vector) """
        return Vector(c.c_void_p(lib.getlstepPpath(self.__data__)))

    @lstep.setter
    def lstep(self, val):
        self.lstep.set(val)

    @property
    def end_pos(self):
        """ End position (Vector) """
        return Vector(c.c_void_p(lib.getend_posPpath(self.__data__)))

    @end_pos.setter
    def end_pos(self, val):
        self.end_pos.set(val)

    @property
    def end_los(self):
        """ End line-of-sight (Vector) """
        return Vector(c.c_void_p(lib.getend_losPpath(self.__data__)))

    @end_los.setter
    def end_los(self, val):
        self.end_los.set(val)

    @property
    def end_lstep(self):
        """ The distance between end pos and the first position in pos (Numeric) """
        return Numeric(c.c_void_p(lib.getend_lstepPpath(self.__data__)))

    @end_lstep.setter
    def end_lstep(self, val):
        self.end_lstep.set(val)

    @property
    def nreal(self):
        """ The real part of the refractive index at each path position (Vector) """
        return Vector(c.c_void_p(lib.getnrealPpath(self.__data__)))

    @nreal.setter
    def nreal(self, val):
        self.nreal.set(val)

    @property
    def ngroup(self):
        """ The group index of refraction (Vector) """
        return Vector(c.c_void_p(lib.getngroupPpath(self.__data__)))

    @ngroup.setter
    def ngroup(self, val):
        self.ngroup.set(val)

    @property
    def gp_p(self):
        """ Index position with respect to the pressure grid (ArrayOfGridPos) """
        return ArrayOfGridPos(c.c_void_p(lib.getgp_pPpath(self.__data__)))

    @gp_p.setter
    def gp_p(self, val):
        self.gp_p.set(val)

    @property
    def gp_lat(self):
        """ Index position with respect to the latitude grid (ArrayOfGridPos) """
        return ArrayOfGridPos(c.c_void_p(lib.getgp_latPpath(self.__data__)))

    @gp_lat.setter
    def gp_lat(self, val):
        self.gp_lat.set(val)

    @property
    def gp_lon(self):
        """ Index position with respect to the longitude grid (ArrayOfGridPos) """
        return ArrayOfGridPos(c.c_void_p(lib.getgp_lonPpath(self.__data__)))

    @gp_lon.setter
    def gp_lon(self, val):
        self.gp_lon.set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printPpath(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deletePpath(self.__data__)

    def __repr__(self):
        return "ARTS Ppath"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Ppath):
            self.dim = other.dim
            self.np = other.np
            self.constant = other.constant
            self.background = other.background
            self.start_pos = other.start_pos
            self.start_los = other.start_los
            self.start_lstep = other.start_lstep
            self.pos = other.pos
            self.los = other.los
            self.r = other.r
            self.lstep = other.lstep
            self.end_pos = other.end_pos
            self.end_los = other.end_los
            self.end_lstep = other.end_lstep
            self.nreal = other.nreal
            self.ngroup = other.ngroup
            self.gp_p = other.gp_p
            self.gp_lat = other.gp_lat
            self.gp_lon = other.gp_lon
        else:
            raise TypeError("Expects Ppath")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadPpath(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsavePpath(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Ppath) and \
                self.dim == other.dim and \
                self.np == other.np and \
                self.constant == other.constant and \
                self.background == other.background and \
                self.start_pos == other.start_pos and \
                self.start_los == other.start_los and \
                self.start_lstep == other.start_lstep and \
                self.pos == other.pos and \
                self.los == other.los and \
                self.r == other.r and \
                self.lstep == other.lstep and \
                self.end_pos == other.end_pos and \
                self.end_los == other.end_los and \
                self.end_lstep == other.end_lstep and \
                self.nreal == other.nreal and \
                self.ngroup == other.ngroup and \
                self.gp_p == other.gp_p and \
                self.gp_lat == other.gp_lat and \
                self.gp_lon == other.gp_lon:
            return True
        else:
            return False

    def __bool__(self):
        return self.np > 0


exec(array_base(Ppath))


lib.createPpath.restype = c.c_void_p
lib.createPpath.argtypes = []

lib.deletePpath.restype = None
lib.deletePpath.argtypes = [c.c_void_p]

lib.printPpath.restype = None
lib.printPpath.argtypes = [c.c_void_p]

lib.xmlreadPpath.restype = c.c_long
lib.xmlreadPpath.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsavePpath.restype = c.c_long
lib.xmlsavePpath.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getdimPpath.restype = c.c_void_p
lib.getdimPpath.argtypes = [c.c_void_p]

lib.getnpPpath.restype = c.c_void_p
lib.getnpPpath.argtypes = [c.c_void_p]

lib.getconstantPpath.restype = c.c_void_p
lib.getconstantPpath.argtypes = [c.c_void_p]

lib.getbackgroundPpath.restype = c.c_void_p
lib.getbackgroundPpath.argtypes = [c.c_void_p]

lib.getstart_posPpath.restype = c.c_void_p
lib.getstart_posPpath.argtypes = [c.c_void_p]

lib.getstart_losPpath.restype = c.c_void_p
lib.getstart_losPpath.argtypes = [c.c_void_p]

lib.getstart_lstepPpath.restype = c.c_void_p
lib.getstart_lstepPpath.argtypes = [c.c_void_p]

lib.getposPpath.restype = c.c_void_p
lib.getposPpath.argtypes = [c.c_void_p]

lib.getlosPpath.restype = c.c_void_p
lib.getlosPpath.argtypes = [c.c_void_p]

lib.getrPpath.restype = c.c_void_p
lib.getrPpath.argtypes = [c.c_void_p]

lib.getlstepPpath.restype = c.c_void_p
lib.getlstepPpath.argtypes = [c.c_void_p]

lib.getend_posPpath.restype = c.c_void_p
lib.getend_posPpath.argtypes = [c.c_void_p]

lib.getend_losPpath.restype = c.c_void_p
lib.getend_losPpath.argtypes = [c.c_void_p]

lib.getend_lstepPpath.restype = c.c_void_p
lib.getend_lstepPpath.argtypes = [c.c_void_p]

lib.getnrealPpath.restype = c.c_void_p
lib.getnrealPpath.argtypes = [c.c_void_p]

lib.getngroupPpath.restype = c.c_void_p
lib.getngroupPpath.argtypes = [c.c_void_p]

lib.getgp_pPpath.restype = c.c_void_p
lib.getgp_pPpath.argtypes = [c.c_void_p]

lib.getgp_latPpath.restype = c.c_void_p
lib.getgp_latPpath.argtypes = [c.c_void_p]

lib.getgp_lonPpath.restype = c.c_void_p
lib.getgp_lonPpath.argtypes = [c.c_void_p]
