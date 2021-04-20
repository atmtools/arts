import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import Numeric, String
from pyarts.classes.Matrix import Matrix
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class Star:
    """ARTS Star data

    Properties:
        description:
            Star description (String)

        spectrum:
            Star spectrum, monochrmatic radiance spectrum at the surface of the
            star (Matrix)

        radius:
            Star radius (Numeric)

        distance:
            Star distance from center of planet to center of star (Numeric)

        latitude:
            Latitude of the star in the sky of the planet (Numeric)

        longitude:
            Longitude of the star in the sky of the planet (Numeric)
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createStar())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "Star"

    @property
    def description(self):
        """Star description (String)"""
        return String(c.c_void_p(lib.getdescriptionStar(self.__data__)))

    @description.setter
    def description(self, val):
        self.description.set(val)

    @property
    def spectrum(self):
        """
        Star spectrum, monochrmatic radiance spectrum at the surface of the
        star (Matrix)
        """
        return Matrix(c.c_void_p(lib.getspectrumStar(self.__data__)))

    @spectrum.setter
    def spectrum(self, val):
        self.spectrum.set(val)

    @property
    def radius(self):
        """Star radius (Numeric)"""
        return Numeric(c.c_void_p(lib.getradiusStar(self.__data__)))

    @radius.setter
    def radius(self, val):
        self.radius.set(val)

    @property
    def distance(self):
        """Star distance from center of planet to center of star (Numeric)"""
        return Numeric(c.c_void_p(lib.getdistanceStar(self.__data__)))

    @distance.setter
    def distance(self, val):
        self.distance.set(val)

    @property
    def latitude(self):
        """Latitude of the star in the sky of the planet (Numeric)"""
        return Numeric(c.c_void_p(lib.getlatitudeStar(self.__data__)))

    @latitude.setter
    def latitude(self, val):
        self.latitude.set(val)

    @property
    def longitude(self):
        """Longitude of the star in the sky of the planet (Numeric)"""
        return Numeric(c.c_void_p(lib.getlongitudeStar(self.__data__)))

    @longitude.setter
    def longitude(self, val):
        self.longitude.set(val)

    def print(self):
        """Print to cout the ARTS representation of the class"""
        lib.printStar(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteStar(self.__data__)

    def __repr__(self):
        return "ARTS Star"

    def set(self, other):
        """Sets this class according to another python instance of itself"""
        if isinstance(other, Star):
            self.description = other.description
            self.spectrum = other.spectrum
            self.radius = other.radius
            self.distance = other.distance
            self.latitude = other.latitude
            self.longitude = other.longitude
        else:
            raise TypeError("Expects Star")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadStar(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveStar(self.__data__,
                           *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Star) and \
                self.description == other.description and \
                self.spectrum == other.spectrum and \
                self.radius == other.radius and \
                self.distance == other.distance and \
                self.latitude == other.latitude and \
                self.longitude == other.longitude:
            return True
        else:
            return False

    def __bool__(self):
        return not self.radius <= 0


exec(array_base(Star))

lib.createStar.restype = c.c_void_p
lib.createStar.argtypes = []

lib.deleteStar.restype = None
lib.deleteStar.argtypes = [c.c_void_p]

lib.printStar.restype = None
lib.printStar.argtypes = [c.c_void_p]

lib.xmlreadStar.restype = c.c_long
lib.xmlreadStar.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveStar.restype = c.c_long
lib.xmlsaveStar.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getdescriptionStar.restype = c.c_void_p
lib.getdescriptionStar.argtypes = [c.c_void_p]

lib.getspectrumStar.restype = c.c_void_p
lib.getspectrumStar.argtypes = [c.c_void_p]

lib.getradiusStar.restype = c.c_void_p
lib.getradiusStar.argtypes = [c.c_void_p]

lib.getdistanceStar.restype = c.c_void_p
lib.getdistanceStar.argtypes = [c.c_void_p]

lib.getlatitudeStar.restype = c.c_void_p
lib.getlatitudeStar.argtypes = [c.c_void_p]

lib.getlongitudeStar.restype = c.c_void_p
lib.getlongitudeStar.argtypes = [c.c_void_p]
