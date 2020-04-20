import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments

from fractions import Fraction

class Rational:
    """ ARTS Rational data

    Basic math is implemented

    Properties:
        nom:
            Nominator (Index)

        denom:
            Denominator (Index)
        """
    def __init__(self, nom=0, denom=1):
        if isinstance(nom, c.c_void_p):
            self.__delete__ = False
            self.__data__ = nom
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createRational())
            self.denom = denom
            self.nom = nom
            self.simplify()

    def simplify(self):
        """ Simplify the Rational in place """
        lib.simplifyRational(self.__data__)

    @property
    def nom(self):
        """ Nominator (Index) """
        return lib.getNomRational(self.__data__)

    @nom.setter
    def nom(self, x):
        lib.setNomRational(self.__data__, int(x))

    @property
    def denom(self):
        """ Denominator (Index) """
        return lib.getDenomRational(self.__data__)

    @denom.setter
    def denom(self, x):
        lib.setDenomRational(self.__data__, int(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printRational(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteRational(self.__data__)

    def __repr__(self):
        return "ARTS Rational {}/{}".format(self.nom, self.denom)

    def __abs__(self):
        return Rational(abs(self.nom), abs(self.denom))

    def __bool__(self):
        return self.nom != 0 and self.denom != 0

    def __add__(self, val):
        if isinstance(val, Rational):
            return Rational(self.nom * val.denom + self.denom * val.nom, self.denom*val.denom)
        else:
            r = Rational()
            r.set(val)
            return self + r

    def __sub__(self, val):
        if isinstance(val, Rational):
            return Rational(self.nom * val.denom - self.denom * val.nom, self.denom*val.denom)
        else:
            r = Rational()
            r.set(val)
            return self - r

    def __mul__(self, val):
        if isinstance(val, Rational):
            return Rational(self.nom*val.nom, self.denom*val.denom)
        else:
            r = Rational()
            r.set(val)
            return self * r

    def __pow__(self, val):
        if isinstance(val, int):
            return Rational(self.nom ** val, self.denom ** val)
        else:
            return float(self)**val

    def __truediv__(self, val):
        if isinstance(val, Rational):
            return Rational(self.nom*val.denom, self.denom*val.nom)
        else:
            r = Rational()
            r.set(val)
            return self / r

    def __float__(self):
        return self.nom / self.denom

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Rational):
            self.nom = other.nom
            self.denom = other.denom
        else:
            x = Fraction(other)
            self.nom = x.numerator
            self.denom = x.denominator

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadRational(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveRational(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, Rational) and other.nom * self.denom == self.nom * other.denom:
            return True
        else:
            return False


lib.createRational.restype = c.c_void_p
lib.createRational.argtypes = []

lib.deleteRational.restype = None
lib.deleteRational.argtypes = [c.c_void_p]

lib.printRational.restype = None
lib.printRational.argtypes = [c.c_void_p]

lib.xmlreadRational.restype = c.c_long
lib.xmlreadRational.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveRational.restype = c.c_long
lib.xmlsaveRational.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getNomRational.restype = c.c_long
lib.getNomRational.argtypes = [c.c_void_p]

lib.getDenomRational.restype = c.c_long
lib.getDenomRational.argtypes = [c.c_void_p]

lib.setNomRational.restype = None
lib.setNomRational.argtypes = [c.c_void_p, c.c_long]

lib.setDenomRational.restype = None
lib.setDenomRational.argtypes = [c.c_void_p, c.c_long]

lib.simplifyRational.restype = None
lib.simplifyRational.argtypes = [c.c_void_p]
