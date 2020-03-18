import ctypes as c
from pyarts.workspace.api import arts_api as lib


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
            self.__data__ = c.cast(lib.createRational(), c.c_void_p)
            self.denom = denom
            self.nom = nom
            self.simplify()

    def simplify(self):
        """ Simplify the Rational in place """
        lib.simplifyRational(self.__data__)

    @property
    def nom(self):
        """ Nominator (Index) """
        return lib.getRationalNom(self.__data__)

    @nom.setter
    def nom(self, x):
        lib.setRationalNom(self.__data__, c.c_long(x))

    @property
    def denom(self):
        """ Denominator (Index) """
        return lib.getRationalDenom(self.__data__)

    @denom.setter
    def denom(self, x):
        lib.setRationalDenom(self.__data__, c.c_long(x))

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
        if isinstance(val, int):
            return Rational(self.nom + val * self.denom, self.denom)
        elif isinstance(val, Rational):
            return Rational(self.nom * val.denom + self.denom * val.nom, self.denom*val.denom)
        else:
            raise ValueError("Cannot add {} to Rational".format(val))

    def __sub__(self, val):
        if isinstance(val, int):
            return Rational(self.nom - val * self.denom, self.denom)
        elif isinstance(val, Rational):
            return Rational(self.nom * val.denom - self.denom * val.nom, self.denom*val.denom)
        else:
            raise ValueError("Cannot sub {} to Rational".format(val))

    def __mul__(self, val):
        if isinstance(val, int):
            return Rational(self.nom * val, self.denom)
        elif isinstance(val, Rational):
            return Rational(self.nom*val.nom, self.denom*val.denom)
        else:
            raise ValueError("Cannot mul {} to Rational".format(val))

    def __pow__(self, val):
        if isinstance(val, int):
            return Rational(self.nom ** val, self.denom ** val)
        else:
            raise ValueError("Cannot pow {} to Rational".format(val))

    def __truediv__(self, val):
        if isinstance(val, int):
            return Rational(self.nom, self.denom*val)
        elif isinstance(val, Rational):
            return Rational(self.nom*val.denom, self.denom*val.nom)
        else:
            raise ValueError("Cannot div {} to Rational".format(val))

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Rational):
            self.nom = other.nom
            self.denom = other.denom
        else:
            raise TypeError("Expects Rational")

lib.createRational.restype = c.c_void_p
lib.createRational.argtypes = []

lib.deleteRational.restype = None
lib.deleteRational.argtypes = [c.c_void_p]

lib.printRational.restype = None
lib.printRational.argtypes = [c.c_void_p]

lib.getRationalNom.restype = c.c_long
lib.getRationalNom.argtypes = [c.c_void_p]

lib.getRationalDenom.restype = c.c_long
lib.getRationalDenom.argtypes = [c.c_void_p]

lib.setRationalNom.restype = None
lib.setRationalNom.argtypes = [c.c_void_p, c.c_long]

lib.setRationalDenom.restype = None
lib.setRationalDenom.argtypes = [c.c_void_p, c.c_long]

lib.simplifyRational.restype = None
lib.simplifyRational.argtypes = [c.c_void_p]
