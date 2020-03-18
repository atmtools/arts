import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Rational import Rational


class QuantumNumbers:
    """ ARTS QuantumNumbers data

    Getter and setter is implemented for valid quantum numbers

    Example:
        >>> qns = QuantumNumbers()
        >>> qns["J"] = Rational(3, 2)
        >>> qns.print()
        J 3/2
        >>> qns["SomeFunnyTest"] = Rational(3, 2)
        IndexError: Out of bounds
        >>> qns[0] += Rational(3, 2)
        >>> qns[qns.size-1] = Rational(1)
        >>> qns.print()
        J 3 Hund 1
        >>> qns[qns.size] = Rational(1)
        IndexError: Out of bounds

    Properties:
        size:
            Number of defined quantum numbers (constexpr Index)

        data:
            The data (list of Rational)
    """
    def __init__(self, qns={}):
        if isinstance(qns, c.c_void_p):
            self.__delete__ = False
            self.__data__ = qns
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumNumbers())
            for key in qns:
                self[key] = qns[key]

    @property
    def size(self):
        """ Number of defined quantum numbers (constexpr Index) """
        return lib.getQuantumNumbersMaxNumber()

    @property
    def data(self):
        """ The data (list of Rational) """
        x = []
        n = self.size
        for i in range(n):
            x.append(self[i])
        return x

    @data.setter
    def data(self, val):
        n = self.size
        if isinstance(val, Sized) and len(val) == self.size:
            for i in range(n):
                self[i] = val[i]
        else:
            raise TypeError("Invalid input")

    def __getitem__(self, ind):
        return Rational(c.c_void_p(lib.getQuantumNumbersNumber(self.to_index(ind), self.__data__)))

    def __setitem__(self, ind, val):
        self[ind].set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumNumbers(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumNumbers(self.__data__)

    def __repr__(self):
        return "ARTS QuantumNumbers"

    def __len__(self):
        return self.size

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, QuantumNumbers):
            self.data = other.data
        else:
            raise TypeError("Expects QuantumNumbers")

    @staticmethod
    def to_index(ind):
        """ Returns a quantum number index after checking it

        If the index is bad, raises an IndexError

        Input:
            ind:
                Quantum number (Index or str)

        Output:
            An index
        """
        if isinstance(ind, str):
            ind = ind.encode("ascii")
            ind = lib.string2quantumnumbersindex(ind)
        else:
            ind = int(ind)

        if ind >= 0 and ind < lib.getQuantumNumbersMaxNumber():
            return ind
        else:
            raise IndexError("Out of bounds")


lib.createQuantumNumbers.restype = c.c_void_p
lib.createQuantumNumbers.argtypes = []

lib.deleteQuantumNumbers.restype = None
lib.deleteQuantumNumbers.argtypes = [c.c_void_p]

lib.printQuantumNumbers.restype = None
lib.printQuantumNumbers.argtypes = [c.c_void_p]

lib.getQuantumNumbersMaxNumber.restype = c.c_long
lib.getQuantumNumbersMaxNumber.argtypes = []

lib.getQuantumNumbersNumber.restype = c.c_void_p
lib.getQuantumNumbersNumber.argtypes = [c.c_long, c.c_void_p]

lib.string2quantumnumbersindex.restype = c.c_long
lib.string2quantumnumbersindex.argtypes = [c.c_char_p]
