import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Rational import Rational
from pyarts.classes.BasicTypes import String

class QuantumNumberType:
    """ ARTS QuantumNumberType type data
    
    Properties:
        name:
            Name of type (String)
    """
    def __init__(self, data):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumNumberType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getQuantumNumberTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setQuantumNumberTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad QuantumNumberType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumNumberType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumNumberType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, QuantumNumberType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, QuantumNumberType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createQuantumNumberType.restype = c.c_void_p
lib.createQuantumNumberType.argtypes = []

lib.deleteQuantumNumberType.restype = None
lib.deleteQuantumNumberType.argtypes = [c.c_void_p]

lib.printQuantumNumberType.restype = None
lib.printQuantumNumberType.argtypes = [c.c_void_p]

lib.getQuantumNumberTypeString.restype = c.c_void_p
lib.getQuantumNumberTypeString.argtypes = [c.c_void_p]

lib.setQuantumNumberTypeString.restype = c.c_int
lib.setQuantumNumberTypeString.argtypes = [c.c_void_p, c.c_char_p]


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
        return lib.sizeQuantumNumbers()

    @property
    def data(self):
        """ The data (list of Rational) """
        x = []
        n = self.size
        for i in range(n):
            x.append(self[i])
        return x

    @data.setter
    def data(self, other):
        n = self.size
        for i in range(n):
            self[i] = other[i]

    def __getitem__(self, ind):
        return Rational(c.c_void_p(lib.getelemQuantumNumbers(self.to_index(ind), self.__data__)))

    def __setitem__(self, ind, val):
        self[ind].set(val)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumNumbers(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumNumbers(self.__data__)

    @property
    def as_string(self):
        return String(c.c_void_p(lib.getQuantumNumbersString(self.__data__)),
                      delete=True)

    def __repr__(self):
        return f"{self.as_string}"

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

        if ind >= 0 and ind < lib.sizeQuantumNumbers():
            return ind
        else:
            raise IndexError("Out of bounds")

    def __eq__(self, other):
        if isinstance(other, QuantumNumbers) and self.data == other.data:
            return True
        else:
            return False

lib.createQuantumNumbers.restype = c.c_void_p
lib.createQuantumNumbers.argtypes = []

lib.deleteQuantumNumbers.restype = None
lib.deleteQuantumNumbers.argtypes = [c.c_void_p]

lib.printQuantumNumbers.restype = None
lib.printQuantumNumbers.argtypes = [c.c_void_p]

lib.getelemQuantumNumbers.restype = c.c_void_p
lib.getelemQuantumNumbers.argtypes = [c.c_long, c.c_void_p]

lib.sizeQuantumNumbers.restype = c.c_long
lib.sizeQuantumNumbers.argtypes = []

lib.string2quantumnumbersindex.restype = c.c_long
lib.string2quantumnumbersindex.argtypes = [c.c_char_p]

lib.getQuantumNumbersString.restype = c.c_void_p
lib.getQuantumNumbersString.argtypes = [c.c_void_p]
