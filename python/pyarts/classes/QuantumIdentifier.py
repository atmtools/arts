import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.QuantumNumbers import QuantumNumbers
from pyarts.classes.SpeciesIsotopeRecord import SpeciesIsotopeRecord
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base
from pyarts.classes.BasicTypes import Index, String


class QuantumIdentifierType:
    """ ARTS QuantumIdentifier type data
    
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
            self.__data__ = c.c_void_p(lib.createQuantumIdentifierType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getQuantumIdentifierTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setQuantumIdentifierTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad QuantumIdentifierType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumIdentifierType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumIdentifierType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, QuantumIdentifierType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, QuantumIdentifierType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createQuantumIdentifierType.restype = c.c_void_p
lib.createQuantumIdentifierType.argtypes = []

lib.deleteQuantumIdentifierType.restype = None
lib.deleteQuantumIdentifierType.argtypes = [c.c_void_p]

lib.printQuantumIdentifierType.restype = None
lib.printQuantumIdentifierType.argtypes = [c.c_void_p]

lib.getQuantumIdentifierTypeString.restype = c.c_void_p
lib.getQuantumIdentifierTypeString.argtypes = [c.c_void_p]

lib.setQuantumIdentifierTypeString.restype = c.c_int
lib.setQuantumIdentifierTypeString.argtypes = [c.c_void_p, c.c_char_p]


class QuantumIdentifier:
    """ ARTS QuantumIdentifier data

    Properties:
        type:
            The type of identifier (QuantumIdentifierType)
            
        spec_ind:
            The index of the isotopologue (Index)
        
        upp:
            The upper quantum numbers (QuantumNumbers)
            
        low:
            The lower quantum numbers (QuantumNumbers)
            
        lvl:
            The level quantum numbers (QuantumNumbers)
            
        isot:
            The ARTS isotopologue (const SpeciesIsotopeRecord)
    """
    def __init__(self, type=None):
        if isinstance(type, c.c_void_p):
            self.__delete__ = False
            self.__data__ = type
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumIdentifier())
            if type is not None:
                if lib.fromstringQuantumIdentifier(self.__data__, str(type).encode('utf-8')):
                    raise RuntimeError(f"Bad QuantumIdentifier: {type}")
    
    @property
    def type(self):
        return QuantumIdentifierType(c.c_void_p(lib.gettypeQuantumIdentifier(self.__data__)))
    
    @type.setter
    def type(self, x):
        self.type.set(x)
    
    @property
    def spec_ind(self):
        return Index(c.c_void_p(lib.getspec_indQuantumIdentifier(self.__data__)))
    
    @spec_ind.setter
    def spec_ind(self, x):
        self.spec_ind.set(x)
    
    @property
    def upp(self):
        if self.type != "Transition":
            raise RuntimeError(f"Bad state, {self.type} should be Transition")
        return QuantumNumbers(c.c_void_p(lib.getuppQuantumIdentifier(self.__data__)))
    
    @upp.setter
    def upp(self, x):
        self.upp.set(x)
    
    @property
    def lvl(self):
        if self.type != "EnergyLevel":
            raise RuntimeError(f"Bad state, {self.type} should be EnergyLevel")
        return QuantumNumbers(c.c_void_p(lib.getuppQuantumIdentifier(self.__data__)))
    
    @lvl.setter
    def lvl(self, x):
        self.lvl.set(x)
    
    @property
    def low(self):
        if self.type != "Transition":
            raise RuntimeError(f"Bad state, {self.type} should be Transition")
        return QuantumNumbers(c.c_void_p(lib.getlowQuantumIdentifier(self.__data__)))
    
    @low.setter
    def low(self, x):
        self.low.set(x)
    
    @property
    def isot(self):
        return SpeciesIsotopeRecord.from_index(self.spec_ind)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumIdentifier(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumIdentifier(self.__data__)
    
    def __repr__(self):
        if self.type == "Transition":
            return f"{self.isot} TR UP {self.upp} LO {self.low}"
        elif self.type == "EnergyLevel":
            return f"{self.isot} EN {self.lvl}"
        elif self.type == "All":
            return f"{self.isot} ALL"
        else:
            return "NONE"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, QuantumIdentifier):
            self.type = other.type
            self.spec_ind = other.spec_ind
            if type == "EnergyLevel":
                self.lvl = other.lvl
            elif type == "Transtition":
                self.upp = other.upp
                self.low = other.low
        else:
            raise TypeError("Expects QuantumIdentifier")

    @staticmethod
    def name():
        return "QuantumIdentifier"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadQuantumIdentifier(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveQuantumIdentifier(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, QuantumIdentifier):
            if self.type == "Transition":
                return self.type == other.type and self.isot == other.isot and other.upp == self.upp and other.low == self.low
            elif self.type == "EnergyLevel":
                return self.type == other.type and self.isot == other.isot and other.lvl == self.lvl
            elif self.type == "All":
                return self.type == other.type and self.isot == other.isot
            else:
                return self.type == other.type
        else:
            return False


# Generate ArrayOfQuantumIdentifier
exec(array_base(QuantumIdentifier))


lib.createQuantumIdentifier.restype = c.c_void_p
lib.createQuantumIdentifier.argtypes = []

lib.deleteQuantumIdentifier.restype = None
lib.deleteQuantumIdentifier.argtypes = [c.c_void_p]

lib.printQuantumIdentifier.restype = None
lib.printQuantumIdentifier.argtypes = [c.c_void_p]

lib.xmlreadQuantumIdentifier.restype = c.c_long
lib.xmlreadQuantumIdentifier.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveQuantumIdentifier.restype = c.c_long
lib.xmlsaveQuantumIdentifier.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.gettypeQuantumIdentifier.restype = c.c_void_p
lib.gettypeQuantumIdentifier.argtypes = [c.c_void_p]

lib.getspec_indQuantumIdentifier.restype = c.c_void_p
lib.getspec_indQuantumIdentifier.argtypes = [c.c_void_p]

lib.getuppQuantumIdentifier.restype = c.c_void_p
lib.getuppQuantumIdentifier.argtypes = [c.c_void_p]

lib.getlowQuantumIdentifier.restype = c.c_void_p
lib.getlowQuantumIdentifier.argtypes = [c.c_void_p]

lib.fromstringQuantumIdentifier.restype = c.c_long
lib.fromstringQuantumIdentifier.argtypes = [c.c_void_p, c.c_char_p]
