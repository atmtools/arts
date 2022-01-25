import ctypes as c
from pyarts.workspace.api import arts_api as lib
from pyarts.classes.ArrayBase import array_base

from pyarts.classes.io import correct_save_arguments, correct_read_arguments

from pyarts.classes.macros import BasicInterfaceCAPI, EnumMacroInterfaceCAP, VoidStructGetterCAPI
from pyarts.classes.BasicTypes import String, Index
from pyarts.classes.Rational import Rational
from pyarts.classes.SpeciesIsotopeRecord import SpeciesIsotopeRecord

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
        return String(c.c_void_p(lib.getQuantumNumberType(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setQuantumNumberType(self.__data__, str(x).encode('utf-8')):
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


BasicInterfaceCAPI(lib, "QuantumNumberType")
EnumMacroInterfaceCAP(lib, "QuantumNumberType")


class QuantumNumberValue:
    """ ARTS QuantumNumberValue type data
    
    Use set("QN VAL1 VAL2") to change values
    
    Properties:
        type:
            A const QuantumNumberType
        upp:
            A const Rational
        low:
            A const Rational
        str_upp:
            A const Rational
        str_low:
            A const Rational
    """
    def __init__(self, value, delete=False):
        if isinstance(value, c.c_void_p):
            self.__delete__ = delete
            self.__data__ = value
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumNumberValue())
            self.set(value)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumNumberValue(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumNumberValue(self.__data__)
        
    @property
    def type(self):
        return QuantumNumberType(c.c_void_p(lib.gettypeQuantumNumberValue(self.__data__)))
    
    @property
    def upp(self):
        return Rational(c.c_void_p(lib.getuppQuantumNumberValue(self.__data__, False)))
    
    @property
    def str_upp(self):
        return String(c.c_void_p(lib.getuppQuantumNumberValue(self.__data__, True), delete=True))
    
    @property
    def low(self):
        return Rational(c.c_void_p(lib.getlowQuantumNumberValue(self.__data__, False)))
    
    @property
    def str_low(self):
        return String(c.c_void_p(lib.getlowQuantumNumberValue(self.__data__, True), delete=True))
    
    def __repr__(self):
        return f"{self.type} {self.str_upp} {self.str_low}"
    
    def set(self, other):
        if isinstance(other, str):
            if lib.setqnQuantumNumberValue(self.__data__, other.encode('utf-8')):
                raise RuntimeError(f"Bad QuantumNumberValue: '{other}'")
        else:
            self.set(str(other))
    
    def __eq__(self, other):
        return isinstance(other, QuantumNumberValue) and other.type == self.type and other.upp_str == self.upp_str and other.low_str == self.low_str
    
BasicInterfaceCAPI(lib, "QuantumNumberValue")
VoidStructGetterCAPI(lib, "QuantumNumberValue", "type")

lib.setqnQuantumNumberValue.restype = c.c_int
lib.setqnQuantumNumberValue.argtypes = [c.c_void_p, c.c_char_p]

lib.getuppQuantumNumberValue.restype = c.c_void_p
lib.getuppQuantumNumberValue.argtypes = [c.c_void_p, c.c_bool]

lib.getlowQuantumNumberValue.restype = c.c_void_p
lib.getlowQuantumNumberValue.argtypes = [c.c_void_p, c.c_bool]

class QuantumNumberValueList:
    def __init__(self, value):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumNumberValueList())
            self.set(value)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumNumberValueList(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumNumberValueList(self.__data__)
    
    def __getitem__(self, key):
        x = lib.getQuantumNumberValueList(self.__data__, str(key).encode("utf-8"))
        if x:
            return QuantumNumberValue(c.c_void_p(x), delete=True)
        else:
            raise RuntimeError(f"Cannot access key: {key}")
    
    def __setitem__(self, key, val):
        if lib.setQuantumNumberValueList(self.__data__, str(key).encode("utf-8"), str(val).encode("utf-8")):
            raise RuntimeError(f"Cannot set key '{key}' to value '{val}'")
    
    def __repr__(self):
        return f"{String(c.c_void_p(lib.getQuantumNumberValueListString(self.__data__)), delete=True)}"
    
    def set(self, val):
        if lib.setQuantumNumberValueListString(self.__data__, str(val).encode('utf-8')):
            raise RuntimeError(f"Bad QuantumNumberValueList: {val}")
    
    def __eq__(self, o):
        if isinstance(o, QuantumNumberValueList):
            return str(self) == str(o)
        else:
            return self == QuantumNumberValueList(o)

BasicInterfaceCAPI(lib, "QuantumNumberValueList")

lib.getQuantumNumberValueList.restype = c.c_void_p
lib.getQuantumNumberValueList.argtypes = [c.c_void_p, c.c_char_p]

lib.setQuantumNumberValueList.restype = c.c_int
lib.setQuantumNumberValueList.argtypes = [c.c_void_p, c.c_char_p, c.c_char_p]

lib.getQuantumNumberValueListString.restype = c.c_void_p
lib.getQuantumNumberValueListString.argtypes = [c.c_void_p]

lib.setQuantumNumberValueListString.restype = c.c_int
lib.setQuantumNumberValueListString.argtypes = [c.c_void_p, c.c_char_p]

class QuantumNumberLocalState:
    def __init__(self, value=""):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumNumberLocalState())
            self.val = value

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumNumberLocalState(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumNumberLocalState(self.__data__)
    
    @property
    def val(self):
        return QuantumNumberValueList(c.c_void_p(lib.getvalQuantumNumberLocalState(self.__data__)))
    
    @val.setter
    def val(self, val):
        self.val.set(val)
    
    def __repr__(self):
        return f"LocalState: {self.val}"

    def __eq__(self, o):
        if isinstance(o, QuantumNumberLocalState):
            return self.val == o.val
        else:
            return self == QuantumNumberLocalState(o)

    def set(self, o):
        if isinstance(o, QuantumNumberLocalState):
            self.val = o.val
        else:
            self.set(QuantumNumberLocalState(o))

BasicInterfaceCAPI(lib, "QuantumNumberLocalState")
VoidStructGetterCAPI(lib, "QuantumNumberLocalState", "val")

class QuantumIdentifier:
    def __init__(self, value=None):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createQuantumIdentifier())
            if value is not None and lib.fromstringQuantumIdentifier(self.__data__, str(value).encode("utf-8")):
                raise RuntimeError(f"Invalid QuantumIdentifier: {value}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printQuantumIdentifier(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteQuantumIdentifier(self.__data__)
    
    @property
    def isot(self):
        return SpeciesIsotopeRecord.from_index(Index(c.c_void_p(lib.getisotopologue_indexQuantumIdentifier(self.__data__))))

    @isot.setter
    def isot(self, o):
        if isinstance(o, int) or isinstance(o, Index):
            Index(c.c_void_p(lib.getisotopologue_indexQuantumIdentifier(self.__data__))).set(SpeciesIsotopeRecord.from_index(o).index)
        else:
            Index(c.c_void_p(lib.getisotopologue_indexQuantumIdentifier(self.__data__))).set(SpeciesIsotopeRecord(str(o)).index)

    @property
    def val(self):
        return QuantumNumberValueList(c.c_void_p(lib.getvalQuantumIdentifier(self.__data__)))
    
    @val.setter
    def val(self, val):
        self.val.set(val)
    
    def __repr__(self):
        return f"{self.isot} {self.val}"

    def __eq__(self, o):
        if isinstance(o, QuantumIdentifier):
            return Index(c.c_void_p(lib.getisotopologue_indexQuantumIdentifier(self.__data__))) == Index(c.c_void_p(lib.getisotopologue_indexQuantumIdentifier(o.__data__))) and self.val == o.val
        else:
            return self == QuantumIdentifier(o)

    def set(self, o):
        if isinstance(o, QuantumIdentifier):
            Index(c.c_void_p(lib.getisotopologue_indexQuantumIdentifier(self.__data__))).set(Index(c.c_void_p(lib.getisotopologue_indexQuantumIdentifier(o.__data__))))
            self.val = o.val
        else:
            self.set(QuantumIdentifier(o))
            
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

BasicInterfaceCAPI(lib, "QuantumIdentifier")
VoidStructGetterCAPI(lib, "QuantumIdentifier", "isotopologue_index")
VoidStructGetterCAPI(lib, "QuantumIdentifier", "val")

lib.fromstringQuantumIdentifier.restype = c.c_long
lib.fromstringQuantumIdentifier.argtypes = [c.c_void_p, c.c_char_p]

# Generate ArrayOfQuantumIdentifier
exec(array_base(QuantumIdentifier))
