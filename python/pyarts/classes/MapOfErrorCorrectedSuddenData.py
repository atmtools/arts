import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.BasicTypes import Numeric
from pyarts.classes.SpeciesIsotopeRecord import Species


class SpeciesErrorCorrectedSuddenData:
    """ ARTS SpeciesErrorCorrectedSuddenData data

    Properties:
        spec:
            Species
        
        a:
            Constant 1
        
        b:
            Constant 2
        
        gamma:
            Constant 3
        
        dc:
            Constant 4
        
        mass:
            Constant 5
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSpeciesErrorCorrectedSuddenData())
            assert data is None, "Fail initialize properly"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSpeciesErrorCorrectedSuddenData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesErrorCorrectedSuddenData(self.__data__)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        raise RuntimeError("Not implemented")

    def __repr__(self):
        return f"{self.spec} {self.a} {self.b} {self.gamma} {self.dc} {self.mass}"

    def __eq__(self, other):
        raise RuntimeError("Not implemented")
        
    @property
    def spec(self):
        return Species(c.c_void_p(lib.getspecSpeciesErrorCorrectedSuddenData(self.__data__)))
    
    @spec.setter
    def spec(self, val):
        self.spec.set(val)
        
    @property
    def a(self):
        return Numeric(c.c_void_p(lib.getaSpeciesErrorCorrectedSuddenData(self.__data__)))
    
    @a.setter
    def a(self, val):
        self.a.set(val)
        
    @property
    def gamma(self):
        return Numeric(c.c_void_p(lib.getgammaSpeciesErrorCorrectedSuddenData(self.__data__)))
    
    @gamma.setter
    def gamma(self, val):
        self.gamma.set(val)
        
    @property
    def b(self):
        return Numeric(c.c_void_p(lib.getbSpeciesErrorCorrectedSuddenData(self.__data__)))
    
    @b.setter
    def b(self, val):
        self.b.set(val)
        
    @property
    def dc(self):
        return Numeric(c.c_void_p(lib.getdcSpeciesErrorCorrectedSuddenData(self.__data__)))
    
    @dc.setter
    def dc(self, val):
        self.dc.set(val)
        
    @property
    def mass(self):
        return Numeric(c.c_void_p(lib.getmassSpeciesErrorCorrectedSuddenData(self.__data__)))
        
    @mass.setter
    def mass(self, val):
        self.mass.set(val)

lib.createSpeciesErrorCorrectedSuddenData.restype = c.c_void_p
lib.createSpeciesErrorCorrectedSuddenData.argtypes = []

lib.deleteSpeciesErrorCorrectedSuddenData.restype = None
lib.deleteSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.printSpeciesErrorCorrectedSuddenData.restype = None
lib.printSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getspecSpeciesErrorCorrectedSuddenData.restype = c.c_void_p
lib.getspecSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getaSpeciesErrorCorrectedSuddenData.restype = c.c_void_p
lib.getaSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getgammaSpeciesErrorCorrectedSuddenData.restype = c.c_void_p
lib.getgammaSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getbSpeciesErrorCorrectedSuddenData.restype = c.c_void_p
lib.getbSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getdcSpeciesErrorCorrectedSuddenData.restype = c.c_void_p
lib.getdcSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getmassSpeciesErrorCorrectedSuddenData.restype = c.c_void_p
lib.getmassSpeciesErrorCorrectedSuddenData.argtypes = [c.c_void_p]


class ErrorCorrectedSuddenData:
    """ ARTS ErrorCorrectedSuddenData data

    Properties:
        id:
            Identity of the band/species/isotopologue
        
        Access operator with Species
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createErrorCorrectedSuddenData())
            assert data is None, "Fail initialize properly"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printErrorCorrectedSuddenData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteErrorCorrectedSuddenData(self.__data__)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        raise RuntimeError("Not implemented")

    def __repr__(self):
        x = ""
        x += f"{self.id}\n"
        for i in range(lib.getnelemErrorCorrectedSuddenData(self.__data__)):
            ptr = lib.getSpeciesErrorCorrectedSuddenDataAtErrorCorrectedSuddenData(self.__data__, i)
            x += f"\t{SpeciesErrorCorrectedSuddenData(c.c_void_p(ptr))}\n"
        return x

    def __eq__(self, other):
        raise RuntimeError("Not implemented")
        
    @property
    def id(self):
        return QuantumIdentifier(c.c_void_p(lib.getidErrorCorrectedSuddenData(self.__data__)))
        
    @id.setter
    def id(self, val):
        self.id.set(val)
        
    def __getitem__(self, key):
        if not isinstance(key, Species):
            key = Species(key)
        return SpeciesErrorCorrectedSuddenData(c.c_void_p(
            lib.getErrorCorrectedSuddenData(self.__data__, key.__data__)))
        
    def __setitem__(self, key, val):
        self[key].set(val)

lib.createErrorCorrectedSuddenData.restype = c.c_void_p
lib.createErrorCorrectedSuddenData.argtypes = []

lib.deleteErrorCorrectedSuddenData.restype = None
lib.deleteErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.printErrorCorrectedSuddenData.restype = None
lib.printErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getidErrorCorrectedSuddenData.restype = c.c_void_p
lib.getidErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getErrorCorrectedSuddenData.restype = c.c_void_p
lib.getErrorCorrectedSuddenData.argtypes = [c.c_void_p, c.c_void_p]

lib.getnelemErrorCorrectedSuddenData.restype = c.c_long
lib.getnelemErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getSpeciesErrorCorrectedSuddenDataAtErrorCorrectedSuddenData.restype = c.c_void_p
lib.getSpeciesErrorCorrectedSuddenDataAtErrorCorrectedSuddenData.argtypes = [c.c_void_p, c.c_long]


class MapOfErrorCorrectedSuddenData:
    """ ARTS MapOfErrorCorrectedSuddenData data

    Properties:
        Access operator with QuantumIdentifier
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createMapOfErrorCorrectedSuddenData())
            assert data is None, "Fail initialize properly"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printMapOfErrorCorrectedSuddenData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteMapOfErrorCorrectedSuddenData(self.__data__)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        raise RuntimeError("Not implemented")

    @staticmethod
    def name():
        return "MapOfErrorCorrectedSuddenData"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadMapOfErrorCorrectedSuddenData(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveMapOfErrorCorrectedSuddenData(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        raise RuntimeError("Not implemented")
        
    def __getitem__(self, key):
        if not isinstance(key, QuantumIdentifier):
            key = QuantumIdentifier(key)
        return ErrorCorrectedSuddenData(c.c_void_p(
            lib.getMapOfErrorCorrectedSuddenData(self.__data__, key.__data__)))
        
    def __setitem__(self, key, val):
        self[key].set(val)

    def __repr__(self):
        x = ""
        for i in range(lib.getnelemMapOfErrorCorrectedSuddenData(self.__data__)):
            ptr = lib.getErrorCorrectedSuddenDataAtMapOfErrorCorrectedSuddenData(self.__data__, i)
            x += f"{ErrorCorrectedSuddenData(c.c_void_p(ptr))}\n"
        return x

lib.createMapOfErrorCorrectedSuddenData.restype = c.c_void_p
lib.createMapOfErrorCorrectedSuddenData.argtypes = []

lib.deleteMapOfErrorCorrectedSuddenData.restype = None
lib.deleteMapOfErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.printMapOfErrorCorrectedSuddenData.restype = None
lib.printMapOfErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.xmlreadMapOfErrorCorrectedSuddenData.restype = c.c_long
lib.xmlreadMapOfErrorCorrectedSuddenData.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveMapOfErrorCorrectedSuddenData.restype = c.c_long
lib.xmlsaveMapOfErrorCorrectedSuddenData.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getMapOfErrorCorrectedSuddenData.restype = c.c_void_p
lib.getMapOfErrorCorrectedSuddenData.argtypes = [c.c_void_p, c.c_void_p]

lib.getnelemMapOfErrorCorrectedSuddenData.restype = c.c_long
lib.getnelemMapOfErrorCorrectedSuddenData.argtypes = [c.c_void_p]

lib.getErrorCorrectedSuddenDataAtMapOfErrorCorrectedSuddenData.restype = c.c_void_p
lib.getErrorCorrectedSuddenDataAtMapOfErrorCorrectedSuddenData.argtypes = [c.c_void_p, c.c_long]
