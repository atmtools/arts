import ctypes as c
from pyarts.workspace.api import arts_api as lib
from pyarts.classes.io import correct_save_arguments, correct_read_arguments

from pyarts.classes.Matrix import Matrix
from pyarts.classes.BasicTypes import String


class PartitionFunctionsType:
    """ ARTS PartitionFunctionsData type data
    
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
            self.__data__ = c.c_void_p(lib.createPartitionFunctionsType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getPartitionFunctionsTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setPartitionFunctionsTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad PartitionFunctionsType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printPartitionFunctionsType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deletePartitionFunctionsType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, PartitionFunctionsType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, PartitionFunctionsType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createPartitionFunctionsType.restype = c.c_void_p
lib.createPartitionFunctionsType.argtypes = []

lib.deletePartitionFunctionsType.restype = None
lib.deletePartitionFunctionsType.argtypes = [c.c_void_p]

lib.printPartitionFunctionsType.restype = None
lib.printPartitionFunctionsType.argtypes = [c.c_void_p]

lib.getPartitionFunctionsTypeString.restype = c.c_void_p
lib.getPartitionFunctionsTypeString.argtypes = [c.c_void_p]

lib.setPartitionFunctionsTypeString.restype = c.c_int
lib.setPartitionFunctionsTypeString.argtypes = [c.c_void_p, c.c_char_p]


class PartitionFunctionsData:
    """ ARTS PartitionFunctionsData data

    Properties:
            
        data:
            The data for the partition function (Matrix)
        type:
            The type of identifier (PartitionFunctionsType)
        
        The shape of the data depends on the type
        
            Interp:
                The row gives the number of temperatures
                
                The col is 2-long, with data[:, 0] giving the temperature and
                data[:, 1] gives the data
                
                The valid row-count must be more than 1
                
            Coeff:
                The row gives the coefficient for the polynominal fit.  The
                first value is the constant, the second the linear term and so
                on...
                
                The col is 1-long, it has no meaning
    """
    def __init__(self, type="Interp", data=Matrix()):
        if isinstance(type, c.c_void_p):
            self.__delete__ = False
            self.__data__ = type
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createPartitionFunctionsData())
            self.type = type
            self.data = data
    
    @property
    def type(self):
        return PartitionFunctionsType(c.c_void_p(lib.gettypePartitionFunctionsData(self.__data__)))
    
    @type.setter
    def type(self, x):
        self.type.set(x)
    
    @property
    def data(self):
        return Matrix(c.c_void_p(lib.getdataPartitionFunctionsData(self.__data__)))
    
    @data.setter
    def data(self, x):
        self.data.set(x)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printPartitionFunctionsData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deletePartitionFunctionsData(self.__data__)
    
    def __repr__(self):
        return "ARTS PartitionFunctionsData"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, PartitionFunctionsData):
            self.type = other.type
            self.data = other.data
        else:
            raise TypeError("Expects PartitionFunctionsData")

    @staticmethod
    def name():
        return "PartitionFunctionsData"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadPartitionFunctionsData(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsavePartitionFunctionsData(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        return self.type == other.type and self.data == other.data


lib.createPartitionFunctionsData.restype = c.c_void_p
lib.createPartitionFunctionsData.argtypes = []

lib.deletePartitionFunctionsData.restype = None
lib.deletePartitionFunctionsData.argtypes = [c.c_void_p]

lib.printPartitionFunctionsData.restype = None
lib.printPartitionFunctionsData.argtypes = [c.c_void_p]

lib.xmlreadPartitionFunctionsData.restype = c.c_long
lib.xmlreadPartitionFunctionsData.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsavePartitionFunctionsData.restype = c.c_long
lib.xmlsavePartitionFunctionsData.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.gettypePartitionFunctionsData.restype = c.c_void_p
lib.gettypePartitionFunctionsData.argtypes = [c.c_void_p]

lib.getdataPartitionFunctionsData.restype = c.c_void_p
lib.getdataPartitionFunctionsData.argtypes = [c.c_void_p]
