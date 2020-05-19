import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class Index:
    """ ARTS Index data

    Properties:
        val:
            a value (int)
    """
    def __init__(self, value):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createIndex())
            self.val = value

    @staticmethod
    def name():
        return "Index"

    @property
    def val(self):
        """ a value (int) """
        return lib.getIndex(self.__data__)

    @val.setter
    def val(self, x):
        lib.setIndex(self.__data__, int(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printIndex(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteIndex(self.__data__)

    def __repr__(self):
        return "{}".format(self.val)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Index):
            self.val = other.val
        else:
            self.val = other

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadIndex(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveIndex(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))
    
    def netcdfify(self):
        """ Create the NETCDF4 information required for writing this data
        
        Output: list that can be processed by netcdf.py, False arraytype
        """
        return [["val", self.val, int, {}]], False
    
    def denetcdf(self, group):
        """ Sets this based on a netcdf group
        
        Input:
            Group of data that can be interpreted as this's information
        """
        self.val = group.val

    def __int__(self):
        return self.val

    def __float__(self):
        return float(self.val)

    def __eq__(self, other):
        return self.val == int(other)

    def __lt__(self, other):
        return self.val < int(other)

    def __le__(self, other):
        return self.val <= int(other)

    def __gt__(self, other):
        return self.val > int(other)

    def __ge__(self, other):
        return self.val >= int(other)

    def __iadd__(self, val):
        self.val += int(val)
        return self

    def __isub__(self, val):
        self.val -= int(val)
        return self

    def __imul__(self, val):
        self.val *= int(val)
        return self

    def __itruediv__(self, val):
        self.val /= int(val)
        return self

    def __ipow__(self, val):
        self.val **= int(val)
        return self

    def __add__(self, val):
        return Index(self.val + int(val))

    def __sub__(self, val):
        return Index(self.val - int(val))

    def __mul__(self, val):
        return Index(self.val * int(val))

    def __truediv__(self, val):
        return Index(self.val / int(val))

    def __pow__(self, val):
        return Index(self.val ** int(val))

    def __radd__(self, val):
        return Index(int(val) + self.val)

    def __rsub__(self, val):
        return Index(int(val) - self.val)

    def __rmul__(self, val):
        return Index(int(val) * self.val)

    def __rtruediv__(self, val):
        return Index(int(val) / self.val)

    def __rpow__(self, val):
        return Index(int(val) ** self.val)


class Numeric:
    """ ARTS Numeric data

    Properties:
        val:
            a value (float)
    """
    def __init__(self, value):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createNumeric())
            self.val = value

    @staticmethod
    def name():
        return "Numeric"

    @property
    def val(self):
        """ a value (float) """
        return lib.getNumeric(self.__data__)

    @val.setter
    def val(self, x):
        lib.setNumeric(self.__data__, float(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printNumeric(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteNumeric(self.__data__)

    def __repr__(self):
        return "{}".format(self.val)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Numeric):
            self.val = other.val
        else:
            self.val = other

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadNumeric(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveNumeric(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __int__(self):
        return int(self.val)

    def __float__(self):
        return self.val

    def __eq__(self, other):
        return self.val == float(other)

    def __lt__(self, other):
        return self.val < float(other)

    def __le__(self, other):
        return self.val <= float(other)

    def __gt__(self, other):
        return self.val > float(other)

    def __ge__(self, other):
        return self.val >= float(other)

    def __iadd__(self, val):
        self.val += float(val)
        return self

    def __isub__(self, val):
        self.val -= float(val)
        return self

    def __imul__(self, val):
        self.val *= float(val)
        return self

    def __itruediv__(self, val):
        self.val /= float(val)
        return self

    def __ipow__(self, val):
        self.val **= float(val)
        return self

    def __add__(self, val):
        return Numeric(self.val + float(val))

    def __sub__(self, val):
        return Numeric(self.val - float(val))

    def __mul__(self, val):
        return Numeric(self.val * float(val))

    def __truediv__(self, val):
        return Numeric(self.val / float(val))

    def __pow__(self, val):
        return Numeric(self.val ** float(val))

    def __radd__(self, val):
        return Numeric(float(val) + self.val)

    def __rsub__(self, val):
        return Numeric(float(val) - self.val)

    def __rmul__(self, val):
        return Numeric(float(val) * self.val)

    def __rtruediv__(self, val):
        return Numeric(float(val) / self.val)

    def __rpow__(self, val):
        return Numeric(float(val) ** self.val)


class String:
    """ ARTS String data

    Properties:
        val:
            a value (str)
    """
    def __init__(self, value):
        if isinstance(value, c.c_void_p):
            self.__delete__ = False
            self.__data__ = value
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createString())
            self.val = value

    @staticmethod
    def name():
        return "String"

    @property
    def val(self):
        """ a value (str) """
        x = lib.getString(self.__data__)
        return x.decode("utf-8") if x else ""

    @val.setter
    def val(self, x):
        lib.setString(self.__data__, str(x).encode("utf-8"))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printString(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteString(self.__data__)

    def __repr__(self):
        return "{}".format(self.val)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, String):
            self.val = other.val
        else:
            self.val = other

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadString(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveString(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))
    
    def netcdfify(self):
        """ Create the NETCDF4 information required for writing this data
        
        Output: list that can be processed by netcdf.py, False arraytype
        """
        return [["val", self.val, str, {}]], False
    
    def denetcdf(self, group):
        """ Sets this based on a netcdf group
        
        Input:
            Group of data that can be interpreted as this's information
        """
        self.val = group.val

    def __eq__(self, other):
        return self.val == str(other)

    def __iadd__(self, val):
        self.val += str(val)
        return self

    def __add__(self, val):
        return String(self.val + str(val))

    def __radd__(self, val):
        return String(str(val) + self.val)


exec(array_base(Index))


exec(array_base(ArrayOfIndex))


exec(array_base(String))


exec(array_base(ArrayOfString))


lib.createIndex.restype = c.c_void_p
lib.createIndex.argtypes = []

lib.deleteIndex.restype = None
lib.deleteIndex.argtypes = [c.c_void_p]

lib.printIndex.restype = None
lib.printIndex.argtypes = [c.c_void_p]

lib.getIndex.restype = c.c_long
lib.getIndex.argtypes = [c.c_void_p]

lib.setIndex.restype = None
lib.setIndex.argtypes = [c.c_void_p, c.c_long]

lib.xmlreadIndex.restype = c.c_long
lib.xmlreadIndex.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveIndex.restype = c.c_long
lib.xmlsaveIndex.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]


lib.createNumeric.restype = c.c_void_p
lib.createNumeric.argtypes = []

lib.deleteNumeric.restype = None
lib.deleteNumeric.argtypes = [c.c_void_p]

lib.printNumeric.restype = None
lib.printNumeric.argtypes = [c.c_void_p]

lib.getNumeric.restype = c.c_double
lib.getNumeric.argtypes = [c.c_void_p]

lib.setNumeric.restype = None
lib.setNumeric.argtypes = [c.c_void_p, c.c_double]

lib.xmlreadNumeric.restype = c.c_long
lib.xmlreadNumeric.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveNumeric.restype = c.c_long
lib.xmlsaveNumeric.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]


lib.createString.restype = c.c_void_p
lib.createString.argtypes = []

lib.deleteString.restype = None
lib.deleteString.argtypes = [c.c_void_p]

lib.printString.restype = None
lib.printString.argtypes = [c.c_void_p]

lib.getString.restype = c.c_char_p
lib.getString.argtypes = [c.c_void_p]

lib.setString.restype = None
lib.setString.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlreadString.restype = c.c_long
lib.xmlreadString.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveString.restype = c.c_long
lib.xmlsaveString.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
