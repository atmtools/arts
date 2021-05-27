import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base
from pyarts.classes.BasicTypes import Index, Numeric, String
from pyarts.classes.SpeciesIsotopeRecord import Species, SpeciesIsotopeRecord


class SpeciesTagType:
    """ ARTS SpeciesTagType data
    
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
            self.__data__ = c.c_void_p(lib.createSpeciesTagType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getSpeciesTagTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setSpeciesTagTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad SpeciesTagType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSpeciesTagType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesTagType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, SpeciesTagType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, SpeciesTagType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createSpeciesTagType.restype = c.c_void_p
lib.createSpeciesTagType.argtypes = []

lib.deleteSpeciesTagType.restype = None
lib.deleteSpeciesTagType.argtypes = [c.c_void_p]

lib.printSpeciesTagType.restype = None
lib.printSpeciesTagType.argtypes = [c.c_void_p]

lib.getSpeciesTagTypeString.restype = c.c_void_p
lib.getSpeciesTagTypeString.argtypes = [c.c_void_p]

lib.setSpeciesTagTypeString.restype = c.c_int
lib.setSpeciesTagTypeString.argtypes = [c.c_void_p, c.c_char_p]

class SpeciesTag:
    """ ARTS SpeciesTag data

    Properties:
        spec_ind:
            Species of tag (Index)
            
        lower_freq:
            Lower frequency range (Numeric)

        upper_freq:
            Upper frequency range (Numeric)
        
        type:
            Tag type (SpeciesTagType)
        
        cia_2nd_species:
            Second species of CIA broadening (Species; if defined)
        
        cia_dataset_index:
            CIA model flag (Index)
            
        isot:
            The ARTS isotopologue (const SpeciesIsotopeRecord)
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSpeciesTag())
            if data is not None:
                self.setFromString(data)
                assert self.OK, "Fail initialize properly"

    def setFromString(self, data):
        """ Set the species tag from a string

        Input:
            data:
                Species tag data (str)
        """
        if isinstance(data, str):
            data = data.encode("utf-8")
            if lib.setSpeciesTag(self.__data__, data):
                raise ValueError(f"Invalid species tag {data}")
        else:
            raise TypeError("Expects str input")
    
    @property
    def OK(self):
        if self.spec_ind < 0: return False
        if self.spec_ind >= SpeciesIsotopeRecord.max_len(): return False
        return True
    
    @property
    def spec_ind(self):
        return Index(c.c_void_p(lib.getspec_indSpeciesTag(self.__data__)))
    
    @spec_ind.setter
    def spec_ind(self, x):
        x = int(x)
        if x >= 0 and x < SpeciesIsotopeRecord.max_len():
            self.spec_ind.set(x)
    
    @property
    def lower_freq(self):
        return Numeric(c.c_void_p(lib.getlower_freqSpeciesTag(self.__data__)))
    
    @lower_freq.setter
    def lower_freq(self, x):
        self.lower_freq.set(x)
    
    @property
    def upper_freq(self):
        return Numeric(c.c_void_p(lib.getupper_freqSpeciesTag(self.__data__)))
    
    @upper_freq.setter
    def upper_freq(self, x):
        self.upper_freq.set(x)
    
    @property
    def type(self):
        return SpeciesTagType(c.c_void_p(lib.gettypeSpeciesTag(self.__data__)))
    
    @type.setter
    def type(self, x):
        self.type.set(x)
    
    @property
    def cia_2nd_species(self):
        return Species(c.c_void_p(lib.getcia_2nd_speciesSpeciesTag(self.__data__)))
    
    @cia_2nd_species.setter
    def cia_2nd_species(self, x):
        self.cia_2nd_species.set(x)
    
    @property
    def cia_dataset_index(self):
        return Index(c.c_void_p(lib.getcia_dataset_indexSpeciesTag(self.__data__)))
    
    @cia_dataset_index.setter
    def cia_dataset_index(self, x):
        self.cia_dataset_index.set(x)
    
    @property
    def isot(self):
        return SpeciesIsotopeRecord.from_index(self.spec_ind)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        if self.OK:
            lib.printSpeciesTag(self.__data__)
        else:
            raise RuntimeError("Cannot print SpeciesTag in bad state")

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesTag(self.__data__)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, SpeciesTag):
            self.spec_ind = other.spec_ind
            self.type = other.type
            self.lower_freq = other.lower_freq
            self.upper_freq = other.upper_freq
            self.type = other.type
            if self.type == "Cia":
                self.cia_2nd_species = other.cia_2nd_species
                self.cia_dataset_index = other.cia_dataset_index
        else:
            raise TypeError("Expects SpeciesTag")
    
    @property
    def as_string(self):
        """ Returns the name as a String """
        if self.OK:
            return String(c.c_void_p(lib.getNameSpeciesTag(self.__data__)), delete=True)
        else:
            return String("Tag is in bad shape")

    def __repr__(self):
        return f"{self.as_string}"

    @staticmethod
    def name():
        return "SpeciesTag"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadSpeciesTag(self.__data__, correct_read_arguments(file)):
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
        if not self.OK:
            raise RuntimeError("Cannot save SpeciesTag in bad state")

        if lib.xmlsaveSpeciesTag(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, SpeciesTag) and \
            self.spec_ind == other.spec_ind and \
            self.lower_freq == other.lower_freq and \
            self.upper_freq == other.upper_freq and \
            self.type == other.type and \
            self.cia_2nd_species == other.cia_2nd_species and \
            self.cia_dataset_index == other.cia_dataset_index:
            return True
        else:
            return False


# Generate ArrayOfSpeciesTag
exec(array_base(SpeciesTag))


# Generate ArrayOfArrayOfSpeciesTag
exec(array_base(ArrayOfSpeciesTag))


class ArrayOfSpeciesTagAdapted(ArrayOfSpeciesTag):
    @property
    def data(self):
        return super().data

    @data.setter
    def data(self, val):
        if isinstance(val, str) and val == "":
            return
        if isinstance(val, str):
            val = [s.strip() for s in val.split(",")]
        super(self.__class__, self.__class__).data.fset(self, val)

    def __str__(self):
        return '"' + ", ".join([str(s.as_string) for s in self.data]) + '"'


ArrayOfSpeciesTag = ArrayOfSpeciesTagAdapted
ArrayOfSpeciesTag.__name__ = "ArrayOfSpeciesTag"

lib.createSpeciesTag.restype = c.c_void_p
lib.createSpeciesTag.argtypes = []

lib.deleteSpeciesTag.restype = None
lib.deleteSpeciesTag.argtypes = [c.c_void_p]

lib.printSpeciesTag.restype = None
lib.printSpeciesTag.argtypes = [c.c_void_p]

lib.xmlreadSpeciesTag.restype = c.c_long
lib.xmlreadSpeciesTag.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveSpeciesTag.restype = c.c_long
lib.xmlsaveSpeciesTag.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.setSpeciesTag.restype = c.c_long
lib.setSpeciesTag.argtypes = [c.c_void_p, c.c_char_p]

lib.getNameSpeciesTag.restype = c.c_void_p
lib.getNameSpeciesTag.argtypes = [c.c_void_p]

lib.getspec_indSpeciesTag.restype = c.c_void_p
lib.getspec_indSpeciesTag.argtypes = [c.c_void_p]

lib.getlower_freqSpeciesTag.restype = c.c_void_p
lib.getlower_freqSpeciesTag.argtypes = [c.c_void_p]

lib.getupper_freqSpeciesTag.restype = c.c_void_p
lib.getupper_freqSpeciesTag.argtypes = [c.c_void_p]

lib.gettypeSpeciesTag.restype = c.c_void_p
lib.gettypeSpeciesTag.argtypes = [c.c_void_p]

lib.getcia_2nd_speciesSpeciesTag.restype = c.c_void_p
lib.getcia_2nd_speciesSpeciesTag.argtypes = [c.c_void_p]

lib.getcia_dataset_indexSpeciesTag.restype = c.c_void_p
lib.getcia_dataset_indexSpeciesTag.argtypes = [c.c_void_p]
