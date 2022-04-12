import ctypes as c
from pyarts.workspace.api import arts_api as lib
from pyarts.classes.ArrayBase import array_base

from pyarts.classes.BasicTypes import String

class Species:
    """ ARTS species
    
    Properties:
        short_name:
            Short name of the species (str)
        
        full_name:
            Full name of the species (str)
    """

    def __init__(self, spec="", short_name=True, delete=False):
        if isinstance(spec, c.c_void_p):
            self.__delete__ = delete
            self.__data__ = spec
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSpecies())
            if short_name:
                self.short_name = spec
            else:
                self.full_name = spec
    
    @property
    def short_name(self):
        return String(c.c_void_p(lib.getSpeciesShortName(self.__data__)),
                      delete=True)
    
    @short_name.setter
    def short_name(self, x):
        if lib.setSpeciesShortName(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad name {x}")
    
    @property
    def full_name(self):
        return String(c.c_void_p(lib.getSpeciesLongName(self.__data__)),
                      delete=True)
    
    @full_name.setter
    def full_name(self, x):
        if lib.setSpeciesLongName(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad name {x}")
    
    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSpecies(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSpecies(self.__data__)

    def __repr__(self):
        return f"{self.short_name}"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Species):
            self.short_name = other.short_name
        else:
            raise TypeError("Expects Species")

    def __eq__(self, other):
        return self.short_name == other.short_name
    
    @staticmethod
    def name():
        return "Species"
    

lib.createSpecies.restype = c.c_void_p
lib.createSpecies.argtypes = []

lib.deleteSpecies.restype = None
lib.deleteSpecies.argtypes = [c.c_void_p]

lib.printSpecies.restype = None
lib.printSpecies.argtypes = [c.c_void_p]

lib.getSpeciesShortName.restype = c.c_void_p
lib.getSpeciesShortName.argtypes = [c.c_void_p]

lib.setSpeciesShortName.restype = c.c_int
lib.setSpeciesShortName.argtypes = [c.c_void_p, c.c_void_p]

lib.getSpeciesLongName.restype = c.c_void_p
lib.getSpeciesLongName.argtypes = [c.c_void_p]

lib.setSpeciesLongName.restype = c.c_int
lib.setSpeciesLongName.argtypes = [c.c_void_p, c.c_void_p]


# Generate ArrayOfSpecies
exec(array_base(Species))


class SpeciesIsotopeRecord:
    """ ARTS isotopologues

    Properties:
        spec:
            Species species (const Species)

        isotname:
            Species isotopologue (const String)
            
        mass:
            Mass of isotope (const Numeric)
            Requires that spec_ind >= 0
        
        gi:
            Degeneracy of isotope (const Index)
            Requires that spec_ind >= 0
        
        spec_ind:
            Equivalent species index in ARTS (Index)
        """

    def __init__(self, spec="", isotname="*"):
        if isinstance(spec, c.c_void_p):
            self.__delete__ = False
            self.__data__ = spec
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSpeciesIsotopeRecord())
            if '-' in spec:
                self.index = self.get_index(spec)
            else:
                self.index = self.get_index_split(spec, isotname)
    
    @staticmethod
    def from_index(ind):
        """ Creates a SpeciesIsotopeRecord from index """
        x = SpeciesIsotopeRecord("H2O")
        x.index = ind
        return x

    @property
    def spec(self):
        """ Species species (const Species) """
        return Species(c.c_void_p(lib.getSpeciesSpeciesIsotopeRecord(self.__data__)),
                       delete=True)

    @property
    def isotname(self):
        return String(c.c_void_p(lib.getIsotnameSpeciesIsotopeRecord(self.__data__)),
                       delete=True)
    
    @property
    def index(self):
        """ Equivalent species index in ARTS (Index) """
        return lib.getIndexSpeciesIsotopeRecordFromData(self.__data__)
    
    @index.setter
    def index(self, index):
        if lib.setSpeciesIsotopeRecordToIndex(self.__data__, int(index)):
            raise RuntimeError(f"Index {index} out of bounds of [{0}, {self.max_len()})")
    
    @property
    def mass(self):
        """ Mass of isotope (const Numeric) """
        return lib.getMassSpeciesIsotopeRecord(self.__data__)
    
    @property
    def gi(self):
        """ Degeneracy of isotope (const Index) """
        return lib.getGSpeciesIsotopeRecord(self.__data__)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSpeciesIsotopeRecord(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesIsotopeRecord(self.__data__)

    def __repr__(self):
        return f"{self.spec}-{self.isotname}"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        self.index = other.index

    def __eq__(self, other):
        return self.index == other.index
    
    @staticmethod
    def max_len():
        return lib.nelemSpeciesIsotopeRecordDefined()
    
    @staticmethod
    def get_index_split(spec, isot):
        return lib.getIndexSpeciesIsotopeRecordFromNames(str(spec).encode("utf-8"),
                                                         str(isot).encode("utf-8"))
    @staticmethod
    def get_index(spec):
        return lib.getIndexSpeciesIsotopeRecordFromFullName(str(spec).encode("utf-8"))

lib.createSpeciesIsotopeRecord.restype = c.c_void_p
lib.createSpeciesIsotopeRecord.argtypes = []

lib.deleteSpeciesIsotopeRecord.restype = None
lib.deleteSpeciesIsotopeRecord.argtypes = [c.c_void_p]

lib.printSpeciesIsotopeRecord.restype = None
lib.printSpeciesIsotopeRecord.argtypes = [c.c_void_p]

lib.getIndexSpeciesIsotopeRecordFromNames.restype = c.c_long
lib.getIndexSpeciesIsotopeRecordFromNames.argtypes = [c.c_char_p, c.c_char_p]

lib.getIndexSpeciesIsotopeRecordFromFullName.restype = c.c_long
lib.getIndexSpeciesIsotopeRecordFromFullName.argtypes = [c.c_char_p]

lib.getIndexSpeciesIsotopeRecordFromData.restype = c.c_long
lib.getIndexSpeciesIsotopeRecordFromData.argtypes = [c.c_void_p]

lib.setSpeciesIsotopeRecordToIndex.restype = c.c_int
lib.setSpeciesIsotopeRecordToIndex.argtypes = [c.c_void_p, c.c_long]

lib.getSpeciesSpeciesIsotopeRecord.restype = c.c_void_p
lib.getSpeciesSpeciesIsotopeRecord.argtypes = [c.c_void_p]

lib.getIsotnameSpeciesIsotopeRecord.restype = c.c_void_p
lib.getIsotnameSpeciesIsotopeRecord.argtypes = [c.c_void_p]

lib.getMassSpeciesIsotopeRecord.restype = c.c_double
lib.getMassSpeciesIsotopeRecord.argtypes = [c.c_void_p]

lib.getGSpeciesIsotopeRecord.restype = c.c_long
lib.getGSpeciesIsotopeRecord.argtypes = [c.c_void_p]

lib.nelemSpeciesIsotopeRecordDefined.restype = c.c_long
lib.nelemSpeciesIsotopeRecordDefined.argtypes = []
