import ctypes as c
from pyarts.workspace.api import arts_api as lib


class SpeciesTag:
    """ ARTS SpeciesTag data

    Can only be set from string using setFromString(str) method or from
    other SpeciesTag using the set(SpeciesTag) method.  Using other methods
    can easily lead to a poorly defined tag

    Properties:
        species:
            Species of tag (const Index)

        isotopologue:
            Isotopologue of tag (const Index)

        type:
            Type of tag (const Index)

        cia_second:
            Secondary species when CIA mode (const Index)

        cia_dataset:
            CIA dataset when CIA mode (const Index)

        lf:
            Lower frequency range (const Numeric)

        uf:
            Upper frequency range (const Numeric)
    """
    def __init__(self, data):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSpeciesTag())
            self.setFromString(data)

    def setFromString(self, data):
        """ Set the species tag from a string

        Note that this is the only valid way to set a species tag but that
        the bare set() function can copy another valid species tag

        Input:
            data:
                Species tag data (str)
        """
        if isinstance(data, str):
            data = data.encode("ascii")
            if lib.setSpeciesTag(self.__data__, data):
                raise ValueError("Invalid species tag")
        else:
            raise TypeError("Expects str input")

    @property
    def species(self):
        """ Species of tag (const Index) """
        return lib.getSpeciesTagSpecies(self.__data__)

    @property
    def isotopologue(self):
        """ Isotopologue of tag (const Index) """
        return lib.getSpeciesTagIsotopologue(self.__data__)

    @property
    def type(self):
        """ Type of tag (const Index) """
        return lib.getSpeciesTagType(self.__data__)

    @property
    def cia_second(self):
        """ Secondary species when CIA mode (const Index) """
        return lib.getSpeciesTagCIASecond(self.__data__)

    @property
    def cia_dataset(self):
        """ CIA dataset when CIA mode (const Index) """
        return lib.getSpeciesTagCIADataset(self.__data__)

    @property
    def lf(self):
        """ Lower frequency range (const Numeric) """
        return lib.getSpeciesTagLowerFrequency(self.__data__)

    @property
    def uf(self):
        """ Upper frequency range (const Numeric) """
        return lib.getSpeciesTagUpperFrequency(self.__data__)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSpeciesTag(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesTag(self.__data__)

    def __repr__(self):
        return "ARTS SpeciesTag"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, SpeciesTag):
            lib.setSpeciesTagType(self.__data__, other.type)
            lib.setSpeciesTagSpecies(self.__data__, other.species)
            lib.setSpeciesTagIsotopologue(self.__data__, other.isotopologue)
            lib.setSpeciesTagCIASecond(self.__data__, other.cia_second)
            lib.setSpeciesTagCIADataset(self.__data__, other.cia_dataset)
            lib.setSpeciesTagLowerFrequency(self.__data__, other.lf)
            lib.setSpeciesTagUpperFrequency(self.__data__, other.uf)
        else:
            raise TypeError("Expects SpeciesTag")

    @staticmethod
    def validSpecies(spec):
        """ Returns whether species is a valid species according to ARTS

        Input:
            spec:
                Species (Index)

        Output:
            Boolean True or False
        """
        spec = int(spec)
        if lib.validSpecies(spec) == 0:
            return True
        else:
            return False

    @staticmethod
    def validIsotopologue(spec, isot):
        """ Returns whether species is a valid species and an isotopologue
        is a valid isotopologue to ARTS

        Input:
            spec:
                Species (Index)

            isot:
                Isotopologue (Index)

        Output:
            Boolean True or False
        """
        spec = int(spec)
        isot = int(isot)
        if SpeciesTag.validSpecies(spec):
            if lib.validIsotopologue(spec, isot) == 0:
                return True
            else:
                return False
        else:
            return False


lib.createSpeciesTag.restype = c.c_void_p
lib.createSpeciesTag.argtypes = []

lib.deleteSpeciesTag.restype = None
lib.deleteSpeciesTag.argtypes = [c.c_void_p]

lib.printSpeciesTag.restype = None
lib.printSpeciesTag.argtypes = [c.c_void_p]

lib.setSpeciesTag.restype = None
lib.setSpeciesTag.argtypes = [c.c_void_p, c.c_char_p]

lib.getSpeciesTagSpecies.restype = c.c_long
lib.getSpeciesTagSpecies.argtypes = [c.c_void_p]

lib.getSpeciesTagIsotopologue.restype = c.c_long
lib.getSpeciesTagIsotopologue.argtypes = [c.c_void_p]

lib.getSpeciesTagType.restype = c.c_long
lib.getSpeciesTagType.argtypes = [c.c_void_p]

lib.getSpeciesTagCIASecond.restype = c.c_long
lib.getSpeciesTagCIASecond.argtypes = [c.c_void_p]

lib.getSpeciesTagCIADataset.restype = c.c_long
lib.getSpeciesTagCIADataset.argtypes = [c.c_void_p]

lib.getSpeciesTagLowerFrequency.restype = c.c_double
lib.getSpeciesTagLowerFrequency.argtypes = [c.c_void_p]

lib.getSpeciesTagUpperFrequency.restype = c.c_double
lib.getSpeciesTagUpperFrequency.argtypes = [c.c_void_p]

lib.setSpeciesTagSpecies.restype = None
lib.setSpeciesTagSpecies.argtypes = [c.c_void_p, c.c_long]

lib.setSpeciesTagIsotopologue.restype = None
lib.setSpeciesTagIsotopologue.argtypes = [c.c_void_p, c.c_long]

lib.setSpeciesTagType.restype = None
lib.setSpeciesTagType.argtypes = [c.c_void_p, c.c_long]

lib.setSpeciesTagCIASecond.restype = None
lib.setSpeciesTagCIASecond.argtypes = [c.c_void_p, c.c_long]

lib.setSpeciesTagCIADataset.restype = None
lib.setSpeciesTagCIADataset.argtypes = [c.c_void_p, c.c_long]

lib.setSpeciesTagLowerFrequency.restype = None
lib.setSpeciesTagLowerFrequency.argtypes = [c.c_void_p, c.c_double]

lib.setSpeciesTagUpperFrequency.restype = None
lib.setSpeciesTagUpperFrequency.argtypes = [c.c_void_p, c.c_double]

lib.validSpecies.restype = c.c_long
lib.validSpecies.argtypes = [c.c_long]

lib.validIsotopologue.restype = c.c_long
lib.validIsotopologue.argtypes = [c.c_long, c.c_long]
