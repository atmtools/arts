import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class SpeciesTag:
    """ ARTS SpeciesTag data

    Properties:
        spec:
            Species of tag (Index)

        isot:
            Isotopologue of tag (Index)

        type:
            Type of tag (Index)

        cia_second:
            Secondary species when CIA mode (Index)

        cia_dataset:
            CIA dataset when CIA mode (Index)

        lf:
            Lower frequency range (Numeric)

        uf:
            Upper frequency range (Numeric)
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
            data = data.encode("ascii")
            if lib.setSpeciesTag(self.__data__, data):
                raise ValueError("Invalid species tag")
        else:
            raise TypeError("Expects str input")

    @property
    def spec(self):
        """ Species of tag (Index) """
        return lib.getSpeciesSpeciesTag(self.__data__)

    @spec.setter
    def spec(self, spec):
        if self.validSpecies(spec):
            return lib.setSpeciesSpeciesTag(self.__data__, int(spec))
        else:
            raise ValueError("Invalid species")

    @property
    def isot(self):
        """ Isotopologue of tag (Index) """
        return lib.getIsotopologueSpeciesTag(self.__data__)

    @isot.setter
    def isot(self, isot):
        if self.validIsotopologue(self.spec, isot) or self.validContinuum(self.spec, isot) or self.validAllIsotopologues(self.spec, isot):
            return lib.setIsotopologueSpeciesTag(self.__data__, int(isot))
        else:
            raise ValueError("Invalid isotopologue")

    @property
    def type(self):
        """ Type of tag (Index) """
        return lib.getTypeSpeciesTag(self.__data__)

    @type.setter
    def type(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexTypeSpeciesTag(self.__data__, type.encode("ascii")))
        else:
            type = int(type)
            if lib.setTypeSpeciesTag(self.__data__, type):
                raise ValueError("Invalid type")

    @property
    def cia_second(self):
        """ Secondary species when CIA mode (Index) """
        return lib.getCIASecondSpeciesTag(self.__data__)

    @cia_second.setter
    def cia_second(self, spec):
        if self.validSpecies(spec):
            lib.setCIASecondSpeciesTag(self.__data__, int(spec))
        elif self.cia_type != self.type:
            lib.setCIASecondSpeciesTag(self.__data__, int(spec))
        else:
            raise ValueError("Invalid species")

    @property
    def cia_dataset(self):
        """ CIA dataset when CIA mode (Index) """
        return lib.getCIADatasetSpeciesTag(self.__data__)

    @cia_dataset.setter
    def cia_dataset(self, data):
        lib.setCIADatasetSpeciesTag(self.__data__, int(data))

    @property
    def lf(self):
        """ Lower frequency range (Numeric) """
        return lib.getLfSpeciesTag(self.__data__)

    @lf.setter
    def lf(self, x):
        x = float(x)
        if x <= self.uf:
            lib.setLfSpeciesTag(self.__data__, x)
        else:
            raise ValueError("Too high cf upper frequency")

    @property
    def uf(self):
        """ Upper frequency range (Numeric) """
        return lib.getLfSpeciesTag(self.__data__)

    @uf.setter
    def uf(self, x):
        x = float(x)
        if x >= self.lf:
            lib.setUfSpeciesTag(self.__data__, x)
        else:
            raise ValueError("Too low cf lower frequency")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        if self.OK:
            lib.printSpeciesTag(self.__data__)
        else:
            raise RuntimeError("Cannot print SpeciesTag in bad state")

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesTag(self.__data__)

    def __repr__(self):
        return "ARTS SpeciesTag"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, SpeciesTag):
            self.type = other.type
            self.spec = other.spec
            self.isot = other.isot
            self.cia_second = other.cia_second
            self.cia_dataset = other.cia_dataset
            self.uf = float("inf")
            self.lf = other.lf
            self.uf = other.uf
        else:
            raise TypeError("Expects SpeciesTag")

    @property
    def OK(self):
        """ Returns true if the class is OK """
        return self.validCIASpecies() and (self.validIsotopologue(self.spec, self.isot) or self.validContinuum(self.spec, self.isot) or self.validAllIsotopologues(self.spec, self.isot))

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

    @property
    def cia_type(self):
        """ Returns the CIA type """
        return lib.string2indexTypeSpeciesTag(self.__data__, "TYPE_CIA".encode("ascii"))

    def validCIASpecies(self):
        """ Returns whether cia species is a valid species according to ARTS

        Input:
            spec:
                Species (Index)

        Output:
            Boolean True or False
        """
        if self.cia_type != self.type or lib.validSpecies(self.cia_second) == 0:
            return True
        else:
            return False

    @staticmethod
    def validAllIsotopologues(spec, isot):
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
            if lib.validAllIsotopologues(spec, isot) == 0:
                return True
            else:
                return False
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

    @staticmethod
    def validContinuum(spec, isot):
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
            if lib.validContinuum(spec, isot) == 0:
                return True
            else:
                return False
        else:
            return False

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
                self.type == other.type and \
                self.spec == other.spec and \
                self.isot == other.isot and \
                self.cia_second == other.cia_second and \
                self.cia_dataset == other.cia_dataset and \
                self.lf == other.lf and \
                self.uf == other.uf:
            return True
        else:
            return False


# Generate ArrayOfSpeciesTag
exec(array_base(SpeciesTag))


# Generate ArrayOfArrayOfSpeciesTag
exec(array_base(ArrayOfSpeciesTag))


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

lib.getSpeciesSpeciesTag.restype = c.c_long
lib.getSpeciesSpeciesTag.argtypes = [c.c_void_p]

lib.getIsotopologueSpeciesTag.restype = c.c_long
lib.getIsotopologueSpeciesTag.argtypes = [c.c_void_p]

lib.getTypeSpeciesTag.restype = c.c_long
lib.getTypeSpeciesTag.argtypes = [c.c_void_p]

lib.string2indexTypeSpeciesTag.restype = c.c_long
lib.string2indexTypeSpeciesTag.argtypes = [c.c_void_p, c.c_char_p]

lib.getCIASecondSpeciesTag.restype = c.c_long
lib.getCIASecondSpeciesTag.argtypes = [c.c_void_p]

lib.getCIADatasetSpeciesTag.restype = c.c_long
lib.getCIADatasetSpeciesTag.argtypes = [c.c_void_p]

lib.getLfSpeciesTag.restype = c.c_double
lib.getLfSpeciesTag.argtypes = [c.c_void_p]

lib.getUfSpeciesTag.restype = c.c_double
lib.getUfSpeciesTag.argtypes = [c.c_void_p]

lib.setSpeciesSpeciesTag.restype = None
lib.setSpeciesSpeciesTag.argtypes = [c.c_void_p, c.c_long]

lib.setIsotopologueSpeciesTag.restype = None
lib.setIsotopologueSpeciesTag.argtypes = [c.c_void_p, c.c_long]

lib.setTypeSpeciesTag.restype = c.c_long
lib.setTypeSpeciesTag.argtypes = [c.c_void_p, c.c_long]

lib.setCIASecondSpeciesTag.restype = None
lib.setCIASecondSpeciesTag.argtypes = [c.c_void_p, c.c_long]

lib.setCIADatasetSpeciesTag.restype = None
lib.setCIADatasetSpeciesTag.argtypes = [c.c_void_p, c.c_long]

lib.setLfSpeciesTag.restype = None
lib.setLfSpeciesTag.argtypes = [c.c_void_p, c.c_double]

lib.setUfSpeciesTag.restype = None
lib.setUfSpeciesTag.argtypes = [c.c_void_p, c.c_double]

lib.validSpecies.restype = c.c_long
lib.validSpecies.argtypes = [c.c_long]

lib.validIsotopologue.restype = c.c_long
lib.validIsotopologue.argtypes = [c.c_long, c.c_long]

lib.validContinuum.restype = c.c_long
lib.validContinuum.argtypes = [c.c_long, c.c_long]

lib.validAllIsotopologues.restype = c.c_long
lib.validAllIsotopologues.argtypes = [c.c_long, c.c_long]
