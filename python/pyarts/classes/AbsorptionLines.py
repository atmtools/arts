import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.AbsorptionSingleLine import AbsorptionSingleLine
from pyarts.classes.LineShapeModel import LineShapeModel
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.QuantumNumbers import QuantumNumbers
from pyarts.classes.Rational import Rational
from pyarts.classes.SpeciesTag import ArrayOfSpeciesTag
from pyarts.classes.BasicTypes import Index, String
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class AbsorptionLines:
    """ ARTS AbsorptionLines data

    Type variables can only be set from strings or by copy from another
    AbsorptionLines-instance

    Properties:
        selfbroadening:
            Self broadening (bool)

        bathbroadening:
            Air broadening (bool)

        cutoff:
            Type of cutoff (get: Index; set: str)

        mirroring:
            Type of mirroring (get: Index; set: str)

        population:
            Type of population (get: Index; set: str)

        normalization:
            Type of normalization (get: Index; set: str)

        lineshapetype:
            Type of lineshapetype (get: Index; set: str)

        t0:
            Reference temperature (Numeric)

        cutofffreq:
            Cutoff frequency (Numeric)

        linemixinglimit:
            Line mixing pressure limit (Numeric)

        quantumidentity:
            Quantum identity (QuantumIdentifier)

        sizelocalquantumnumbers:
            Number of local quantum numbers (Index)

        localquantumnumbers:
            Local quantum numbers (list of Index)

        sizebroadeningspecies:
            Number of broadening species (Index)

        broadeningspecies:
            Broadening species (list of SpeciesTag)

        sizelines:
            Number of absorption lines (Index)

        lines:
            Absorption lines (list of AbsorptionSingleLine)
    """
    def __init__(self, selfbroadening=False, bathbroadening=False, nlines=0,
                 cutoff="None", mirroring="None", population="LTE", normalization="None",
                 lineshapetype="VP", t0=296, cutofffreq=-1, linemixinglimit=-1,
                 quantumidentity=QuantumIdentifier(), localquantumnumbers=[],
                 broadeningspecies=ArrayOfSpeciesTag(), lsm=LineShapeModel()):
        if isinstance(selfbroadening, c.c_void_p):
            self.__delete__ = False
            self.__data__ = selfbroadening
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createAbsorptionLines())
            self.selfbroadening = selfbroadening
            self.bathbroadening = bathbroadening
            self.cutoff = cutoff
            self.mirroring = mirroring
            self.population = population
            self.normalization = normalization
            self.lineshapetype = lineshapetype
            self.t0 = t0
            self.cutofffreq = cutofffreq
            self.linemixinglimit = linemixinglimit
            self.quantumidentity = quantumidentity
            self.localquantumnumbers = localquantumnumbers
            self.broadeningspecies = broadeningspecies

            n = self.sizelocalquantumnumbers
            x = AbsorptionSingleLine(lsm=lsm, qupp=[Rational()]*n, qlow=[Rational()]*n)
            self.sizelines = nlines
            y = self.lines
            for line in y:
                line.set(x)

            if not self.OK:
                raise ValueError("Bad initialization of class")

    @property
    def OK(self):
        """ Returns whether the class is OK or not """
        return bool(lib.isAbsorptionLinesOK(self.__data__))

    @property
    def species_name(self):
        """ Return the name of the species """
        if self.OK:
            species = c.c_void_p(lib.getSpeciesNameAbsorptionLines(self.__data__))
            return String(species, delete=True).val
        else:
            raise RuntimeError("Cannot access SpeciesName; class is not OK")

    @property
    def selfbroadening(self):
        """ Self broadening (bool) """
        return lib.getSelfAbsorptionLines(self.__data__)

    @selfbroadening.setter
    def selfbroadening(self, val):
        lib.setSelfAbsorptionLines(self.__data__, bool(val))

    @property
    def bathbroadening(self):
        """ Air broadening (bool) """
        return lib.getBathAbsorptionLines(self.__data__)

    @bathbroadening.setter
    def bathbroadening(self, val):
        lib.setBathAbsorptionLines(self.__data__, bool(val))

    @property
    def cutoff(self):
        """ Type of cutoff (get: Index; set: str) """
        return lib.getCutoffAbsorptionLines(self.__data__)

    @cutoff.setter
    def cutoff(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexCutoffAbsorptionLines(self.__data__, type.encode("ascii")))
        else:
            if lib.setCutoffAbsorptionLines(self.__data__, int(type)):
                raise ValueError("Invalid type")

    @property
    def mirroring(self):
        """ Type of mirroring (get: Index; set: str) """
        return lib.getMirroringAbsorptionLines(self.__data__)

    @mirroring.setter
    def mirroring(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexMirroringAbsorptionLines(self.__data__, type.encode("ascii")))
        else:
            if lib.setMirroringAbsorptionLines(self.__data__, int(type)):
                raise ValueError("Invalid type")

    @property
    def population(self):
        """ Type of population (get: Index; set: str) """
        return lib.getPopulationAbsorptionLines(self.__data__)

    @population.setter
    def population(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexPopulationAbsorptionLines(self.__data__, type.encode("ascii")))
        else:
            if lib.setPopulationAbsorptionLines(self.__data__, int(type)):
                raise ValueError("Invalid type")

    @property
    def normalization(self):
        """ Type of normalization (get: Index; set: str) """
        return lib.getNormalizationAbsorptionLines(self.__data__)

    @normalization.setter
    def normalization(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexPopulationAbsorptionLines(self.__data__, type.encode("ascii")))
        else:
            if lib.setNormalizationAbsorptionLines(self.__data__, int(type)):
                raise ValueError("Invalid type")

    @property
    def lineshapetype(self):
        """ Type of lineshapetype (get: Index; set: str) """
        return lib.getLineShapeTypeAbsorptionLines(self.__data__)

    @lineshapetype.setter
    def lineshapetype(self, type):
        if isinstance(type, str):
            self.type = int(lib.string2indexLineShapeTypeAbsorptionLines(self.__data__, type.encode("ascii")))
        else:
            if lib.setLineShapeTypeAbsorptionLines(self.__data__, int(type)):
                raise ValueError("Invalid type")

    @property
    def t0(self):
        """ Reference temperature (Numeric) """
        return lib.getT0AbsorptionLines(self.__data__)

    @t0.setter
    def t0(self, data):
        lib.setT0AbsorptionLines(self.__data__, float(data))

    @property
    def cutofffreq(self):
        """ Cutoff frequency (Numeric) """
        return lib.getCutoffFreqValueAbsorptionLines(self.__data__)

    @cutofffreq.setter
    def cutofffreq(self, data):
        lib.setCutoffFreqValueAbsorptionLines(self.__data__, float(data))

    @property
    def linemixinglimit(self):
        """ Line mixing pressure limit (Numeric) """
        return lib.getLinemixingLimitAbsorptionLines(self.__data__)

    @linemixinglimit.setter
    def linemixinglimit(self, data):
        lib.setLinemixingLimitAbsorptionLines(self.__data__, float(data))

    @property
    def quantumidentity(self):
        """ Quantum identity (QuantumIdentifier) """
        return QuantumIdentifier(c.c_void_p(lib.getQuantumIdentityAbsorptionLines(self.__data__)))

    @quantumidentity.setter
    def quantumidentity(self, val):
        self.quantumidentity.set(val)

    @property
    def sizelocalquantumnumbers(self):
        """ Number of local quantum numbers (Index) """
        return lib.sizeLocalQuantaAbsorptionLines(self.__data__)

    @sizelocalquantumnumbers.setter
    def sizelocalquantumnumbers(self, size):
        size = int(size)
        if size >= 0:
            lib.resizeLocalQuantaAbsorptionLines(size, self.__data__)
        else:
            raise ValueError("Invalid size")

    @property
    def localquantumnumbers(self):
        """ Local quantum numbers (list of Index) """
        n = self.sizelocalquantumnumbers
        x = []
        for i in range(n):
            x.append(Index(QuantumNumbers.to_index(
                    lib.getLocalQuantaAbsorptionLines(i, self.__data__))))
        return x

    @localquantumnumbers.setter
    def localquantumnumbers(self, val):
        if isinstance(val, Sized):
            self.sizelocalquantumnumbers = len(val)
            n = self.sizelocalquantumnumbers
            for i in range(n):
                lib.setLocalQuantaAbsorptionLines(i, self.__data__,
                                                  QuantumNumbers.to_index(val[i]))
        else:
            raise TypeError("Only accepts array-like input")

    @property
    def broadeningspecies(self):
        """ Broadening species (ArrayOfSpeciesTag) """
        return ArrayOfSpeciesTag(c.c_void_p(
                lib.getBroadeningSpeciesAbsorptionLines(self.__data__)))

    @broadeningspecies.setter
    def broadeningspecies(self, val):
        self.broadeningspecies.set(val)

    @property
    def sizelines(self):
        """ Number of absorption lines (Index) """
        return lib.sizeAllLinesAbsorptionLines(self.__data__)

    @sizelines.setter
    def sizelines(self, size):
        size = int(size)
        if size >= 0:
            lib.resizeAllLinesAbsorptionLines(size, self.__data__)
        else:
            raise ValueError("Invalid size")

    @property
    def lines(self):
        """ Absorption lines (list of AbsorptionSingleLine) """
        n = self.sizelines
        x = []
        for i in range(n):
            x.append(AbsorptionSingleLine(c.c_void_p(
                    lib.getelemAllLinesAbsorptionLines(i, self.__data__))))
        return x

    @lines.setter
    def lines(self, val):
        if isinstance(val, Sized):
            self.sizelines = len(val)
            x = self.lines
            n = self.sizelines
            for i in range(n):
                x[i].set(val[i])
        else:
            raise TypeError("Only accepts array-like input")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        if self.OK:
            lib.printAbsorptionLines(self.__data__)
        else:
            raise ValueError("Cannot print; class is not OK")

    def printmeta(self):
        """ Print to cout the ARTS meta representation of the class

        Will check if the class is OK or raise a ValueError
        """
        if self.OK:
            lib.printmetaAbsorptionLines(self.__data__)
        else:
            raise ValueError("Cannot understand metadata; class is not OK")

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionLines(self.__data__)

    def __repr__(self):
        return "ARTS AbsorptionLines"

    @staticmethod
    def name():
        """ Name of the class as a string.  Required for arrayification """
        return "AbsorptionLines"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, AbsorptionLines):
            self.selfbroadening = other.selfbroadening
            self.bathbroadening = other.bathbroadening
            self.cutoff = other.cutoff
            self.mirroring = other.mirroring
            self.population = other.population
            self.normalization = other.normalization
            self.lineshapetype = other.lineshapetype
            self.t0 = other.t0
            self.cutofffreq = other.cutofffreq
            self.linemixinglimit = other.linemixinglimit
            self.quantumidentity = other.quantumidentity
            self.localquantumnumbers = other.localquantumnumbers
            self.broadeningspecies = other.broadeningspecies
            self.lines = other.lines
        else:
            raise TypeError("Expects AbsorptionLines")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadAbsorptionLines(self.__data__, correct_read_arguments(file)):
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
            raise ValueError("Cannot store; class is not OK")

        if lib.xmlsaveAbsorptionLines(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, AbsorptionLines) and \
                self.selfbroadening == other.selfbroadening and \
                self.bathbroadening == other.bathbroadening and \
                self.cutoff == other.cutoff and \
                self.mirroring == other.mirroring and \
                self.population == other.population and \
                self.normalization == other.normalization and \
                self.lineshapetype == other.lineshapetype and \
                self.t0 == other.t0 and \
                self.cutofffreq == other.cutofffreq and \
                self.linemixinglimit == other.linemixinglimit and \
                self.quantumidentity == other.quantumidentity and \
                self.localquantumnumbers == other.localquantumnumbers and \
                self.broadeningspecies == other.broadeningspecies and \
                self.lines == other.lines:
            return True
        else:
            return False


# Generate ArrayOfAbsorptionLines
exec(array_base(AbsorptionLines))


# Generate ArrayOfArrayOfAbsorptionLines
exec(array_base(ArrayOfAbsorptionLines))


lib.createAbsorptionLines.restype = c.c_void_p
lib.createAbsorptionLines.argtypes = []

lib.deleteAbsorptionLines.restype = None
lib.deleteAbsorptionLines.argtypes = [c.c_void_p]

lib.printAbsorptionLines.restype = None
lib.printAbsorptionLines.argtypes = [c.c_void_p]

lib.xmlreadAbsorptionLines.restype = c.c_long
lib.xmlreadAbsorptionLines.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveAbsorptionLines.restype = c.c_long
lib.xmlsaveAbsorptionLines.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.printmetaAbsorptionLines.restype = None
lib.printmetaAbsorptionLines.argtypes = [c.c_void_p]

lib.getSelfAbsorptionLines.restype = c.c_bool
lib.getSelfAbsorptionLines.argtypes = [c.c_void_p]

lib.setSelfAbsorptionLines.restype = None
lib.setSelfAbsorptionLines.argtypes = [c.c_void_p, c.c_bool]

lib.getBathAbsorptionLines.restype = c.c_bool
lib.getBathAbsorptionLines.argtypes = [c.c_void_p]

lib.setBathAbsorptionLines.restype = None
lib.setBathAbsorptionLines.argtypes = [c.c_void_p, c.c_bool]

lib.getCutoffAbsorptionLines.restype = c.c_long
lib.getCutoffAbsorptionLines.argtypes = [c.c_void_p]

lib.setCutoffAbsorptionLines.restype = c.c_long
lib.setCutoffAbsorptionLines.argtypes = [c.c_void_p, c.c_long]

lib.string2indexCutoffAbsorptionLines.restype = c.c_long
lib.string2indexCutoffAbsorptionLines.argtypes = [c.c_void_p, c.c_char_p]

lib.getMirroringAbsorptionLines.restype = c.c_long
lib.getMirroringAbsorptionLines.argtypes = [c.c_void_p]

lib.setMirroringAbsorptionLines.restype = c.c_long
lib.setMirroringAbsorptionLines.argtypes = [c.c_void_p, c.c_long]

lib.string2indexMirroringAbsorptionLines.restype = c.c_long
lib.string2indexMirroringAbsorptionLines.argtypes = [c.c_void_p, c.c_char_p]

lib.getLineShapeTypeAbsorptionLines.restype = c.c_long
lib.getLineShapeTypeAbsorptionLines.argtypes = [c.c_void_p]

lib.setLineShapeTypeAbsorptionLines.restype = c.c_long
lib.setLineShapeTypeAbsorptionLines.argtypes = [c.c_void_p, c.c_long]

lib.string2indexLineShapeTypeAbsorptionLines.restype = c.c_long
lib.string2indexLineShapeTypeAbsorptionLines.argtypes = [c.c_void_p, c.c_char_p]

lib.getPopulationAbsorptionLines.restype = c.c_long
lib.getPopulationAbsorptionLines.argtypes = [c.c_void_p]

lib.setPopulationAbsorptionLines.restype = c.c_long
lib.setPopulationAbsorptionLines.argtypes = [c.c_void_p, c.c_long]

lib.string2indexPopulationAbsorptionLines.restype = c.c_long
lib.string2indexPopulationAbsorptionLines.argtypes = [c.c_void_p, c.c_char_p]

lib.getNormalizationAbsorptionLines.restype = c.c_long
lib.getNormalizationAbsorptionLines.argtypes = [c.c_void_p]

lib.setNormalizationAbsorptionLines.restype = c.c_long
lib.setNormalizationAbsorptionLines.argtypes = [c.c_void_p, c.c_long]

lib.string2indexNormalizationAbsorptionLines.restype = c.c_long
lib.string2indexPopulationAbsorptionLines.argtypes = [c.c_void_p, c.c_char_p]

lib.getT0AbsorptionLines.restype = c.c_double
lib.getT0AbsorptionLines.argtypes = [c.c_void_p]

lib.setT0AbsorptionLines.restype = None
lib.setT0AbsorptionLines.argtypes = [c.c_void_p, c.c_double]

lib.getCutoffFreqValueAbsorptionLines.restype = c.c_double
lib.getCutoffFreqValueAbsorptionLines.argtypes = [c.c_void_p]

lib.setCutoffFreqValueAbsorptionLines.restype = None
lib.setCutoffFreqValueAbsorptionLines.argtypes = [c.c_void_p, c.c_double]

lib.getLinemixingLimitAbsorptionLines.restype = c.c_double
lib.getLinemixingLimitAbsorptionLines.argtypes = [c.c_void_p]

lib.setLinemixingLimitAbsorptionLines.restype = None
lib.setLinemixingLimitAbsorptionLines.argtypes = [c.c_void_p, c.c_double]

lib.getQuantumIdentityAbsorptionLines.restype = c.c_void_p
lib.getQuantumIdentityAbsorptionLines.argtypes = [c.c_void_p]

lib.sizeLocalQuantaAbsorptionLines.restype = c.c_long
lib.sizeLocalQuantaAbsorptionLines.argtypes = [c.c_void_p]

lib.resizeLocalQuantaAbsorptionLines.restype = None
lib.resizeLocalQuantaAbsorptionLines.argtypes = [c.c_long, c.c_void_p]

lib.getLocalQuantaAbsorptionLines.restype = c.c_long
lib.getLocalQuantaAbsorptionLines.argtypes = [c.c_long, c.c_void_p]

lib.setLocalQuantaAbsorptionLines.restype = None
lib.setLocalQuantaAbsorptionLines.argtypes = [c.c_long, c.c_void_p, c.c_long]

lib.getBroadeningSpeciesAbsorptionLines.restype = c.c_void_p
lib.getBroadeningSpeciesAbsorptionLines.argtypes = [c.c_void_p]

lib.sizeAllLinesAbsorptionLines.restype = c.c_long
lib.sizeAllLinesAbsorptionLines.argtypes = [c.c_void_p]

lib.resizeAllLinesAbsorptionLines.restype = None
lib.resizeAllLinesAbsorptionLines.argtypes = [c.c_long, c.c_void_p]

lib.getelemAllLinesAbsorptionLines.restype = c.c_void_p
lib.getelemAllLinesAbsorptionLines.argtypes = [c.c_long, c.c_void_p]

lib.isAbsorptionLinesOK.restype = c.c_long
lib.isAbsorptionLinesOK.argtypes = [c.c_void_p]

lib.getSpeciesNameAbsorptionLines.restype = c.c_void_p
lib.getSpeciesNameAbsorptionLines.argtypes = [c.c_void_p]
