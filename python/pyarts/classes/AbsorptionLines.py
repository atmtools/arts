import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.AbsorptionSingleLine import AbsorptionSingleLine
from pyarts.classes.LineShapeModel import LineShapeModel
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.QuantumNumbers import QuantumNumbers
from pyarts.classes.Rational import Rational
from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base, define_array_lib


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
                 quantumidentity=QuantumIdentifier(),
                 localquantumnumbers=[], broadeningspecies=[], lsm=LineShapeModel()):
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

            if not self.OK():
                raise ValueError("Bad initialization of class")

    def OK(self):
        """ Returns whether the class is OK or not """
        return bool(lib.isAbsorptionLinesOK(self.__data__))

    @property
    def selfbroadening(self):
        """ Self broadening (bool) """
        return lib.getAbsorptionLinesSelfBroadening(self.__data__)

    @selfbroadening.setter
    def selfbroadening(self, val):
        val = bool(val)
        lib.setAbsorptionLinesSelfBroadening(self.__data__, val)

    @property
    def bathbroadening(self):
        """ Air broadening (bool) """
        return lib.getAbsorptionLinesBathBroadening(self.__data__)

    @bathbroadening.setter
    def bathbroadening(self, val):
        val = bool(val)
        lib.setAbsorptionLinesBathBroadening(self.__data__, val)

    @property
    def cutoff(self):
        """ Type of cutoff (get: Index; set: str) """
        return lib.getAbsorptionLinesCutoffType(self.__data__)

    @cutoff.setter
    def cutoff(self, data):
        if isinstance(data, str):
            data = data.encode("ascii")
            if lib.setAbsorptionLinesCutoffType(self.__data__, data):
                raise ValueError("Invalid cutoff type")
        else:
            raise TypeError("Expects str input")

    @property
    def mirroring(self):
        """ Type of mirroring (get: Index; set: str) """
        return lib.getAbsorptionLinesMirroringType(self.__data__)

    @mirroring.setter
    def mirroring(self, data):
        if isinstance(data, str):
            data = data.encode("ascii")
            if lib.setAbsorptionLinesMirroringType(self.__data__, data):
                raise ValueError("Invalid mirroring type")
        else:
            raise TypeError("Expects str input")

    @property
    def population(self):
        """ Type of population (get: Index; set: str) """
        return lib.getAbsorptionLinesPopulationType(self.__data__)

    @population.setter
    def population(self, data):
        if isinstance(data, str):
            data = data.encode("ascii")
            if lib.setAbsorptionLinesPopulationType(self.__data__, data):
                raise ValueError("Invalid population type")
        else:
            raise TypeError("Expects str input")

    @property
    def normalization(self):
        """ Type of normalization (get: Index; set: str) """
        return lib.getAbsorptionLinesNormalizationType(self.__data__)

    @normalization.setter
    def normalization(self, data):
        if isinstance(data, str):
            data = str(data).encode("ascii")
            if lib.setAbsorptionLinesNormalizationType(self.__data__, data):
                raise ValueError("Invalid normalization type")
        else:
            raise TypeError("Expects str input")

    @property
    def lineshape(self):
        """ Type of lineshapetype (get: Index; set: str) """
        return lib.getAbsorptionLinesLineShapeType(self.__data__)

    @lineshape.setter
    def lineshape(self, data):
        if isinstance(data, str):
            data = str(data).encode("ascii")
            if lib.setAbsorptionLinesLineShapeType(self.__data__, data):
                raise ValueError("Invalid lineshape type")
        else:
            raise TypeError("Expects str input")

    @property
    def t0(self):
        """ Reference temperature (Numeric) """
        return lib.getAbsorptionLinesT0(self.__data__)

    @t0.setter
    def t0(self, data):
        data = float(data)
        lib.setAbsorptionLinesT0(self.__data__, data)

    @property
    def cutofffreq(self):
        """ Cutoff frequency (Numeric) """
        return lib.getAbsorptionLinesCutoffFrequency(self.__data__)

    @cutofffreq.setter
    def cutofffreq(self, data):
        data = float(data)
        lib.setAbsorptionLinesCutoffFrequency(self.__data__, data)

    @property
    def linemixinglimit(self):
        """ Line mixing pressure limit (Numeric) """
        return lib.getAbsorptionLinesLinemixingLimit(self.__data__)

    @linemixinglimit.setter
    def linemixinglimit(self, data):
        data = float(data)
        lib.setAbsorptionLinesLinemixingLimit(self.__data__, data)

    @property
    def quantumidentity(self):
        """ Quantum identity (QuantumIdentifier) """
        return QuantumIdentifier(c.c_void_p(lib.getAbsorptionLinesQuantumIdentifier(self.__data__)))

    @quantumidentity.setter
    def quantumidentity(self, val):
        self.quantumidentity.set(val)

    @property
    def sizelocalquantumnumbers(self):
        """ Number of local quantum numbers (Index) """
        return lib.getAbsorptionLinesLocalQuantumNumberCount(self.__data__)

    @sizelocalquantumnumbers.setter
    def sizelocalquantumnumbers(self, size):
        size = int(size)
        if size >= 0:
            lib.resizeAbsorptionLinesLocalQuantumNumber(size, self.__data__)
        else:
            raise ValueError("Invalid size")

    @property
    def localquantumnumbers(self):
        """ Local quantum numbers (list of Index) """
        n = self.sizelocalquantumnumbers
        x = []
        for i in range(n):
            x.append(lib.getAbsorptionLinesLocalQuantumNumber(i, self.__data__))
        return x

    @localquantumnumbers.setter
    def localquantumnumbers(self, val):
        if isinstance(val, Sized):
            self.sizelocalquantumnumbers = len(val)
            n = self.sizelocalquantumnumbers
            for i in range(n):
                lib.setAbsorptionLinesLocalQuantumNumber(i, self.__data__, QuantumNumbers.to_index(val[i]))
        else:
            raise TypeError("Only accepts array-like input")

    @property
    def sizebroadeningspecies(self):
        """ Number of broadening species (Index) """
        return lib.getAbsorptionLinesSpeciesTagCount(self.__data__)

    @sizebroadeningspecies.setter
    def sizebroadeningspecies(self, size):
        size = int(size)
        if size >= 0:
            lib.resizeAbsorptionLinesSpeciesTag(size, self.__data__)
        else:
            raise ValueError("Invalid size")

    @property
    def broadeningspecies(self):
        """ Broadening species (list of SpeciesTag) """
        n = self.sizebroadeningspecies
        x = []
        for i in range(n):
            x.append(SpeciesTag(c.c_void_p(lib.getAbsorptionLinesSpeciesTag(i, self.__data__))))
        return x

    @broadeningspecies.setter
    def broadeningspecies(self, val):
        if isinstance(val, Sized):
            self.sizebroadeningspecies = len(val)
            x = self.broadeningspecies
            n = self.sizebroadeningspecies
            for i in range(n):
                x[i].set(val[i])
        else:
            raise TypeError("Only accepts array-like input")

    @property
    def sizelines(self):
        """ Number of absorption lines (Index) """
        return lib.getAbsorptionLinesSingleLineCount(self.__data__)

    @sizelines.setter
    def sizelines(self, size):
        size = int(size)
        if size >= 0:
            lib.resizeAbsorptionLinesSingleLine(size, self.__data__)
        else:
            raise ValueError("Invalid size")

    @property
    def lines(self):
        """ Absorption lines (list of AbsorptionSingleLine) """
        n = self.sizelines
        x = []
        for i in range(n):
            x.append(AbsorptionSingleLine(c.c_void_p(lib.getAbsorptionLinesSingleLine(i, self.__data__))))
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
        lib.printAbsorptionLines(self.__data__)

    def printmeta(self):
        """ Print to cout the ARTS meta representation of the class

        Will check if the class is OK or raise a ValueError
        """
        if self.OK():
            lib.printmetaAbsorptionLines(self.__data__)
        else:
            raise ValueError("Class is in a bad state")

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionLines(self.__data__)

    def __repr__(self):
        return "ARTS AbsorptionLines"

    @staticmethod
    def name():
        """ Name of the class as a string.  Required for arrayification """
        return __class__.__name__

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, AbsorptionLines):
            self.selfbroadening = other.selfbroadening
            self.bathbroadening = other.bathbroadening
            lib.setAbsorptionLinesCutoffTypeByIndex(self.__data__, other.cutoff)
            lib.setAbsorptionLinesMirroringTypeByIndex(self.__data__, other.mirroring)
            lib.setAbsorptionLinesPopulationTypeByIndex(self.__data__, other.population)
            lib.setAbsorptionLinesNormalizationTypeByIndex(self.__data__, other.normalization)
            lib.setAbsorptionLinesLineShapeTypeByIndex(self.__data__, other.lineshapetype)
            self.t0 = other.t0
            self.cutofffreq = other.cutofffreq
            self.linemixinglimit = other.linemixinglimit
            self.quantumidentity = other.quantumidentity
            self.localquantumnumbers = other.localquantumnumbers
            self.broadeningspecies = other.broadeningspecies
            self.lsm = other.lsm
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
        if lib.xmlsaveAbsorptionLines(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


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

lib.getAbsorptionLinesSelfBroadening.restype = c.c_bool
lib.getAbsorptionLinesSelfBroadening.argtypes = [c.c_void_p]

lib.setAbsorptionLinesSelfBroadening.restype = None
lib.setAbsorptionLinesSelfBroadening.argtypes = [c.c_void_p, c.c_bool]

lib.getAbsorptionLinesCutoffType.restype = c.c_long
lib.getAbsorptionLinesCutoffType.argtypes = [c.c_void_p]

lib.setAbsorptionLinesCutoffType.restype = c.c_long
lib.setAbsorptionLinesCutoffType.argtypes = [c.c_void_p, c.c_char_p]

lib.setAbsorptionLinesCutoffTypeByIndex.restype = None
lib.setAbsorptionLinesCutoffTypeByIndex.argtypes = [c.c_void_p, c.c_long]

lib.getAbsorptionLinesMirroringType.restype = c.c_long
lib.getAbsorptionLinesMirroringType.argtypes = [c.c_void_p]

lib.setAbsorptionLinesMirroringType.restype = c.c_long
lib.setAbsorptionLinesMirroringType.argtypes = [c.c_void_p, c.c_char_p]

lib.setAbsorptionLinesMirroringTypeByIndex.restype = None
lib.setAbsorptionLinesMirroringTypeByIndex.argtypes = [c.c_void_p, c.c_long]

lib.getAbsorptionLinesPopulationType.restype = c.c_long
lib.getAbsorptionLinesPopulationType.argtypes = [c.c_void_p]

lib.setAbsorptionLinesPopulationType.restype = c.c_long
lib.setAbsorptionLinesPopulationType.argtypes = [c.c_void_p, c.c_char_p]

lib.setAbsorptionLinesPopulationTypeByIndex.restype = None
lib.setAbsorptionLinesPopulationTypeByIndex.argtypes = [c.c_void_p, c.c_long]

lib.getAbsorptionLinesNormalizationType.restype = c.c_long
lib.getAbsorptionLinesNormalizationType.argtypes = [c.c_void_p]

lib.setAbsorptionLinesNormalizationType.restype = c.c_long
lib.setAbsorptionLinesNormalizationType.argtypes = [c.c_void_p, c.c_char_p]

lib.setAbsorptionLinesNormalizationTypeByIndex.restype = None
lib.setAbsorptionLinesNormalizationTypeByIndex.argtypes = [c.c_void_p, c.c_long]

lib.getAbsorptionLinesLineShapeType.restype = c.c_long
lib.getAbsorptionLinesLineShapeType.argtypes = [c.c_void_p]

lib.setAbsorptionLinesLineShapeType.restype = c.c_long
lib.setAbsorptionLinesLineShapeType.argtypes = [c.c_void_p, c.c_char_p]

lib.setAbsorptionLinesLineShapeTypeByIndex.restype = None
lib.setAbsorptionLinesLineShapeTypeByIndex.argtypes = [c.c_void_p, c.c_long]

lib.getAbsorptionLinesT0.restype = c.c_double
lib.getAbsorptionLinesT0.argtypes = [c.c_void_p]

lib.setAbsorptionLinesT0.restype = None
lib.setAbsorptionLinesT0.argtypes = [c.c_void_p, c.c_double]

lib.getAbsorptionLinesCutoffFrequency.restype = c.c_double
lib.getAbsorptionLinesCutoffFrequency.argtypes = [c.c_void_p]

lib.setAbsorptionLinesCutoffFrequency.restype = None
lib.setAbsorptionLinesCutoffFrequency.argtypes = [c.c_void_p, c.c_double]

lib.getAbsorptionLinesLinemixingLimit.restype = c.c_double
lib.getAbsorptionLinesLinemixingLimit.argtypes = [c.c_void_p]

lib.setAbsorptionLinesLinemixingLimit.restype = None
lib.setAbsorptionLinesLinemixingLimit.argtypes = [c.c_void_p, c.c_double]

lib.getAbsorptionLinesQuantumIdentifier.restype = c.c_void_p
lib.getAbsorptionLinesQuantumIdentifier.argtypes = [c.c_void_p]

lib.resizeAbsorptionLinesLocalQuantumNumber.restype = None
lib.resizeAbsorptionLinesLocalQuantumNumber.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionLinesLocalQuantumNumber.restype = c.c_long
lib.getAbsorptionLinesLocalQuantumNumber.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionLinesLocalQuantumNumberCount.restype = c.c_long
lib.getAbsorptionLinesLocalQuantumNumberCount.argtypes = [c.c_void_p]

lib.setAbsorptionLinesLocalQuantumNumber.restype = None
lib.setAbsorptionLinesLocalQuantumNumber.argtypes = [c.c_long, c.c_void_p, c.c_long]

lib.resizeAbsorptionLinesSpeciesTag.restype = None
lib.resizeAbsorptionLinesSpeciesTag.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionLinesSpeciesTag.restype = c.c_void_p
lib.getAbsorptionLinesSpeciesTag.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionLinesSpeciesTagCount.restype = c.c_long
lib.getAbsorptionLinesSpeciesTagCount.argtypes = [c.c_void_p]

lib.resizeAbsorptionLinesSingleLine.restype = None
lib.resizeAbsorptionLinesSingleLine.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionLinesSingleLine.restype = c.c_void_p
lib.getAbsorptionLinesSingleLine.argtypes = [c.c_long, c.c_void_p]

lib.getAbsorptionLinesSingleLineCount.restype = c.c_long
lib.getAbsorptionLinesSingleLineCount.argtypes = [c.c_void_p]

lib.isAbsorptionLinesOK.restype = c.c_long
lib.isAbsorptionLinesOK.argtypes = [c.c_void_p]
