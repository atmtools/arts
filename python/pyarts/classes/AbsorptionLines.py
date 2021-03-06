import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.AbsorptionSingleLine import AbsorptionSingleLine
from pyarts.classes.LineShapeModel import LineShapeModel, LineShapeType
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.QuantumNumbers import QuantumNumberType
from pyarts.classes.Rational import Rational
from pyarts.classes.BasicTypes import String
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base
from pyarts.classes.SpeciesIsotopeRecord import ArrayOfSpecies


class VectorOfQuantumNumberType:
    def __init__(self, data):
        assert isinstance(data, c.c_void_p)
        self.__data__ = data
    
    def __getitem__(self, key):
        x = lib.getQuantumNumberTypeLocalQuantaAbsorptionLines(self.__data__, int(key))
        if x:
            return QuantumNumberType(c.c_void_p(x))
        else:
            raise IndexError("Out of bounds")
    
    def __setitem__(self, key, val):
        self.__getitem__(key).set(val)
    
    def __len__(self):
        return int(lib.sizeLocalQuantaAbsorptionLines(self.__data__))
    
    def resize(self, n):
        n = int(n)
        assert n >= 0
        lib.resizeLocalQuantaAbsorptionLines(n, self.__data__)
    
    def __repr__(self):
        out = "["
        for i in range(len(self)):
            out += str(self[i])
            if (i+1) != len(self):
                out += ', '
        return out+ ']'
    
    def set(self, other):
        self.resize(len(other))
        for i in range(len(self)):
            self[i] = other[i]
    
    def __eq__(self, other):
        if not len(self) == len(other):
            return False
        else:
            for i in range(len(self)):
                if not self[i] == other[i]:
                    return False
        return True
            

lib.sizeLocalQuantaAbsorptionLines.restype = c.c_long
lib.sizeLocalQuantaAbsorptionLines.argtypes = [c.c_void_p]

lib.resizeLocalQuantaAbsorptionLines.restype = None
lib.resizeLocalQuantaAbsorptionLines.argtypes = [c.c_long, c.c_void_p]

lib.getQuantumNumberTypeLocalQuantaAbsorptionLines.restype = c.c_void_p
lib.getQuantumNumberTypeLocalQuantaAbsorptionLines.argtypes = [c.c_void_p, c.c_long]


class AbsorptionPopulationType:
    """ ARTS AbsorptionPopulationType type data
    
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
            self.__data__ = c.c_void_p(lib.createAbsorptionPopulationType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getAbsorptionPopulationTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setAbsorptionPopulationTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad AbsorptionPopulationType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAbsorptionPopulationType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionPopulationType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, AbsorptionPopulationType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, AbsorptionPopulationType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createAbsorptionPopulationType.restype = c.c_void_p
lib.createAbsorptionPopulationType.argtypes = []

lib.deleteAbsorptionPopulationType.restype = None
lib.deleteAbsorptionPopulationType.argtypes = [c.c_void_p]

lib.printAbsorptionPopulationType.restype = None
lib.printAbsorptionPopulationType.argtypes = [c.c_void_p]

lib.getAbsorptionPopulationTypeString.restype = c.c_void_p
lib.getAbsorptionPopulationTypeString.argtypes = [c.c_void_p]

lib.setAbsorptionPopulationTypeString.restype = c.c_int
lib.setAbsorptionPopulationTypeString.argtypes = [c.c_void_p, c.c_char_p]


class AbsorptionNormalizationType:
    """ ARTS AbsorptionNormalizationType data
    
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
            self.__data__ = c.c_void_p(lib.createAbsorptionNormalizationType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getAbsorptionNormalizationTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setAbsorptionNormalizationTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad AbsorptionNormalizationType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAbsorptionNormalizationType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionNormalizationType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, AbsorptionNormalizationType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, AbsorptionNormalizationType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createAbsorptionNormalizationType.restype = c.c_void_p
lib.createAbsorptionNormalizationType.argtypes = []

lib.deleteAbsorptionNormalizationType.restype = None
lib.deleteAbsorptionNormalizationType.argtypes = [c.c_void_p]

lib.printAbsorptionNormalizationType.restype = None
lib.printAbsorptionNormalizationType.argtypes = [c.c_void_p]

lib.getAbsorptionNormalizationTypeString.restype = c.c_void_p
lib.getAbsorptionNormalizationTypeString.argtypes = [c.c_void_p]

lib.setAbsorptionNormalizationTypeString.restype = c.c_int
lib.setAbsorptionNormalizationTypeString.argtypes = [c.c_void_p, c.c_char_p]


class AbsorptionMirroringType:
    """ ARTS AbsorptionMirroringType data
    
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
            self.__data__ = c.c_void_p(lib.createAbsorptionMirroringType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getAbsorptionMirroringTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setAbsorptionMirroringTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad AbsorptionMirroringType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAbsorptionMirroringType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionMirroringType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, AbsorptionMirroringType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, AbsorptionMirroringType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createAbsorptionMirroringType.restype = c.c_void_p
lib.createAbsorptionMirroringType.argtypes = []

lib.deleteAbsorptionMirroringType.restype = None
lib.deleteAbsorptionMirroringType.argtypes = [c.c_void_p]

lib.printAbsorptionMirroringType.restype = None
lib.printAbsorptionMirroringType.argtypes = [c.c_void_p]

lib.getAbsorptionMirroringTypeString.restype = c.c_void_p
lib.getAbsorptionMirroringTypeString.argtypes = [c.c_void_p]

lib.setAbsorptionMirroringTypeString.restype = c.c_int
lib.setAbsorptionMirroringTypeString.argtypes = [c.c_void_p, c.c_char_p]


class AbsorptionCutoffType:
    """ ARTS AbsorptionCutoffType data
    
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
            self.__data__ = c.c_void_p(lib.createAbsorptionCutoffType())
            self.name = data
    
    @property
    def name(self):
        return String(c.c_void_p(lib.getAbsorptionCutoffTypeString(self.__data__)),
                      delete=True)
    
    @name.setter
    def name(self, x):
        if lib.setAbsorptionCutoffTypeString(self.__data__, str(x).encode('utf-8')):
            raise RuntimeError(f"Bad AbsorptionCutoffType: {x}")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printAbsorptionCutoffType(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteAbsorptionCutoffType(self.__data__)
            
    def set(self, other):
        s = other.name if isinstance(other, AbsorptionCutoffType) else other
        self.name = s
    
    def __eq__(self, other):
        s = other.name if isinstance(other, AbsorptionCutoffType) else other
        return self.name == s
    
    def __repr__(self):
        return f"{self.name}"


lib.createAbsorptionCutoffType.restype = c.c_void_p
lib.createAbsorptionCutoffType.argtypes = []

lib.deleteAbsorptionCutoffType.restype = None
lib.deleteAbsorptionCutoffType.argtypes = [c.c_void_p]

lib.printAbsorptionCutoffType.restype = None
lib.printAbsorptionCutoffType.argtypes = [c.c_void_p]

lib.getAbsorptionCutoffTypeString.restype = c.c_void_p
lib.getAbsorptionCutoffTypeString.argtypes = [c.c_void_p]

lib.setAbsorptionCutoffTypeString.restype = c.c_int
lib.setAbsorptionCutoffTypeString.argtypes = [c.c_void_p, c.c_char_p]


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
            Type of cutoff (str)

        mirroring:
            Type of mirroring (str)

        population:
            Type of population (str)

        normalization:
            Type of normalization (str)

        lineshapetype:
            Type of lineshapetype (str)

        t0:
            Reference temperature (Numeric)

        cutofffreq:
            Cutoff frequency (Numeric)

        linemixinglimit:
            Line mixing pressure limit (Numeric)

        quantumidentity:
            Quantum identity (QuantumIdentifier)

        localquantumnumbers:
            Local quantum numbers (list of QuantumNumberType(s))

        sizebroadeningspecies:
            Number of broadening species (Index)

        broadeningspecies:
            Broadening species (list of Species)

        sizelines:
            Number of absorption lines (Index)

        lines:
            Absorption lines (list of AbsorptionSingleLine)
    """
    def __init__(self, selfbroadening=False, bathbroadening=False, nlines=0,
                 cutoff="None", mirroring="None", population="LTE", normalization="None",
                 lineshapetype="VP", t0=296, cutofffreq=-1, linemixinglimit=-1,
                 quantumidentity=QuantumIdentifier(), localquantumnumbers=[],
                 broadeningspecies=ArrayOfSpecies(), lsm=LineShapeModel()):
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

            n = len(self.localquantumnumbers)
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
        """ Type of cutoff """
        return AbsorptionCutoffType(c.c_void_p(lib.getCutoffAbsorptionLines(self.__data__)))

    @cutoff.setter
    def cutoff(self, type):
        self.cutoff.set(type)

    @property
    def mirroring(self):
        """ Type of mirroring """
        return AbsorptionMirroringType(c.c_void_p(lib.getMirroringAbsorptionLines(self.__data__)))

    @mirroring.setter
    def mirroring(self, type):
        self.mirroring.set(type)

    @property
    def population(self):
        """ Type of population """
        return AbsorptionPopulationType(c.c_void_p(lib.getPopulationAbsorptionLines(self.__data__)))

    @population.setter
    def population(self, type):
        self.population.set(type)

    @property
    def normalization(self):
        """ Type of normalization """
        return AbsorptionNormalizationType(c.c_void_p(lib.getNormalizationAbsorptionLines(self.__data__)))

    @normalization.setter
    def normalization(self, type):
        self.normalization.set(type)

    @property
    def lineshapetype(self):
        """ Type of lineshapetype """
        return LineShapeType(c.c_void_p(lib.getLineShapeTypeAbsorptionLines(self.__data__)))

    @lineshapetype.setter
    def lineshapetype(self, type):
        self.lineshapetype.set(type)

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
    def localquantumnumbers(self):
        """ Local quantum numbers (list of QuantumNumberType(s)) """
        return VectorOfQuantumNumberType(c.c_void_p(lib.getLocalQuantaAbsorptionLines(self.__data__)))

    @localquantumnumbers.setter
    def localquantumnumbers(self, val):
        self.localquantumnumbers.set(val)

    @property
    def broadeningspecies(self):
        """ Broadening species (ArrayOfSpecies) """
        return ArrayOfSpecies(c.c_void_p(
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

lib.getCutoffAbsorptionLines.restype = c.c_void_p
lib.getCutoffAbsorptionLines.argtypes = [c.c_void_p]

lib.getMirroringAbsorptionLines.restype = c.c_void_p
lib.getMirroringAbsorptionLines.argtypes = [c.c_void_p]

lib.getLineShapeTypeAbsorptionLines.restype = c.c_void_p
lib.getLineShapeTypeAbsorptionLines.argtypes = [c.c_void_p]

lib.getPopulationAbsorptionLines.restype = c.c_void_p
lib.getPopulationAbsorptionLines.argtypes = [c.c_void_p]

lib.getNormalizationAbsorptionLines.restype = c.c_void_p
lib.getNormalizationAbsorptionLines.argtypes = [c.c_void_p]

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

lib.getLocalQuantaAbsorptionLines.restype = c.c_void_p
lib.getLocalQuantaAbsorptionLines.argtypes = [c.c_void_p]

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
