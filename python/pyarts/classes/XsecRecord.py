import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.GriddedField2 import ArrayOfGriddedField2
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class XsecRecord:
    """ ARTS XsecRecord data

    Properties:
        version:
            Version (Index)

        spec:
            Species (Species)

        coeffs:
            Coefficients (Vector)

        ref_pressure:
            Reference pressure (Vector)

        ref_temperature:
            Reference temperature (Vector)

        fgrids:
            Frequency grids (ArrayOfVector)

        xsecs:
            Cross-sections (ArrayOfVector)

        temperature_slope:
            Slope of temperature fits (ArrayOfVector)

        temperature_intersect:
            Intersect of temperature fits (ArrayOfVector)

        fitminpressures:
            Mininum pressures of fit

        fitmaxpressures:
            Maximum pressures of fit

        fitmintemperatures:
            Mininum temperatures of fit

        fitmaxtemperatures:
            Mininum temperatures of fit

        fitcoeffs:
            Fit coefficients
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createXsecRecord())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def version(self):
        """ Version (Index) """
        return lib.getVersionXsecRecord(self.__data__)

    @version.setter
    def version(self, val):
        lib.setVersionXsecRecord(self.__data__, int(val))

    @property
    def spec(self):
        """ Species (Species) """
        return Species(c.c_void_p(lib.getSpeciesXsecRecord(self.__data__)))

    @spec.setter
    def spec(self, val):
        spec = val if isinstance(val, Species) else Species(val)
        lib.setSpeciesXsecRecord(self.__data__, spec.__data__)

    @property
    def coeffs(self):
        """ Coefficients (Vector) """
        return Vector(c.c_void_p(lib.getCoeffsXsecRecord(self.__data__)))

    @coeffs.setter
    def coeffs(self, val):
        self.coeffs.set(val)

    @property
    def ref_pressure(self):
        """ Reference pressure (Vector) """
        return Vector(c.c_void_p(lib.getRefPressureXsecRecord(self.__data__)))

    @ref_pressure.setter
    def ref_pressure(self, val):
        self.ref_pressure.set(val)

    @property
    def ref_temperature(self):
        """ Reference temperature (Vector) """
        return Vector(c.c_void_p(lib.getRefTemperatureXsecRecord(self.__data__)))

    @ref_temperature.setter
    def ref_temperature(self, val):
        self.ref_temperature.set(val)

    @property
    def fgrids(self):
        """ Frequency grids (ArrayOfVector) """
        return ArrayOfVector(c.c_void_p(lib.getFgridsXsecRecord(self.__data__)))

    @fgrids.setter
    def fgrids(self, val):
        self.fgrids.set(val)

    @property
    def xsecs(self):
        """ Cross-sections (ArrayOfVector) """
        return ArrayOfVector(c.c_void_p(lib.getXsecsXsecRecord(self.__data__)))

    @xsecs.setter
    def xsecs(self, val):
        self.xsecs.set(val)

    @property
    def temperature_slope(self):
        """ Slope of temperature fits (ArrayOfVector) """
        return ArrayOfVector(c.c_void_p(lib.getTemperatureSlopeXsecRecord(self.__data__)))

    @temperature_slope.setter
    def temperature_slope(self, val):
        self.temperature_slope.set(val)

    @property
    def temperature_intersect(self):
        """ Intersect of temperature fits (ArrayOfVector) """
        return ArrayOfVector(c.c_void_p(lib.getTemperatureIntersectXsecRecord(self.__data__)))

    @temperature_intersect.setter
    def temperature_intersect(self, val):
        self.temperature_intersect.set(val)

    @property
    def fitminpressures(self):
        """ Minimum pressures of fit (Vector) """
        return Vector(c.c_void_p(lib.getFitMinPressuresXsecRecord(self.__data__)))

    @fitminpressures.setter
    def fitminpressures(self, val):
        self.fitminpressures.set(val)

    @property
    def fitmaxpressures(self):
        """ Maximum pressures of fit (Vector) """
        return Vector(c.c_void_p(lib.getFitMaxPressuresXsecRecord(self.__data__)))

    @fitmaxpressures.setter
    def fitmaxpressures(self, val):
        self.fitmaxpressures.set(val)

    @property
    def fitmintemperatures(self):
        """ Minimum temperatures of fit (Vector) """
        return Vector(c.c_void_p(lib.getFitMinTemperaturesXsecRecord(self.__data__)))

    @fitmintemperatures.setter
    def fitmintemperatures(self, val):
        self.fitmintemperatures.set(val)

    @property
    def fitmaxtemperatures(self):
        """ Maximum temperatures of fit (Vector) """
        return Vector(c.c_void_p(lib.getFitMaxTemperaturesXsecRecord(self.__data__)))

    @fitmaxtemperatures.setter
    def fitmaxtemperatures(self, val):
        self.fitmaxtemperatures.set(val)

    @property
    def fitcoeffs(self):
        """ Fit coefficients (ArrayOfGriddedField2) """
        return ArrayOfGriddedField2(c.c_void_p(lib.getFitCoeffsXsecRecord(self.__data__)))

    @fitcoeffs.setter
    def fitcoeffs(self, val):
        self.fitcoeffs.set(val)

    @staticmethod
    def name():
        return "XsecRecord"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printXsecRecord(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteXsecRecord(self.__data__)

    def __repr__(self):
        return "ARTS XsecRecord"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, XsecRecord):
            lib.setVersionXsecRecord(self.__data__, int(other.version))
            lib.setSpeciesXsecRecord(self.__data__, int(other.spec))
            self.coeffs = other.coeffs
            self.ref_pressure = other.ref_pressure
            self.ref_temperature = other.ref_temperature
            self.fgrids = other.fgrids
            self.xsecs = other.xsecs
            self.temperature_slope = other.temperature_slope
            self.temperature_intersect = other.temperature_intersect
            self.fitminpressures = other.fitminpressures
            self.fitmaxpressures = other.fitmaxpressures
            self.fitmintemperatures = other.fitmintemperatures
            self.fitmaxtemperatures = other.fitmaxtemperatures
            self.fitcoeffs = other.fitcoeffs
        else:
            raise TypeError("Expects XsecRecord")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadXsecRecord(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveXsecRecord(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, XsecRecord) and \
                self.version == other.version and \
                self.spec == other.spec and \
                self.coeffs == other.coeffs and \
                self.ref_pressure == other.ref_pressure and \
                self.ref_temperature == other.ref_temperature and \
                self.fgrids == other.fgrids and \
                self.xsecs == other.xsecs and \
                self.temperature_slope == other.temperature_slope and \
                self.temperature_intersect == other.temperature_intersect and \
                self.fitminpressures == other.fitminpressures and \
                self.fitmaxpressures == other.fitmaxpressures and \
                self.fitmintemperatures == other.fitmintemperatures and \
                self.fitmaxtemperatures == other.fitmaxtemperatures and \
                self.fitcoeffs == other.fitcoeffs:
            return True
        else:
            return False


exec(array_base(XsecRecord))


lib.createXsecRecord.restype = c.c_void_p
lib.createXsecRecord.argtypes = []

lib.deleteXsecRecord.restype = None
lib.deleteXsecRecord.argtypes = [c.c_void_p]

lib.printXsecRecord.restype = None
lib.printXsecRecord.argtypes = [c.c_void_p]

lib.xmlreadXsecRecord.restype = c.c_long
lib.xmlreadXsecRecord.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveXsecRecord.restype = c.c_long
lib.xmlsaveXsecRecord.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getVersionXsecRecord.restype = c.c_long
lib.getVersionXsecRecord.argtypes = [c.c_void_p]

lib.getSpeciesXsecRecord.restype = c.c_long
lib.getSpeciesXsecRecord.argtypes = [c.c_void_p]

lib.setSpeciesXsecRecord.restype = None
lib.setSpeciesXsecRecord.argtypes = [c.c_void_p, c.c_void_p]

lib.getCoeffsXsecRecord.restype = c.c_void_p
lib.getCoeffsXsecRecord.argtypes = [c.c_void_p]

lib.getRefPressureXsecRecord.restype = c.c_void_p
lib.getRefPressureXsecRecord.argtypes = [c.c_void_p]

lib.getRefTemperatureXsecRecord.restype = c.c_void_p
lib.getRefTemperatureXsecRecord.argtypes = [c.c_void_p]

lib.getFgridsXsecRecord.restype = c.c_void_p
lib.getFgridsXsecRecord.argtypes = [c.c_void_p]

lib.getXsecsXsecRecord.restype = c.c_void_p
lib.getXsecsXsecRecord.argtypes = [c.c_void_p]

lib.getTemperatureSlopeXsecRecord.restype = c.c_void_p
lib.getTemperatureSlopeXsecRecord.argtypes = [c.c_void_p]

lib.getTemperatureIntersectXsecRecord.restype = c.c_void_p
lib.getTemperatureIntersectXsecRecord.argtypes = [c.c_void_p]

lib.getFitMinPressuresXsecRecord.restype = c.c_void_p
lib.getFitMinPressuresXsecRecord.argtypes = [c.c_void_p]

lib.getFitMaxPressuresXsecRecord.restype = c.c_void_p
lib.getFitMaxPressuresXsecRecord.argtypes = [c.c_void_p]

lib.getFitMinTemperaturesXsecRecord.restype = c.c_void_p
lib.getFitMinTemperaturesXsecRecord.argtypes = [c.c_void_p]

lib.getFitMaxTemperaturesXsecRecord.restype = c.c_void_p
lib.getFitMaxTemperaturesXsecRecord.argtypes = [c.c_void_p]

lib.getFitCoeffsXsecRecord.restype = c.c_void_p
lib.getFitCoeffsXsecRecord.argtypes = [c.c_void_p]
