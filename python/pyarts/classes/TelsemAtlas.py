import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import String, Numeric, Index, ArrayOfIndex
from pyarts.classes.Matrix import Matrix
from pyarts.classes.Tensor3 import Tensor3
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class TelsemAtlas:
    """ ARTS TelsemAtlas data

    Properties:
        ndat:
            Number of lines in the Atlas (Index)

        nchan:
            Number of channels in the Atlas (Index)

        atlasname:
            Name of the atlas (including version number) (String)

        month:
            Month of the Atlas (Index)

        dlat:
            Resolution of the Atlas (Numeric)

        ncells:
            Number of cells per lat band (ArrayOfIndex)

        firstcells:
            The first cell number of lat band (ArrayOfIndex)

        emis:
            Emissivities (Matrix)

        emis_err:
            Emissivity uncertainties (Matrix)

        correl:
            Emissivity correlations (Tensor3)

        classes1:
            Surface class (ArrayOfIndex)

        classes2:
            Surface class (ArrayOfIndex)

        cellnums:
            Cellnumber of each of the pixels in the atlas (ArrayOfIndex)

        correspondance:
            Derived from file data (ArrayOfIndex)

        a0_k0:
           Regression coefficients (constexpr list of 30 Numeric)

        a0_k1:
           Regression coefficients (constexpr list of 30 Numeric)

        a0_k2:
           Regression coefficients (constexpr list of 30 Numeric)

        a0_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        a1_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        a2_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        a3_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        b0_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        b1_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        b2_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        b3_eveh:
           Regression coefficients (constexpr list of 30 Numeric)

        rapport43_32:
           Regression coefficients (constexpr list of 4 Numeric)

        rapport54_43:
           Regression coefficients (constexpr list of 4 Numeric)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTelsemAtlas())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def ndat(self):
        """ Number of lines in the Atlas (Index) """
        return Index(c.c_void_p(lib.getDataCountTelsemAtlas(self.__data__)))

    @property
    def nchan(self):
        """ Number of channels in the Atlas (Index) """
        return Index(c.c_void_p(lib.getChannelCountTelsemAtlas(self.__data__)))

    @property
    def atlasname(self):
        """ Name of the atlas (including version number) (String) """
        return String(c.c_void_p(lib.getNameTelsemAtlas(self.__data__)))

    @property
    def month(self):
        """ Month of the Atlas (Index) """
        return Index(c.c_void_p(lib.getMonthTelsemAtlas(self.__data__)))

    @property
    def dlat(self):
        """ Resolution of the Atlas (Numeric) """
        return Numeric(c.c_void_p(lib.getLatTelsemAtlas(self.__data__)))

    @property
    def ncells(self):
        """ Number of cells per lat band (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getCellsTelsemAtlas(self.__data__)))

    @property
    def firstcells(self):
        """ The first cell number of lat band (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getFirstCellsTelsemAtlas(self.__data__)))

    @property
    def emis(self):
        """ Emissivities (Matrix) """
        return Matrix(c.c_void_p(lib.getEmisTelsemAtlas(self.__data__)))

    @property
    def emis_err(self):
        """ Emissivity uncertainties (Matrix) """
        return Matrix(c.c_void_p(lib.getEmis_errTelsemAtlas(self.__data__)))

    @property
    def correl(self):
        """ Emissivity correlations (Tensor3) """
        return Tensor3(c.c_void_p(lib.getCorrelationsTelsemAtlas(self.__data__)))

    @property
    def classes1(self):
        """ Surface class (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getClasses1TelsemAtlas(self.__data__)))

    @property
    def classes2(self):
        """ Surface class (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getClasses2TelsemAtlas(self.__data__)))

    @property
    def cellnums(self):
        """ Cellnumber of each of the pixels in the atlas (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getCellnumberTelsemAtlas(self.__data__)))

    @property
    def correspondance(self):
        """ Derived from file data (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getCorrespondanceTelsemAtlas(self.__data__)))

    @property
    def a0_k0(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getA0_K0TelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def a0_k1(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getA0_K1TelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def a0_k2(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getA0_K2TelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def a0_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getA0_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def a1_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getA1_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def a2_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getA2_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def a3_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getA3_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def b0_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getB0_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def b1_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getB1_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def b2_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getB2_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def b3_eveh(self):
       """ Regression coefficients (constexpr list of 30 Numeric) """
       return [lib.getB3_EVEHTelsemAtlas(int(i), self.__data__) for i in range(30)]

    @property
    def rapport43_32(self):
       """ Regression coefficients (constexpr list of 4 Numeric) """
       return [lib.getRAPPORT43_32TelsemAtlas(int(i), self.__data__) for i in range(4)]

    @property
    def rapport54_43(self):
       """ Regression coefficients (constexpr list of 4 Numeric) """
       return [lib.getRAPPORT54_43TelsemAtlas(int(i), self.__data__) for i in range(4)]

    @ndat.setter
    def ndat(self, x):
        self.ndat.set(x)

    @nchan.setter
    def nchan(self, x):
        self.nchan.set(x)

    @atlasname.setter
    def atlasname(self, x):
        self.atlasname.set(x)

    @month.setter
    def month(self, x):
        self.month.set(x)

    @dlat.setter
    def dlat(self, x):
        self.dlat.set(x)

    @ncells.setter
    def ncells(self, x):
        self.ncells.set(x)

    @firstcells.setter
    def firstcells(self, x):
        self.firstcells.set(x)

    @emis.setter
    def emis(self, x):
        self.emis.set(x)

    @emis_err.setter
    def emis_err(self, x):
        self.emis_err.set(x)

    @correl.setter
    def correl(self, x):
        self.correl.set(x)

    @classes1.setter
    def classes1(self, x):
        self.classes1.set(x)

    @classes2.setter
    def classes2(self, x):
        self.classes2.set(x)

    @cellnums.setter
    def cellnums(self, x):
        self.cellnums.set(x)

    @correspondance.setter
    def correspondance(self, x):
        self.correspondance.set(x)

    @staticmethod
    def name():
        return "TelsemAtlas"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTelsemAtlas(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTelsemAtlas(self.__data__)

    def __repr__(self):
        return "ARTS TelsemAtlas"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, TelsemAtlas):
            self.ndat = other.ndat
            self.nchan = other.nchan
            self.atlasname = other.atlasname
            self.month = other.month
            self.dlat = other.dlat
            self.ncells = other.ncells
            self.firstcells = other.firstcells
            self.emis = other.emis
            self.emis_err = other.emis_err
            self.correl = other.correl
            self.classes1 = other.classes1
            self.classes2 = other.classes2
            self.cellnums = other.cellnums
            self.correspondance = other.correspondance
        else:
            raise TypeError("Expects TelsemAtlas")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTelsemAtlas(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTelsemAtlas(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, TelsemAtlas) and \
                self.ndat == other.ndat and \
                self.nchan == other.nchan and \
                self.atlasname == other.atlasname and \
                self.month == other.month and \
                self.dlat == other.dlat and \
                self.ncells == other.ncells and \
                self.firstcells == other.firstcells and \
                self.emis == other.emis and \
                self.emis_err == other.emis_err and \
                self.correl == other.correl and \
                self.classes1 == other.classes1 and \
                self.classes2 == other.classes2 and \
                self.cellnums == other.cellnums and \
                self.correspondance == other.correspondance:
            return True
        else:
            return False


exec(array_base(TelsemAtlas))


lib.createTelsemAtlas.restype = c.c_void_p
lib.createTelsemAtlas.argtypes = []

lib.deleteTelsemAtlas.restype = None
lib.deleteTelsemAtlas.argtypes = [c.c_void_p]

lib.printTelsemAtlas.restype = None
lib.printTelsemAtlas.argtypes = [c.c_void_p]

lib.xmlreadTelsemAtlas.restype = c.c_long
lib.xmlreadTelsemAtlas.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTelsemAtlas.restype = c.c_long
lib.xmlsaveTelsemAtlas.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getDataCountTelsemAtlas.restype = c.c_void_p
lib.getDataCountTelsemAtlas.argtypes = [c.c_void_p]

lib.getChannelCountTelsemAtlas.restype = c.c_void_p
lib.getChannelCountTelsemAtlas.argtypes = [c.c_void_p]

lib.getNameTelsemAtlas.restype = c.c_void_p
lib.getNameTelsemAtlas.argtypes = [c.c_void_p]

lib.getMonthTelsemAtlas.restype = c.c_void_p
lib.getMonthTelsemAtlas.argtypes = [c.c_void_p]

lib.getLatTelsemAtlas.restype = c.c_void_p
lib.getLatTelsemAtlas.argtypes = [c.c_void_p]

lib.getCellsTelsemAtlas.restype = c.c_void_p
lib.getCellsTelsemAtlas.argtypes = [c.c_void_p]

lib.getFirstCellsTelsemAtlas.restype = c.c_void_p
lib.getFirstCellsTelsemAtlas.argtypes = [c.c_void_p]

lib.getEmisTelsemAtlas.restype = c.c_void_p
lib.getEmisTelsemAtlas.argtypes = [c.c_void_p]

lib.getEmis_errTelsemAtlas.restype = c.c_void_p
lib.getEmis_errTelsemAtlas.argtypes = [c.c_void_p]

lib.getCorrelationsTelsemAtlas.restype = c.c_void_p
lib.getCorrelationsTelsemAtlas.argtypes = [c.c_void_p]

lib.getClasses1TelsemAtlas.restype = c.c_void_p
lib.getClasses1TelsemAtlas.argtypes = [c.c_void_p]

lib.getClasses2TelsemAtlas.restype = c.c_void_p
lib.getClasses2TelsemAtlas.argtypes = [c.c_void_p]

lib.getCellnumberTelsemAtlas.restype = c.c_void_p
lib.getCellnumberTelsemAtlas.argtypes = [c.c_void_p]

lib.getCorrespondanceTelsemAtlas.restype = c.c_void_p
lib.getCorrespondanceTelsemAtlas.argtypes = [c.c_void_p]

lib.getA0_K0TelsemAtlas.restype = c.c_double
lib.getA0_K0TelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getA0_K1TelsemAtlas.restype = c.c_double
lib.getA0_K1TelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getA0_K2TelsemAtlas.restype = c.c_double
lib.getA0_K2TelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getA0_EVEHTelsemAtlas.restype = c.c_double
lib.getA0_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getA1_EVEHTelsemAtlas.restype = c.c_double
lib.getA1_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getA2_EVEHTelsemAtlas.restype = c.c_double
lib.getA2_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getA3_EVEHTelsemAtlas.restype = c.c_double
lib.getA3_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getB0_EVEHTelsemAtlas.restype = c.c_double
lib.getB0_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getB1_EVEHTelsemAtlas.restype = c.c_double
lib.getB1_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getB2_EVEHTelsemAtlas.restype = c.c_double
lib.getB2_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getB3_EVEHTelsemAtlas.restype = c.c_double
lib.getB3_EVEHTelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getRAPPORT43_32TelsemAtlas.restype = c.c_double
lib.getRAPPORT43_32TelsemAtlas.argtypes = [c.c_long, c.c_void_p]

lib.getRAPPORT54_43TelsemAtlas.restype = c.c_double
lib.getRAPPORT54_43TelsemAtlas.argtypes = [c.c_long, c.c_void_p]
