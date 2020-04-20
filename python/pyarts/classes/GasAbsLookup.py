import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.BasicTypes import ArrayOfIndex
from pyarts.classes.GridPosPoly import ArrayOfGridPosPoly
from pyarts.classes.Matrix import Matrix
from pyarts.classes.SpeciesTag import ArrayOfArrayOfSpeciesTag
from pyarts.classes.Tensor4 import Tensor4
from pyarts.classes.Vector import Vector
from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class GasAbsLookup:
    """ ARTS GasAbsLookup data

    Properties:
        specs:
            The species tags for which the table is valid (ArrayOfArrayOfSpeciesTag)

        nonlinspecs:
            The species tags with non-linear treatment (ArrayOfIndex)

        f_grid:
            The frequency grid [Hz] (Vector)

        fgp_default:
            Frequency grid positions (ArrayOfGridPosPoly)

        p_grid:
            The pressure grid for the table [Pa] (Vector)

        log_p_grid:
            The natural log of the pressure grid (Vector)

        vmrs:
            The reference VMR profiles (Matrix)

        t_ref:
            The reference temperature profile [K] (Vector)

        t_pert:
            The vector of temperature perturbations [K] (Vector)

        nls_pert:
            The vector of perturbations for the VMRs of the nonlinear species (Vector)

        xsec:
            Absorption cross sections (Tensor4)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createGasAbsLookup())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def specs(self):
        """ The species tags for which the table is valid (ArrayOfArrayOfSpeciesTag) """
        return ArrayOfArrayOfSpeciesTag(c.c_void_p(lib.getSpeciesGasAbsLookup(self.__data__)))

    @specs.setter
    def specs(self, val):
        self.specs.set(val)

    @property
    def nonlinspecs(self):
        """ The species tags with non-linear treatment (ArrayOfIndex) """
        return ArrayOfIndex(c.c_void_p(lib.getNonLinearSpeciesGasAbsLookup(self.__data__)))

    @nonlinspecs.setter
    def nonlinspecs(self, val):
        self.nonlinspecs.set(val)

    @property
    def f_grid(self):
        """ The frequency grid [Hz] (Vector) """
        return Vector(c.c_void_p(lib.getFgridGasAbsLookup(self.__data__)))

    @f_grid.setter
    def f_grid(self, val):
        self.f_grid.set(val)

    @property
    def fgp_default(self):
        """ Frequency grid positions (ArrayOfGridPosPoly) """
        return ArrayOfGridPosPoly(c.c_void_p(lib.getFGPDefaultGasAbsLookup(self.__data__)))

    @fgp_default.setter
    def fgp_default(self, val):
        self.fgp_default.set(val)

    @property
    def p_grid(self):
        """ The pressure grid for the table [Pa] (Vector) """
        return Vector(c.c_void_p(lib.getPgridGasAbsLookup(self.__data__)))

    @p_grid.setter
    def p_grid(self, val):
        self.p_grid.set(val)

    @property
    def log_p_grid(self):
        """ The natural log of the pressure grid (Vector) """
        return Vector(c.c_void_p(lib.getLogPgridGasAbsLookup(self.__data__)))

    @log_p_grid.setter
    def log_p_grid(self, val):
        self.log_p_grid.set(val)

    @property
    def vmrs(self):
        """ The reference VMR profiles (Matrix) """
        return Matrix(c.c_void_p(lib.getVMRsGasAbsLookup(self.__data__)))

    @vmrs.setter
    def vmrs(self, val):
        self.vmrs.set(val)

    @property
    def t_ref(self):
        """ The reference temperature profile [K] (Vector) """
        return Vector(c.c_void_p(lib.getTrefGasAbsLookup(self.__data__)))

    @t_ref.setter
    def t_ref(self, val):
        self.t_ref.set(val)

    @property
    def t_pert(self):
        """ The vector of temperature perturbations [K] (Vector) """
        return Vector(c.c_void_p(lib.getTpertGasAbsLookup(self.__data__)))

    @t_pert.setter
    def t_pert(self, val):
        self.t_pert.set(val)

    @property
    def nls_pert(self):
        """ The vector of perturbations for the VMRs of the nonlinear species (Vector) """
        return Vector(c.c_void_p(lib.getNLSPertGasAbsLookup(self.__data__)))

    @nls_pert.setter
    def nls_pert(self, val):
        self.nls_pert.set(val)

    @property
    def xsec(self):
        """ Absorption cross sections (Tensor4) """
        return Tensor4(c.c_void_p(lib.getXsecGasAbsLookup(self.__data__)))

    @xsec.setter
    def xsec(self, val):
        self.xsec.set(val)

    @staticmethod
    def name():
        return "GasAbsLookup"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printGasAbsLookup(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteGasAbsLookup(self.__data__)

    def __repr__(self):
        return "ARTS GasAbsLookup"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, GasAbsLookup):
            self.specs = other.specs
            self.nonlinspecs = other.nonlinspecs
            self.f_grid = other.f_grid
            self.fgp_default = other.fgp_default
            self.p_grid = other.p_grid
            self.log_p_grid = other.log_p_grid
            self.vmrs = other.vmrs
            self.t_ref = other.t_ref
            self.t_pert = other.t_pert
            self.nls_pert = other.nls_pert
            self.xsec = other.xsec
        else:
            raise TypeError("Expects GasAbsLookup")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadGasAbsLookup(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveGasAbsLookup(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, GasAbsLookup) and \
                self.specs == other.specs and \
                self.nonlinspecs == other.nonlinspecs and \
                self.f_grid == other.f_grid and \
                self.fgp_default == other.fgp_default and \
                self.p_grid == other.p_grid and \
                self.log_p_grid == other.log_p_grid and \
                self.vmrs == other.vmrs and \
                self.t_ref == other.t_ref and \
                self.t_pert == other.t_pert and \
                self.nls_pert == other.nls_pert and \
                self.xsec == other.xsec:
            return True
        else:
            return False

    def __bool__(self):
        return bool(self.xsec)


lib.createGasAbsLookup.restype = c.c_void_p
lib.createGasAbsLookup.argtypes = []

lib.deleteGasAbsLookup.restype = None
lib.deleteGasAbsLookup.argtypes = [c.c_void_p]

lib.printGasAbsLookup.restype = None
lib.printGasAbsLookup.argtypes = [c.c_void_p]

lib.xmlreadGasAbsLookup.restype = c.c_long
lib.xmlreadGasAbsLookup.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveGasAbsLookup.restype = c.c_long
lib.xmlsaveGasAbsLookup.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getSpeciesGasAbsLookup.restype = c.c_void_p
lib.getSpeciesGasAbsLookup.argtypes = [c.c_void_p]

lib.getNonLinearSpeciesGasAbsLookup.restype = c.c_void_p
lib.getNonLinearSpeciesGasAbsLookup.argtypes = [c.c_void_p]

lib.getFgridGasAbsLookup.restype = c.c_void_p
lib.getFgridGasAbsLookup.argtypes = [c.c_void_p]

lib.getFGPDefaultGasAbsLookup.restype = c.c_void_p
lib.getFGPDefaultGasAbsLookup.argtypes = [c.c_void_p]

lib.getPgridGasAbsLookup.restype = c.c_void_p
lib.getPgridGasAbsLookup.argtypes = [c.c_void_p]

lib.getLogPgridGasAbsLookup.restype = c.c_void_p
lib.getLogPgridGasAbsLookup.argtypes = [c.c_void_p]

lib.getVMRsGasAbsLookup.restype = c.c_void_p
lib.getVMRsGasAbsLookup.argtypes = [c.c_void_p]

lib.getTrefGasAbsLookup.restype = c.c_void_p
lib.getTrefGasAbsLookup.argtypes = [c.c_void_p]

lib.getTpertGasAbsLookup.restype = c.c_void_p
lib.getTpertGasAbsLookup.argtypes = [c.c_void_p]

lib.getNLSPertGasAbsLookup.restype = c.c_void_p
lib.getNLSPertGasAbsLookup.argtypes = [c.c_void_p]

lib.getXsecGasAbsLookup.restype = c.c_void_p
lib.getXsecGasAbsLookup.argtypes = [c.c_void_p]
