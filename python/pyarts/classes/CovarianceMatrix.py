import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Block import Block
from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class CovarianceMatrix:
    """ ARTS CovarianceMatrix data

    Properties:
        sizeblocks:
            Size of blocks (Index)

        blocks:
            List of blocks of the covariance matrix (list of Block)

        sizeinverse_blocks:
            Size of inverse_blocks (Index)

        inverse_blocks:
            List of inverse blocks of the covariance matrix (list of Block)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createCovarianceMatrix())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @staticmethod
    def name():
        return "CovarianceMatrix"

    @property
    def sizeblocks(self):
        """ Size of blocks (Index) """
        return lib.sizeget_blocksCovarianceMatrix(self.__data__)

    @sizeblocks.setter
    def sizeblocks(self, x):
        assert x >= 0
        lib.resizeget_blocksCovarianceMatrix(int(x), self.__data__)

    @property
    def blocks(self):
        """ List of blocks of the covariance matrix (list of Block) """
        x = []
        n = self.sizeblocks
        for i in range(n):
            x.append(Block(c.c_void_p(lib.getelemget_blocksCovarianceMatrix(i, self.__data__))))
        return x

    @blocks.setter
    def blocks(self, val):
        if isinstance(val, Sized):
            self.sizeblocks = len(val)
            n = self.sizeblocks
            x = self.blocks
            for i in range(n):
                x[i].set(val[i])
        else:
            raise TypeError("Only accepts array-like input")

    @property
    def sizeinverse_blocks(self):
        """ Size of inverse_blocks (Index) """
        return lib.sizeget_inverse_blocksCovarianceMatrix(self.__data__)

    @sizeinverse_blocks.setter
    def sizeinverse_blocks(self, x):
        assert x >= 0
        lib.resizeget_inverse_blocksCovarianceMatrix(int(x), self.__data__)

    @property
    def inverse_blocks(self):
        """ List of inverse blocks of the covariance matrix (list of Block) """
        x = []
        n = self.sizeinverse_blocks
        for i in range(n):
            x.append(Block(c.c_void_p(lib.getelemget_inverse_blocksCovarianceMatrix(i, self.__data__))))
        return x

    @inverse_blocks.setter
    def inverse_blocks(self, val):
        if isinstance(val, Sized):
            self.sizeinverse_blocks = len(val)
            n = self.sizeinverse_blocks
            x = self.inverse_blocks
            for i in range(n):
                x[i].set(val[i])
        else:
            raise TypeError("Only accepts array-like input")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printCovarianceMatrix(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteCovarianceMatrix(self.__data__)

    def __repr__(self):
        return "ARTS CovarianceMatrix"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, CovarianceMatrix):
              self.blocks = other.blocks
              self.inverse_blocks = other.inverse_blocks
        else:
            raise TypeError("Expects CovarianceMatrix")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadCovarianceMatrix(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveCovarianceMatrix(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))
    
    def netcdfify(self):
        """ Create the NETCDF4 information required for writing this data
        
            Output: list that can be processed by netcdf.py, True arraytype
        """
        return [x.netcdfify() for x in self.blocks], True

    def denetcdf(self, group):
        """ Sets this based on a netcdf group
        
        Input:
            Group of data that can be interpreted as this's information
        """
        self.sizeinverse_blocks = 0
        self.sizeblocks = group.nelem
        for i in range(group.nelem):
            self.blocks[i].denetcdf(group["pos{}".format(i)])


lib.createCovarianceMatrix.restype = c.c_void_p
lib.createCovarianceMatrix.argtypes = []

lib.deleteCovarianceMatrix.restype = None
lib.deleteCovarianceMatrix.argtypes = [c.c_void_p]

lib.printCovarianceMatrix.restype = None
lib.printCovarianceMatrix.argtypes = [c.c_void_p]

lib.xmlreadCovarianceMatrix.restype = c.c_long
lib.xmlreadCovarianceMatrix.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveCovarianceMatrix.restype = c.c_long
lib.xmlsaveCovarianceMatrix.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.sizeget_blocksCovarianceMatrix.restype = c.c_long
lib.sizeget_blocksCovarianceMatrix.argtypes = [c.c_void_p]

lib.resizeget_blocksCovarianceMatrix.restype = None
lib.resizeget_blocksCovarianceMatrix.argtypes = [c.c_long, c.c_void_p]

lib.getelemget_blocksCovarianceMatrix.restype = c.c_void_p
lib.getelemget_blocksCovarianceMatrix.argtypes = [c.c_long, c.c_void_p]

lib.sizeget_inverse_blocksCovarianceMatrix.restype = c.c_long
lib.sizeget_inverse_blocksCovarianceMatrix.argtypes = [c.c_void_p]

lib.resizeget_inverse_blocksCovarianceMatrix.restype = None
lib.resizeget_inverse_blocksCovarianceMatrix.argtypes = [c.c_long, c.c_void_p]

lib.getelemget_inverse_blocksCovarianceMatrix.restype = c.c_void_p
lib.getelemget_inverse_blocksCovarianceMatrix.argtypes = [c.c_long, c.c_void_p]
