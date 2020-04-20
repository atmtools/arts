import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Range import Range
from pyarts.classes.Matrix import Matrix
from pyarts.classes.Sparse import Sparse

class Block:
    """ ARTS Block data

    Properties:
        type:
            Type of data (const Index, Matrix: 0; Sparse: 1)

        row_range:
            Row range of block (Range)

        col_range:
            Column range of block (Range)

        indices:
            The indices of the retrieval quantities correlated (Index tuple)

        data:
            Data (Matrix/Sparse)
    """
    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = False
            self.__data__ = c.c_void_p(lib.createBlock())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def type(self):
        """ Type of data (const Index, Matrix: 0; Sparse: 1) """
        return lib.get_matrix_typeBlock(self.__data__)

    @property
    def row_range(self):
        """ Row range of block (Range) """
        return Range(c.c_void_p(lib.getget_row_rangeBlock(self.__data__)))

    @property
    def col_range(self):
        """ Column range of block (Range) """
        return Range(c.c_void_p(lib.getget_row_rangeBlock(self.__data__)))

    @row_range.setter
    def row_range(self, other):
        self.row_range.set(other)

    @col_range.setter
    def col_range(self, other):
        self.col_range.set(other)

    @property
    def indices(self):
        """ The indices of the retrieval quantities correlated (Index tuple) """
        return lib.get_index1Block(self.__data__), lib.get_index2Block(self.__data__)

    @indices.setter
    def indices(self, x):
        lib.set_indicesBlock(self.__data__, int(x[0]), int(x[1]))

    @property
    def data(self):
        if self.type == 0:
            x = lib.getget_denseBlock(self.__data__)
            if x:
                return Matrix(c.c_void_p(x))
            else:
                raise RuntimeError("No data")
        elif self.type == 1:
            x = lib.getget_sparseBlock(self.__data__)
            if x:
                return Sparse(c.c_void_p(x))
            else:
                raise RuntimeError("No data")
        else:
            raise TypeError("Cannot understand self's type")

    @data.setter
    def data(self, x):
        if isinstance(x, Matrix):
            lib.set_matrixBlock(self.__data__, x.__data__, True)
        elif isinstance(x, Sparse):
            lib.set_matrixBlock(self.__data__, x.__data__, False)
        else:
            raise TypeError("Expects Sparse or Matrix")

    @staticmethod
    def name():
        return "Block"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printBlock(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteBlock(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Block):
            self.row_range = other.row_range
            self.col_range = other.col_range
            self.indices = other.indices
            self.data = other.data
        else:
            raise TypeError("Expects Block")


lib.createBlock.restype = c.c_void_p
lib.createBlock.argtypes = []

lib.deleteBlock.restype = None
lib.deleteBlock.argtypes = [c.c_void_p]

lib.printBlock.restype = None
lib.printBlock.argtypes = [c.c_void_p]

lib.getget_row_rangeBlock.restype = c.c_void_p
lib.getget_row_rangeBlock.argtypes = [c.c_void_p]

lib.getget_column_rangeBlock.restype = c.c_void_p
lib.getget_column_rangeBlock.argtypes = [c.c_void_p]

lib.getget_denseBlock.restype = c.c_void_p
lib.getget_denseBlock.argtypes = [c.c_void_p]

lib.getget_sparseBlock.restype = c.c_void_p
lib.getget_sparseBlock.argtypes = [c.c_void_p]

lib.get_matrix_typeBlock.restype = c.c_long
lib.get_matrix_typeBlock.argtypes = [c.c_void_p]

lib.get_index1Block.restype = c.c_long
lib.get_index1Block.argtypes = [c.c_void_p]

lib.get_index2Block.restype = c.c_long
lib.get_index2Block.argtypes = [c.c_void_p]

lib.set_indicesBlock.restype = None
lib.set_indicesBlock.argtypes = [c.c_void_p, c.c_long, c.c_long]

lib.set_matrixBlock.restype = None
lib.set_matrixBlock.argtypes = [c.c_void_p, c.c_void_p, c.c_bool]
