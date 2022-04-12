import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base

import numpy as np
import scipy as sp

class Sparse:
    """ ARTS Sparse data

    Properties:
            Dimensionality of Sparse (constexpr Index)

        data:
            The data (scipy-sparse-csr_matrix)

        shape:
            Shape of the Sparse (tuple)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSparse())
            if data is not None:
                self.data = data

    @property
    def dims(self):
        return 2

    @property
    def flat(self):
        """ flat iterator of raw data """
        return self.raw.flat

    @property
    def raw(self):
        """ Raw data (numpy-array) """
        return np.ctypeslib.as_array(lib.getDataSparse(self.__data__), (self.size,))

    @property
    def col_ind(self):
        """ Column index (numpy-array) """
        return np.ctypeslib.as_array(lib.colsptrSparse(self.__data__), (self.size,))

    @property
    def row_ind_ptr(self):
        """ Row pointers (numpy-array) """
        x = np.ctypeslib.as_array(lib.rowsptrSparse(self.__data__), (self.shape[0]+1,))
        x[-1] = self.size
        return x

    @property
    def data(self):
        """ The data (scipy-sparse-csr_matrix) """
        if self.size == 0:
            return sp.sparse.csr_matrix(self.shape, dtype=float)
        else:
            s = self.shape
            r = self.col_ind
            c = self.row_ind_ptr
            return sp.sparse.csr_matrix((self.raw, r, c), shape=s)

    @data.setter
    def data(self, val):
        x = sp.sparse.csr_matrix(val, copy=True)
        self.shape = x.shape
        rows, cols, vals = sp.sparse.find(x)
        for i in range(len(rows)):
            lib.setDataSparse(self.__data__, int(rows[i]), int(cols[i]), float(vals[i]))

    @property
    def size(self):
        return lib.sizeSparse(self.__data__)

    @property
    def shape(self):
        return lib.rowsSparse(self.__data__), lib.colsSparse(self.__data__),

    @shape.setter
    def shape(self, shape):
        if len(shape) == self.dims:
            lib.resizeSparse(*shape, self.__data__)
        else:
            raise RuntimeError("Cannot resize dimension")

    @staticmethod
    def name():
        return "Sparse"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSparse(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteSparse(self.__data__)

    def __repr__(self):
        return str(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Sparse):
            self.data = other.data
        else:
            raise TypeError("Expects Sparse")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadSparse(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveSparse(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))
    
    def netcdfify(self):
        """ Create the NETCDF4 information required for writing this data
        
        Output: list that can be processed by netcdf.py, False arraytype
        """
        return [["data", self.data.data, float, {"size": self.size}],
                ["row_ind_ptr", self.row_ind_ptr, np.int32, {"nrowsp1": len(self.row_ind_ptr)}],
                ["col_ind", self.col_ind, np.int32, {"size": self.size}],
                ["shape", self.shape, int, {"two": 2}]], False
    
    def denetcdf(self, group):
        """ Sets this based on a netcdf group
        
        Input:
            Group of data that can be interpreted as this's information
        """
        self.data = \
            sp.sparse.csr_matrix((np.array(group["data"]),
                                  np.array(group["col_ind"]),
                                  np.array(group["row_ind_ptr"])),
                                 shape=np.array(group["shape"]))

    def __array__(self):
        return np.array(self.data.todense())

    def __eq__(self, other):
        if isinstance(other, Sparse):
            return bool((self.data == other.data).data.all())
        else:
            return self.data == other

    def __bool__(self):
        return bool(self.size)

    def __iadd__(self, val):
        self.data += np.array(val)
        return self

    def __isub__(self, val):
        self.data -= np.array(val)
        return self

    def __imul__(self, val):
        self.data *= np.array(val)
        return self

    def __imatmul__(self, val):
        self.data = self.data @ np.array(val)
        return self

    def __itruediv__(self, val):
        self.data /= np.array(val)
        return self

    def __ipow__(self, val):
        self.data **= np.array(val)
        return self

    def __add__(self, val):
        return Sparse(self.data + np.array(val))

    def __sub__(self, val):
        return Sparse(self.data - np.array(val))

    def __mul__(self, val):
        return Sparse(self.data * np.array(val))

    def __matmul__(self, val):
        return self.data @ np.array(val)

    def __truediv__(self, val):
        return Sparse(self.data / np.array(val))

    def __pow__(self, val):
        return Sparse(self.data ** np.array(val))

    def __radd__(self, val):
        return Sparse(np.array(val) + self.data)

    def __rsub__(self, val):
        return Sparse(np.array(val) - self.data)

    def __rmul__(self, val):
        return Sparse(np.array(val) * self.data)

    def __rmatmul__(self, val):
        return np.array(val) @ self.data

    def __rtruediv__(self, val):
        return Sparse(np.array(val) / self.data)

    def __rpow__(self, val):
        return Sparse(np.array(val) ** self.data)


exec(array_base(Sparse))


lib.createSparse.restype = c.c_void_p
lib.createSparse.argtypes = []

lib.deleteSparse.restype = None
lib.deleteSparse.argtypes = [c.c_void_p]

lib.printSparse.restype = None
lib.printSparse.argtypes = [c.c_void_p]

lib.xmlreadSparse.restype = c.c_long
lib.xmlreadSparse.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveSparse.restype = c.c_long
lib.xmlsaveSparse.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.colsSparse.restype = c.c_long
lib.colsSparse.argtypes = [c.c_void_p]

lib.rowsSparse.restype = c.c_long
lib.rowsSparse.argtypes = [c.c_void_p]

lib.resizeMatrix.restype = None
lib.resizeMatrix.argtypes = [c.c_long, c.c_long, c.c_void_p]

lib.colsptrSparse.restype = c.POINTER(c.c_int)
lib.colsptrSparse.argtypes = [c.c_void_p]

lib.rowsptrSparse.restype = c.POINTER(c.c_int)
lib.rowsptrSparse.argtypes = [c.c_void_p]

lib.sizeSparse.restype = c.c_long
lib.sizeSparse.argtypes = [c.c_void_p]

lib.getDataSparse.restype = c.POINTER(c.c_double)
lib.getDataSparse.argtypes = [c.c_void_p]

lib.setDataSparse.restype = None
lib.setDataSparse.argtypes = [c.c_void_p, c.c_long, c.c_long, c.c_double]
