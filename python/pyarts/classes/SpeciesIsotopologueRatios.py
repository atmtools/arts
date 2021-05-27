import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.SpeciesIsotopeRecord import SpeciesIsotopeRecord

import numpy as np

class SpeciesIsotopologueRatios:
    """ ARTS SpeciesIsotopologueRatios data
    
    Only accessible via the access operator for indices valid according to
    SpeciesIsotopeRecord.max_len()
    """
    def __init__(self, data=None, delete=False):
        if isinstance(data, c.c_void_p):
            self.__delete__ = delete
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createSpeciesIsotopologueRatios())
    
    def __getitem__(self, ind):
        return self.data[ind]
    
    def __setitem__(self, ind, val):
        self.data[ind] = val
        
    @property
    def data(self):
        return np.ctypeslib.as_array(lib.getdataSpeciesIsotopologueRatios(self.__data__), [SpeciesIsotopeRecord.max_len()])

    def __del__(self):
        if self.__delete__:
            lib.deleteSpeciesIsotopologueRatios(self.__data__)

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printSpeciesIsotopologueRatios(self.__data__)
    
    def __repr__(self):
        return f"{self.data}"

lib.createSpeciesIsotopologueRatios.restype = c.c_void_p
lib.createSpeciesIsotopologueRatios.argtypes = []

lib.deleteSpeciesIsotopologueRatios.restype = None
lib.deleteSpeciesIsotopologueRatios.argtypes = [c.c_void_p]

lib.printSpeciesIsotopologueRatios.restype = None
lib.printSpeciesIsotopologueRatios.argtypes = [c.c_void_p]

lib.getdataSpeciesIsotopologueRatios.restype = c.POINTER(c.c_double)
lib.getdataSpeciesIsotopologueRatios.argtypes = [c.c_void_p]
