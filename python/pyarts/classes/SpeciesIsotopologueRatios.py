import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.SpeciesIsotopeRecord import SpeciesIsotopeRecord
from pyarts.classes.io import correct_save_arguments, correct_read_arguments

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

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadSpeciesIsotopologueRatios(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveSpeciesIsotopologueRatios(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

lib.createSpeciesIsotopologueRatios.restype = c.c_void_p
lib.createSpeciesIsotopologueRatios.argtypes = []

lib.deleteSpeciesIsotopologueRatios.restype = None
lib.deleteSpeciesIsotopologueRatios.argtypes = [c.c_void_p]

lib.printSpeciesIsotopologueRatios.restype = None
lib.printSpeciesIsotopologueRatios.argtypes = [c.c_void_p]

lib.getdataSpeciesIsotopologueRatios.restype = c.POINTER(c.c_double)
lib.getdataSpeciesIsotopologueRatios.argtypes = [c.c_void_p]

lib.xmlreadSpeciesIsotopologueRatios.restype = c.c_long
lib.xmlreadSpeciesIsotopologueRatios.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveSpeciesIsotopologueRatios.restype = c.c_long
lib.xmlsaveSpeciesIsotopologueRatios.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]
