import ctypes as c
from collections.abc import Sized
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Tensor4 import Tensor4
from pyarts.classes.io import correct_save_arguments, correct_read_arguments
from pyarts.classes.ArrayBase import array_base


class PropagationMatrix:
    """ ARTS PropagationMatrix data

    Properties:
        data:
            The data (Tensor4)

        stokes:
            Stokes dimension (const Index)

        freqs:
            Frequency dimension (const Index)

        zeniths:
            Zenith dimension (const Index)

        azimuths:
            Azimuth dimension (const Index)

    """
    def __init__(self, data=0.0, freqs=0, stokes=1, za=1, aa=1):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createPropagationMatrix())
            self.setData(freqs, stokes, za, aa, data)

    @property
    def OK(self):
        return lib.getOKPropagationMatrix(self.__data__)

    @property
    def data(self):
        """ The data (Tensor4) """
        return Tensor4(c.c_void_p(lib.getDataPropagationMatrix(self.__data__)))

    @data.setter
    def data(self, value):
        self.data.set(value)

    @property
    def stokes(self):
        """ Stokes dimension (const Index) """
        return lib.stokesPropagationMatrix(self.__data__)

    @property
    def frequencies(self):
        """ Frequency dimension (const Index) """
        return lib.frequenciesPropagationMatrix(self.__data__)

    @property
    def zeniths(self):
        """ Zenith dimension (const Index) """
        return lib.zenithsPropagationMatrix(self.__data__)

    @property
    def azimuths(self):
        """ Azimuth dimension (const Index) """
        return lib.azimuthsPropagationMatrix(self.__data__)

    def setData(self, freqs, stokes, za, aa, data):
        """ Sets the data by reinitialization to new size """
        if lib.setPropagationMatrix(self.__data__, int(freqs), int(stokes), int(za), int(aa), float(data)):
            raise ValueError("Bad input")

    def print(self):
        """ Print to cout the ARTS representation of the class """
        if self.OK:
            lib.printPropagationMatrix(self.__data__)
        else:
            raise RuntimeError("Class is in bad state")

    def __del__(self):
        if self.__delete__:
            lib.deletePropagationMatrix(self.__data__)

    def __repr__(self):
        return repr(self.data)

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, PropagationMatrix):
            tmp = Tensor4(other.data.data)
            self.setData(other.frequencies, other.stokes, other.zeniths, other.azimuths, 0.0)
            self.data = tmp
        else:
            raise TypeError("Expects PropagationMatrix")

    @staticmethod
    def name():
        return "PropagationMatrix"

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadPropagationMatrix(self.__data__, correct_read_arguments(file)):
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
            raise RuntimeError("Class is in bad state")

        if lib.xmlsavePropagationMatrix(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, PropagationMatrix) and \
                self.frequencies == other.frequencies and \
                self.stokes == other.stokes and \
                self.zeniths == other.zeniths and \
                self.azimuths == other.azimuths and \
                self.data == other.data:
            return True
        else:
            return False

    def __bool__(self):
        return self.OK and bool(self.data)


# ArrayOfPropagationMatrix
exec(array_base(PropagationMatrix))


# ArrayOfArrayOfPropagationMatrix
exec(array_base(ArrayOfPropagationMatrix))


lib.createPropagationMatrix.restype = c.c_void_p
lib.createPropagationMatrix.argtypes = []

lib.deletePropagationMatrix.restype = None
lib.deletePropagationMatrix.argtypes = [c.c_void_p]

lib.printPropagationMatrix.restype = None
lib.printPropagationMatrix.argtypes = [c.c_void_p]

lib.xmlreadPropagationMatrix.restype = c.c_long
lib.xmlreadPropagationMatrix.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsavePropagationMatrix.restype = c.c_long
lib.xmlsavePropagationMatrix.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.stokesPropagationMatrix.restype = c.c_long
lib.stokesPropagationMatrix.argtypes = [c.c_void_p]

lib.frequenciesPropagationMatrix.restype = c.c_long
lib.frequenciesPropagationMatrix.argtypes = [c.c_void_p]

lib.zenithsPropagationMatrix.restype = c.c_long
lib.zenithsPropagationMatrix.argtypes = [c.c_void_p]

lib.azimuthsPropagationMatrix.restype = c.c_long
lib.azimuthsPropagationMatrix.argtypes = [c.c_void_p]

lib.getDataPropagationMatrix.restype = c.c_void_p
lib.getDataPropagationMatrix.argtypes = [c.c_void_p]

lib.setPropagationMatrix.restype = c.c_long
lib.setPropagationMatrix.argtypes = [c.c_void_p, c.c_long, c.c_long, c.c_long, c.c_long, c.c_double]

lib.getOKPropagationMatrix.restype = c.c_bool
lib.getOKPropagationMatrix.argtypes = [c.c_void_p]
