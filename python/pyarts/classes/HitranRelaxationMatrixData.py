import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.Tensor4 import Tensor4

from pyarts.classes.io import correct_read_arguments, correct_save_arguments

class HitranRelaxationMatrixData:
    """ ARTS HitranRelaxationMatrixData data

    Properties:
        W0pp: (Tensor4)
        B0pp: (Tensor4)
        W0rp: (Tensor4)
        B0rp: (Tensor4)
        W0qp: (Tensor4)
        B0qp: (Tensor4)
        W0pr: (Tensor4)
        B0pr: (Tensor4)
        W0rr: (Tensor4)
        B0rr: (Tensor4)
        W0qr: (Tensor4)
        B0qr: (Tensor4)
        W0pq: (Tensor4)
        B0pq: (Tensor4)
        W0rq: (Tensor4)
        B0rq: (Tensor4)
        W0qq: (Tensor4)
        B0qq: (Tensor4)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createHitranRelaxationMatrixData())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

    @property
    def W0rr(self):
        """ W0rr: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0rrHitranRelaxationMatrixData(self.__data__)))

    @W0rr.setter
    def W0rr(self, other):
        self.W0rr.set(other)

    @property
    def B0rr(self):
        """ B0rr: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0rrHitranRelaxationMatrixData(self.__data__)))

    @B0rr.setter
    def B0rr(self, other):
        self.B0rr.set(other)

    @property
    def W0rq(self):
        """ W0rq: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0rqHitranRelaxationMatrixData(self.__data__)))

    @W0rq.setter
    def W0rq(self, other):
        self.W0rq.set(other)

    @property
    def B0rq(self):
        """ B0rq: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0rqHitranRelaxationMatrixData(self.__data__)))

    @B0rq.setter
    def B0rq(self, other):
        self.B0rq.set(other)

    @property
    def W0rp(self):
        """ W0rp: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0rpHitranRelaxationMatrixData(self.__data__)))

    @W0rp.setter
    def W0rp(self, other):
        self.W0rp.set(other)

    @property
    def B0rp(self):
        """ B0rp: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0rpHitranRelaxationMatrixData(self.__data__)))

    @B0rp.setter
    def B0rp(self, other):
        self.B0rp.set(other)

    @property
    def W0qr(self):
        """ W0qr: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0qrHitranRelaxationMatrixData(self.__data__)))

    @W0qr.setter
    def W0qr(self, other):
        self.W0qr.set(other)

    @property
    def B0qr(self):
        """ B0qr: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0qrHitranRelaxationMatrixData(self.__data__)))

    @B0qr.setter
    def B0qr(self, other):
        self.B0qr.set(other)

    @property
    def W0qq(self):
        """ W0qq: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0qqHitranRelaxationMatrixData(self.__data__)))

    @W0qq.setter
    def W0qq(self, other):
        self.W0qq.set(other)

    @property
    def B0qq(self):
        """ B0qq: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0qqHitranRelaxationMatrixData(self.__data__)))

    @B0qq.setter
    def B0qq(self, other):
        self.B0qq.set(other)

    @property
    def W0qp(self):
        """ W0qp: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0qpHitranRelaxationMatrixData(self.__data__)))

    @W0qp.setter
    def W0qp(self, other):
        self.W0qp.set(other)

    @property
    def B0qp(self):
        """ B0qp: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0qpHitranRelaxationMatrixData(self.__data__)))

    @B0qp.setter
    def B0qp(self, other):
        self.B0qp.set(other)

    @property
    def W0pr(self):
        """ W0pr: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0prHitranRelaxationMatrixData(self.__data__)))

    @W0pr.setter
    def W0pr(self, other):
        self.W0pr.set(other)

    @property
    def B0pr(self):
        """ B0pr: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0prHitranRelaxationMatrixData(self.__data__)))

    @B0pr.setter
    def B0pr(self, other):
        self.B0pr.set(other)

    @property
    def W0pq(self):
        """ W0pq: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0pqHitranRelaxationMatrixData(self.__data__)))

    @W0pq.setter
    def W0pq(self, other):
        self.W0pq.set(other)

    @property
    def B0pq(self):
        """ B0pq: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0pqHitranRelaxationMatrixData(self.__data__)))

    @B0pq.setter
    def B0pq(self, other):
        self.B0pq.set(other)

    @property
    def W0pp(self):
        """ W0pp: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getW0ppHitranRelaxationMatrixData(self.__data__)))

    @W0pp.setter
    def W0pp(self, other):
        self.W0pp.set(other)

    @property
    def B0pp(self):
        """ B0pp: (Tensor4) """
        return Tensor4(c.c_void_p(lib.getB0ppHitranRelaxationMatrixData(self.__data__)))

    @B0pp.setter
    def B0pp(self, other):
        self.B0pp.set(other)

    @staticmethod
    def name():
        return "HitranRelaxationMatrixData"

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printHitranRelaxationMatrixData(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteHitranRelaxationMatrixData(self.__data__)

    def __repr__(self):
        return "ARTS HitranRelaxationMatrixData"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, HitranRelaxationMatrixData):
            self.W0rr = other.W0rr
            self.B0rr = other.B0rr
            self.W0rq = other.W0rq
            self.B0rq = other.B0rq
            self.W0rp = other.W0rp
            self.B0rp = other.B0rp
            self.W0qr = other.W0qr
            self.B0qr = other.B0qr
            self.W0qq = other.W0qq
            self.B0qq = other.B0qq
            self.W0qp = other.W0qp
            self.B0qp = other.B0qp
            self.W0pr = other.W0pr
            self.B0pr = other.B0pr
            self.W0pq = other.W0pq
            self.B0pq = other.B0pq
            self.W0pp = other.W0pp
            self.B0pp = other.B0pp
        else:
            raise TypeError("Expects HitranRelaxationMatrixData")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadHitranRelaxationMatrixData(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveHitranRelaxationMatrixData(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))

    def __eq__(self, other):
        if isinstance(other, HitranRelaxationMatrixData) and \
                self.W0rr == other.W0rr and \
                self.B0rr == other.B0rr and \
                self.W0rq == other.W0rq and \
                self.B0rq == other.B0rq and \
                self.W0rp == other.W0rp and \
                self.B0rp == other.B0rp and \
                self.W0qr == other.W0qr and \
                self.B0qr == other.B0qr and \
                self.W0qq == other.W0qq and \
                self.B0qq == other.B0qq and \
                self.W0qp == other.W0qp and \
                self.B0qp == other.B0qp and \
                self.W0pr == other.W0pr and \
                self.B0pr == other.B0pr and \
                self.W0pq == other.W0pq and \
                self.B0pq == other.B0pq and \
                self.W0pp == other.W0pp and \
                self.B0pp == other.B0pp:
            return True
        else:
            return False


lib.createHitranRelaxationMatrixData.restype = c.c_void_p
lib.createHitranRelaxationMatrixData.argtypes = []

lib.deleteHitranRelaxationMatrixData.restype = None
lib.deleteHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.printHitranRelaxationMatrixData.restype = None
lib.printHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.xmlreadHitranRelaxationMatrixData.restype = c.c_long
lib.xmlreadHitranRelaxationMatrixData.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveHitranRelaxationMatrixData.restype = c.c_long
lib.xmlsaveHitranRelaxationMatrixData.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getW0rrHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0rrHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0rrHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0rrHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0rqHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0rqHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0rqHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0rqHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0rpHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0rpHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0rpHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0rpHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0qrHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0qrHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0qrHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0qrHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0qqHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0qqHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0qqHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0qqHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0qpHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0qpHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0qpHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0qpHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0prHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0prHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0prHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0prHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0pqHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0pqHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0pqHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0pqHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getW0ppHitranRelaxationMatrixData.restype = c.c_void_p
lib.getW0ppHitranRelaxationMatrixData.argtypes = [c.c_void_p]

lib.getB0ppHitranRelaxationMatrixData.restype = c.c_void_p
lib.getB0ppHitranRelaxationMatrixData.argtypes = [c.c_void_p]
