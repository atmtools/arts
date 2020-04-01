import ctypes as c
from pyarts.workspace.api import arts_api as lib


def define_array_lib(var):
    """ Set all the interface to ARTS Array

    Do not use manually

    Input:
        var:
            The baseclass of an array
    """
    global lib
    string = var.name()

    create = eval ("lib.createArrayOf{}".format(string))
    create.restype = c.c_void_p
    create.argtypes = []

    delete = eval("lib.deleteArrayOf{}".format(string))
    delete.restype = None
    delete.argtypes = [c.c_void_p]

    print = eval("lib.printArrayOf{}".format(string))
    print.restype = None
    print.argtypes = [c.c_void_p]

    size = eval("lib.sizeArrayOf{}".format(string))
    size.restype = c.c_long
    size.argtypes = [c.c_void_p]

    resize = eval("lib.resizeArrayOf{}".format(string))
    resize.restype = None
    resize.argtypes = [c.c_long, c.c_void_p]

    getelem = eval("lib.getelemArrayOf{}".format(string))
    getelem.restype = c.c_void_p
    getelem.argtypes = [c.c_long, c.c_void_p]

    xmlread = eval("lib.xmlreadArrayOf{}".format(string))
    xmlread.restype = c.c_long
    xmlread.argtypes = [c.c_void_p, c.c_char_p]

    xmlsave = eval("lib.xmlsaveArrayOf{}".format(string))
    xmlsave.restype = c.c_long
    xmlsave.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]


def array_base(var):
    """ A full Array<>-class implementation as a string

    Do not use manually

    Input:
        var:
            Class to become an Array<var> (type)

    return:
        A string to be exec to generate the class
    """

    define_array_lib(var)


    return '''class ArrayOfBASENAME:
        """ ArrayOfBASENAME from ARTS

            Getter and setter returns underlying type

            Properties:
                size:
                    Size of the array (Index)

                type:
                    Underlying type (constexpr BASENAME)

                data:
                    All the data (list of BASENAME with size-number of elements)
        """

        def __init__(self, data=None):
            if isinstance(data, c.c_void_p):
                self.__delete__ = False
                self.__data__ = data
            else:
                self.__delete__ = True
                self.__data__ = c.c_void_p(lib.createArrayOfBASENAME())
                if data is not None:
                    self.data = data

        @property
        def type(self):
            return BASENAME

        @property
        def size(self):
            """ Size of the array (Index) """
            return lib.sizeArrayOfBASENAME(self.__data__)

        @size.setter
        def size(self, size):
            lib.resizeArrayOfBASENAME(int(size), self.__data__)

        @property
        def data(self):
            """ All the data (list of BASENAME with size-number of elements)"""
            x = []
            n = self.size
            for i in range(n):
                x.append(self[i])
            return x

        @data.setter
        def data(self, val):
            if isinstance(val, Sized):
                self.size = len(val)
                n = self.size
                for i in range(n):
                    self[i] = val[i]
            else:
                raise TypeError("Expects list of values")

        def __getitem__(self, ind):
            if isinstance(ind, slice):
                return [self[i] for i in range(ind)]
            elif isinstance(ind, int) and ind < self.size:
                return self.type(c.c_void_p(lib.getelemArrayOfBASENAME(ind % self.size, self.__data__)))
            else:
                raise IndexError("Out of bounds")

        def __setitem__(self, ind, val):
            old = self[ind]
            if isinstance(ind, slice):
                if isinstance(val, Sized):
                    for i in range(len(old)):
                        old[i].set(val[i])
                else:
                    for i in range(len(old)):
                        old[i].set(val)
            elif isinstance(val, self.type):
                old.set(val)
            else:
                raise TypeError("Expect BASENAME")

        def append(self, val):
            if isinstance(val, self.type):
                self.size += 1
                self[self.size-1].set(val)
            else:
                raise TypeError("Expect BASENAME")

        def __len__(self):
            return self.size

        def print(self):
            """ Print to cout the ARTS representation of the class """
            lib.printArrayOfBASENAME(self.__data__)

        def __del__(self):
            if self.__delete__:
                lib.deleteArrayOfBASENAME(self.__data__)

        def __repr__(self):
            return repr(self.data)

        def set(self, other):
            """ Sets this class according to another python instance of itself """
            if isinstance(other, ArrayOfBASENAME):
                self.data = other.data
            else:
                raise TypeError("Expects ArrayOfBASENAME")

        @staticmethod
        def name():
            return "ArrayOfBASENAME"

        def readxml(self, file):
            """ Reads the XML file

            Input:
                file:
                    Filename to valid class-file (str)
            """
            if lib.xmlreadArrayOfBASENAME(self.__data__, correct_read_arguments(file)):
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
            if lib.xmlsaveArrayOfBASENAME(self.__data__, *correct_save_arguments(file, type, clobber)):
                raise OSError("Cannot save {}".format(file))

        def __eq__(self, other):
            n = len(self)
            if len(other) != n:
                return False

            for i in range(n):
                if self[i] != other[i]:
                    return False
            return True

        def __lt__(self, other):
            n = len(self)
            if len(other) != n:
                return False

            for i in range(n):
                if self[i] >= other[i]:
                    return False
            return True

    '''.replace("BASENAME", var.name())