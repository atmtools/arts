import ctypes as c
from pyarts.workspace.api import arts_api as lib


def define_array_lib(var):
    """ Ugly way to return and set all the interface to ARTS Array

    Do not use manually

    Input:
        var:
            The baseclass of an array

    Output:
        Creatoin functions for the array type
        create, delete, print, size, resize, getelem, xmlread, xmlsave, var, var.name() (tuple)
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

    return (create, delete, print, size, resize, getelem, xmlread, xmlsave, var, string)

def array_base(var):
    """ A full Array<>-class implementation as a string

    Do not use manually

    Input:
        var:
            Class to become an Array<var> (type)

    return:
        A string to be exec to generate the class
    """
    return '''class ArrayOfBASENAME:
        __funcreate, __fundel, __funprint, __funsize, __funresize, \\
            __fungetelem, __funxmlread, __funxmlsave, __type, __strtype = \\
                define_array_lib(BASENAME)

        __doc__ = """ Array of {} from ARTS

            Getter and setter returns underlying type

            Properties:
                size:
                    Size of the array (Index)

                type:
                    Underlying type (const type)

                data:
                    All the data (list of type with size elements)
        """.format(__strtype)

        def __init__(self, data=None):
            if isinstance(data, c.c_void_p):
                self.__delete__ = False
                self.__data__ = data
            elif data is None:
                self.__delete__ = True
                self.__data__ = c.c_void_p(__class__.__funcreate())
            else:
                raise TypeError("Invalid initialization")

        @property
        def type(self):
            return self.__type

        @property
        def size(self):
            """ Size of the array (Index) """
            return __class__.__funsize(self.__data__)

        @size.setter
        def size(self, size):
            size = int(size)
            __class__.__funresize(size, self.__data__)

        @property
        def data(self):
            """ All the data (list of type with size elements)"""
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
            ind = int(ind)
            if ind >= 0 and ind < self.size:
                return self.type(c.c_void_p(__class__.__fungetelem(ind, self.__data__)))
            else:
                raise IndexError("Out of bounds")

        def __setitem__(self, ind, val):
            old = self[ind]
            if isinstance(val, self.type):
                old.set(val)
            else:
                raise TypeError("Expect {}".format(__class__.__strtype))

        def print(self):
            """ Print to cout the ARTS representation of the class """
            __class__.__funprint(self.__data__)

        def __del__(self):
            if self.__delete__:
                __class__.__fundel(self.__data__)

        def __repr__(self):
            return "ARTS {} with {} elements".format(__class__.__name__, self.size)

        def set(self, other):
            """ Sets this class according to another python instance of itself """
            if isinstance(other, type(self)):
                pass
            else:
                raise TypeError("Expects {}".format(self.name()))

        @staticmethod
        def name():
            return __class__.__name__

        def readxml(self, file):
            """ Reads the XML file

            Input:
                file:
                    Filename to valid class-file (str)
            """
            if __class__.__funxmlread(self.__data__, correct_read_arguments(file)):
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
            if __class__.__funxmlsave(self.__data__, *correct_save_arguments(file, type, clobber)):
                raise OSError("Cannot save {}".format(file))
    '''.replace("BASENAME", var.name())