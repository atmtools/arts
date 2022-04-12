import ctypes as c
from pyarts.workspace.api import arts_api as lib


def correct_read_arguments(filename):
    """ Return ARTS-like arguments for reading file IO

    Do not use manually

    Input:
        filename:
            Filename (str)

    Output:
        Filename (bytes)
    """

    if isinstance(filename, str):
        return filename.encode("ascii")
    else:
        raise TypeError("Expects str")


def correct_save_arguments(filename, type, clobber):
    """ Return ARTS-like arguments for writing file IO

    Do not use manually

    Input:
        filename:
            Filename (str)

        type:
            Filetype.  Valid are "ascii", "zascii", or "binary" (str)

        clobber:
            Truth value wheether a file should be clobbered (any boolean)

    Output:
        (Index, Index) input for xmlsave*-calls
    """
    if not isinstance(type, str):
        raise TypeError("Expects str")

    typefileascii = type.encode("ascii")

    typeval = lib.string2filetypeindex(typefileascii)
    if typeval < 0:
        raise ValueError("Bad filetype {}".format(type))

    return (correct_read_arguments(filename), typeval, 1 if clobber else 0)


lib.string2filetypeindex.restype = c.c_long
lib.string2filetypeindex.argtypes = [c.c_char_p]
