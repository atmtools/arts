import ctypes as c
from pyarts.workspace.api import arts_api as lib

from pyarts.classes.io import correct_save_arguments, correct_read_arguments


class Timer:
    """ ARTS Timer data

    Depends on if the system defined TIMER_SUPPORT

    Converts between clock_t and Index

    Properties:
        running:
            (bool)

        finished:
            (bool)

        cputime_start:
            (list of 4 Index)

        realtime_start:
            (Index)

        cputime_end:
            (list of 4 Index)

        realtime_end:
            (Index)

        supported:
            Is TIMER_SUPPORT OK? (constexpr bool)

        tick:
            Tick time of system (constexpr Index)
    """

    def __init__(self, data=None):
        if isinstance(data, c.c_void_p):
            self.__delete__ = False
            self.__data__ = data
        else:
            self.__delete__ = True
            self.__data__ = c.c_void_p(lib.createTimer())
            if data is not None:
                raise RuntimeError("Only supports void initialization")

        if not self.supported:
            raise RuntimeWarning("System oes not support Timer")

    @staticmethod
    def name():
        return "Timer"

    @property
    def supported(self):
        """ Is TIMER_SUPPORT OK? (constexpr bool) """
        return lib.supportTimer()

    @property
    def tick(self):
        """ Tick time of system (constexpr Index)"""
        return lib.tickTimer()

    @property
    def realtime_end(self):
        """ (Index) """
        return lib.getrealtime_endTimer(self.__data__)

    @realtime_end.setter
    def realtime_end(self, x):
        lib.setrealtime_endTimer(self.__data__, int(x))

    @property
    def realtime_start(self):
        """ (Index) """
        return lib.getrealtime_startTimer(self.__data__)

    @realtime_start.setter
    def realtime_start(self, x):
        lib.setrealtime_startTimer(self.__data__, int(x))

    @property
    def cputime_end(self):
        """ (list of 4 Index) """
        return lib.getcputime_end_utimeTimer(self.__data__),  \
               lib.getcputime_end_stimeTimer(self.__data__), \
               lib.getcputime_end_cutimeTimer(self.__data__), \
               lib.getcputime_end_cstimeTimer(self.__data__)

    @cputime_end.setter
    def cputime_end(self, x):
        lib.getcputime_end_utimeTimer(self.__data__, int(x[0]))
        lib.getcputime_end_stimeTimer(self.__data__, int(x[1]))
        lib.getcputime_end_cutimeTimer(self.__data__, int(x[2]))
        lib.getcputime_end_cstimeTimer(self.__data__, int(x[3]))

    @property
    def cputime_start(self):
        """ (list of 4 Index) """
        return lib.getcputime_start_utimeTimer(self.__data__),  \
               lib.getcputime_start_stimeTimer(self.__data__), \
               lib.getcputime_start_cutimeTimer(self.__data__), \
               lib.getcputime_start_cstimeTimer(self.__data__)

    @cputime_start.setter
    def cputime_start(self, x):
        lib.getcputime_start_utimeTimer(self.__data__, int(x[0]))
        lib.getcputime_start_stimeTimer(self.__data__, int(x[1]))
        lib.getcputime_start_cutimeTimer(self.__data__, int(x[2]))
        lib.getcputime_start_cstimeTimer(self.__data__, int(x[3]))

    @property
    def running(self):
        """ (bool) """
        return lib.getrunningTimer(self.__data__)

    @running.setter
    def running(self, x):
        lib.setrunningTimer(self.__data__, bool(x))

    @property
    def finished(self):
        """ (bool) """
        return lib.getfinishedTimer(self.__data__)

    @finished.setter
    def finished(self, x):
        lib.setfinishedTimer(self.__data__, bool(x))

    def print(self):
        """ Print to cout the ARTS representation of the class """
        lib.printTimer(self.__data__)

    def __del__(self):
        if self.__delete__:
            lib.deleteTimer(self.__data__)

    def __repr__(self):
        return "ARTS Timer"

    def set(self, other):
        """ Sets this class according to another python instance of itself """
        if isinstance(other, Timer):
            raise RuntimeWarning("Cannot set Timer, remains constant")
        else:
            raise TypeError("Expects Timer")

    def readxml(self, file):
        """ Reads the XML file

        Input:
            file:
                Filename to valid class-file (str)
        """
        if lib.xmlreadTimer(self.__data__, correct_read_arguments(file)):
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
        if lib.xmlsaveTimer(self.__data__, *correct_save_arguments(file, type, clobber)):
            raise OSError("Cannot save {}".format(file))


lib.createTimer.restype = c.c_void_p
lib.createTimer.argtypes = []

lib.deleteTimer.restype = None
lib.deleteTimer.argtypes = [c.c_void_p]

lib.printTimer.restype = None
lib.printTimer.argtypes = [c.c_void_p]

lib.xmlreadTimer.restype = c.c_long
lib.xmlreadTimer.argtypes = [c.c_void_p, c.c_char_p]

lib.xmlsaveTimer.restype = c.c_long
lib.xmlsaveTimer.argtypes = [c.c_void_p, c.c_char_p, c.c_long, c.c_long]

lib.getrunningTimer.restype = c.c_bool
lib.getrunningTimer.argtypes = [c.c_void_p]

lib.getfinishedTimer.restype = c.c_bool
lib.getfinishedTimer.argtypes = [c.c_void_p]

lib.getcputime_start_utimeTimer.restype = c.c_long
lib.getcputime_start_utimeTimer.argtypes = [c.c_void_p]

lib.getcputime_start_stimeTimer.restype = c.c_long
lib.getcputime_start_stimeTimer.argtypes = [c.c_void_p]

lib.getcputime_start_cutimeTimer.restype = c.c_long
lib.getcputime_start_cutimeTimer.argtypes = [c.c_void_p]

lib.getcputime_start_cstimeTimer.restype = c.c_long
lib.getcputime_start_cstimeTimer.argtypes = [c.c_void_p]

lib.getrealtime_startTimer.restype = c.c_long
lib.getrealtime_startTimer.argtypes = [c.c_void_p]

lib.getcputime_end_utimeTimer.restype = c.c_long
lib.getcputime_end_utimeTimer.argtypes = [c.c_void_p]

lib.getcputime_end_stimeTimer.restype = c.c_long
lib.getcputime_end_stimeTimer.argtypes = [c.c_void_p]

lib.getcputime_end_cutimeTimer.restype = c.c_long
lib.getcputime_end_cutimeTimer.argtypes = [c.c_void_p]

lib.getcputime_end_cstimeTimer.restype = c.c_long
lib.getcputime_end_cstimeTimer.argtypes = [c.c_void_p]

lib.getrealtime_endTimer.restype = c.c_long
lib.getrealtime_endTimer.argtypes = [c.c_void_p]

lib.setrunningTimer.restype = None
lib.setrunningTimer.argtypes = [c.c_void_p, c.c_bool]

lib.setfinishedTimer.restype = None
lib.setfinishedTimer.argtypes = [c.c_void_p, c.c_bool]

lib.setcputime_start_utimeTimer.restype = None
lib.setcputime_start_utimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setcputime_start_stimeTimer.restype = None
lib.setcputime_start_stimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setcputime_start_cutimeTimer.restype = None
lib.setcputime_start_cutimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setcputime_start_cstimeTimer.restype = None
lib.setcputime_start_cstimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setrealtime_startTimer.restype = None
lib.setrealtime_startTimer.argtypes = [c.c_void_p, c.c_long]

lib.setcputime_end_utimeTimer.restype = None
lib.setcputime_end_utimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setcputime_end_stimeTimer.restype = None
lib.setcputime_end_stimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setcputime_end_cutimeTimer.restype = None
lib.setcputime_end_cutimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setcputime_end_cstimeTimer.restype = None
lib.setcputime_end_cstimeTimer.argtypes = [c.c_void_p, c.c_long]

lib.setrealtime_endTimer.restype = None
lib.setrealtime_endTimer.argtypes = [c.c_void_p, c.c_long]

lib.supportTimer.restype = c.c_bool
lib.supportTimer.argtypes = []

lib.tickTimer.restype = c.c_long
lib.tickTimer.argtypes = []
