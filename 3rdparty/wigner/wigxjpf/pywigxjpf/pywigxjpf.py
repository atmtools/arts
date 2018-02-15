#
# Copyright 2015 Christian Forssen
#
#  This file is part of WIGXJPF.
#
#  WIGXJPF is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  WIGXJPF is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with WIGXJPF.  If not, see
#  <http://www.gnu.org/licenses/>.
#

"""
Filename: pywigxjpf.py
Python interface to the wigxjpf (dynamic) library with the use of ctypes.
Defines seven functions:

    (1) wig_table_init(max_two_j,wigner_type)
    (2) wig_table_free()
    (3) wig_temp_init(max_two_j)
    (4) wig_temp_free()
    (5) wig3jj(jj1,jj2,jj3,
               mm1,mm2,mm3)
    (6) wig6jj(jj1,jj2,jj3,
               jj4,jj5,jj6)
    (7) wig9jj(jj1,jj2,jj3,
               jj4,jj5,jj6,
               jj7,jj8,jj9)
"""

__version__ = 1.7
__author__ = 'C. Forssen'
__all__ = ['pywigxjpf']

import numpy as nm
import ctypes as ct
import os

# Load the library as _libwigxjpf.
_path = os.path.join(os.path.dirname(__file__), '../lib')
_libwigxjpf = nm.ctypeslib.load_library('libwigxjpf_shared', _path)

# Define argument and result types
_libwigxjpf.wig_table_init.argtypes = [ct.c_int]*2
_libwigxjpf.wig_table_init.restype = ct.c_void_p
_libwigxjpf.wig_table_free.argtypes = []
_libwigxjpf.wig_table_free.restype = ct.c_void_p
_libwigxjpf.wig_temp_init.argtypes = [ct.c_int]
_libwigxjpf.wig_temp_init.restype = ct.c_void_p
_libwigxjpf.wig_temp_free.argtypes = []
_libwigxjpf.wig_temp_free.restype = ct.c_void_p

_libwigxjpf.wig3jj.argtypes = [ct.c_int]*6
_libwigxjpf.wig3jj.restype = ct.c_double
_libwigxjpf.wig6jj.argtypes = [ct.c_int]*6
_libwigxjpf.wig6jj.restype = ct.c_double
_libwigxjpf.wig9jj.argtypes = [ct.c_int]*9
_libwigxjpf.wig9jj.restype = ct.c_double


def wig_table_init(max_two_j,wigner_type=9):
    """Setup of precalculated tables of prime-factorised factorials.
   The tables must be large enough to handle the largest symbol that
   shall be evaluated.  For multi-threaded programs, this is done once
   (globally).  The limit max_two_j is chosen as the largest two_j in
   any symbol calculated, and wigner_type as the (largest of) 3, 6, or
   9 for 3j, 6j, or 9j symbols respectively.  Note that giving
   e.g. max_two_j = 2*100 and wigner_type = 9 only requires < 400 kB
   of memory, such that there is nowadays little point in keeping this
   limit very tight.  The evaluation routines will only use the
   necessary parts of each prime-factorisation list, so there is only
   a slight cache penalty for using larger tables than necessary.

    ARGUMENT(S):
        An integer.

    RESULT(S):
        No output
    """
    return _libwigxjpf.wig_table_init(int(max_two_j),int(wigner_type))

def wig_table_free():
    """Free above table memory.

    ARGUMENT(S):
        None

    RESULT(S):
        No output
    """
    return _libwigxjpf.wig_table_free()

def wig_temp_init(max_two_j):
    """Allocation of temporary storage for the evaluation routines.
   Enough space must be allocated for the largest symbol that shall be
   evaluated.  It also depends on the tables above, which must be
   allocated first.  As above, there is little to gain by giving this
   limit tightly, giving e.g. max_two_j = 2*100 uses less than 50 kB.

    ARGUMENT(S):
        An integer.

    RESULT(S):
        No output
    """
    return _libwigxjpf.wig_temp_init(int(max_two_j))

def wig_temp_free():
    """Free above temporary memory.

    ARGUMENT(S):
        None

    RESULT(S):
        No output
    """
    return _libwigxjpf.wig_temp_free()

def wig3jj(b):
    """Compute the Wigner 3j symbol
    (j1 j2 j3
     m1 m2 m3)
     where j1 = jj1 / 2; m1 = mm1 / 2; etc


    ARGUMENT(S): 
        A list or array consisting of 6 integers.  
        b = [jj1,jj2,jj3,mm1,mm2,mm3]
        Note that the arguments are to be given as integers, with
        twice the numeric value (this is what jj tries to indicate).
        I.e. half-integer arguments will be passed as odd integers.

    SETUP(S):
        Remember to first call wig_table_init and wig_temp_init.

    RESULT(S):
        Wigner 3j symbol in double precision.

    EXAMPLE(S):
    # Evaluate the 3j symbol (3/2  3/2  1 
    #                         3/2 -1/2 -1)
    >>> wig.wig_table_init(100)        # only once
    >>> wig.wig_temp_init(100)         # only once
    >>> wig3jj([3, 3, 2,   3, -1, -2]) 
    -0.316227766016838
    """
    jj1,jj2,jj3,mm1,mm2,mm3 = nm.asarray(b,dtype=nm.intc)
    return _libwigxjpf.wig3jj(jj1,jj2,jj3,mm1,mm2,mm3)

def wig6jj(b):
    """Compute the Wigner 6j symbol
    (j1 j2 j3
     j4 j5 j6), 
     where j1 = jj1 / 2; j2 = jj2 / 2; etc

    ARGUMENT(S): 
        A list or array consisting of 6 integers.  
        b = [jj1,jj2,jj3,jj4,jj5,jj6]
        Note that the arguments are to be given as integers, with
        twice the numeric value (this is what jj tries to indicate).
        I.e. half-integer arguments will be passed as odd integers.

    SETUP(S):
        Remember to first call wig_table_init and wig_temp_init.

    RESULT(S):
        Wigner 6j symbol in double precision.

    EXAMPLE(S): 
    # Evaluate the 6j symbol (2  2  1 
    #                         2  1  1)
    >>> wig.wig_table_init(100)        # only once
    >>> wig.wig_temp_init(100)         # only once
    >>> wig6jj([4, 4, 2,   4, 2, 2]) 
    -0.1000000000000
    """
    jj1,jj2,jj3,jj4,jj5,jj6 = nm.asarray(b,dtype=nm.intc)
    return _libwigxjpf.wig6jj(jj1,jj2,jj3,jj4,jj5,jj6)

def wig9jj(b):
    """Compute the Wigner 9j symbol
    (j1 j2 j3
     j4 j5 j6
     j7 j8 j9
     where j1 = jj1 / 2; j2 = jj2 / 2; etc

    ARGUMENT(S): 
        A list or array consisting of 9 integers.  
        b = [jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9]
        Note that the arguments are to be given as integers, with
        twice the numeric value (this is what jj tries to indicate).
        I.e. half-integer arguments will be passed as odd integers.

    SETUP(S):
        Remember to first call wig_table_init and wig_temp_init.

    RESULT(S):
        Wigner 9j symbol in double precision.

    EXAMPLE(S): 
    # Evaluate the 9j symbol (20  20  40 
    #                         20  20  40
    #                         20  20  40)
    >>> wig.wig_table_init(100)        # only once
    >>> wig.wig_temp_init(100)         # only once
    >>> wig9jj([40, 40, 80, 40, 40, 80, 40, 40, 80]) 
    0.000111633558388972
    """
    jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9 = nm.asarray(b,dtype=nm.intc)
    return _libwigxjpf.wig9jj(jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9)

if __name__ == '__main__':
    print 'WIGXJPF python test facility'

    # Initialize
    wig_table_init(2*100, 9)
    wig_temp_init(2*100)

    # Note that arguments are in two_j = 2*j.
    val3j = wig3jj([2* 10 , 2* 15 , 2* 10 ,\
                        2*(-3), 2* 12 , 2*(-9)])
    print '3J(10  15  10; -3  12  -9):', val3j

    val6j = wig6jj([2* 10 , 2* 15 , 2* 10 ,\
                        2*  7,  2*  7 , 2*  9 ])
    print "6J{10  15  10;  7   7   9}:", val6j

    val9j = wig9jj( [1,  2,  3,\
                         4,  6,  8, \
                         3,  6,  9] )
    print "9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}:", val9j

    # Free memory space
    wig_temp_free()
    wig_table_free()
