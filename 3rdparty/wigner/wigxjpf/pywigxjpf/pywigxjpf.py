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
Python interface to the wigxjpf (dynamic) library with the use of cffi.
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

Note that the arguments are to be given as integers, with
twice the numeric value (this is what jj tries to indicate).
I.e. half-integer arguments will be passed as odd integers.

The two init functions must be called before evaluating any
symbol.

In addition, interfaces that take an array with the arguments
are also provided.  They match the previously provided
interface, but are renamed:

    (8)  wig3jj_array([jj1,jj2,jj3,
                       mm1,mm2,mm3])
    (9)  wig6jj_array([jj1,jj2,jj3,
                       jj4,jj5,jj6])
    (10) wig9jj_array([jj1,jj2,jj3,
                       jj4,jj5,jj6,
                       jj7,jj8,jj9])
"""

__version__ = '1.11'
__author__ = 'C. Forssen and H.T. Johansson'
__all__ = ['pywigxjpf']

import numpy as np

try:
    from pywigxjpf_ffi import ffi, lib
except ImportError:
    from pywigxjpf.pywigxjpf_ffi import ffi, lib

wig_table_init = lib.wig_table_init
wig_table_free = lib.wig_table_free

wig_temp_init  = lib.wig_temp_init
wig_temp_free  = lib.wig_temp_free

wig3jj = lib.wig3jj
wig6jj = lib.wig6jj
wig9jj = lib.wig9jj

def wig3jj_array(b):
    jj1,jj2,jj3,mm1,mm2,mm3 = np.asarray(b,dtype=np.intc)
    return lib.wig3jj(jj1,jj2,jj3,mm1,mm2,mm3)

def wig6jj_array(b):
    jj1,jj2,jj3,jj4,jj5,jj6 = np.asarray(b,dtype=np.intc)
    return lib.wig6jj(jj1,jj2,jj3,jj4,jj5,jj6)

def wig9jj_array(b):
    jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9 = np.asarray(b,dtype=np.intc)
    return lib.wig9jj(jj1,jj2,jj3,jj4,jj5,jj6,jj7,jj8,jj9)

if __name__ == '__main__':
    print('WIGXJPF python test facility')

    # Initialize
    wig_table_init(2*100, 9)
    wig_temp_init(2*100)

    # Note that arguments are in two_j = 2*j.

    # First use the CFFI style direct calls:
    val3j = wig3jj(2* 10 , 2* 15 , 2* 10 ,\
                   2*(-3), 2* 12 , 2*(-9))
    print('3J(10  15  10; -3  12  -9): %.12g' % val3j)

    val6j = wig6jj(2* 10 , 2* 15 , 2* 10 ,\
                   2*  7,  2*  7 , 2*  9 )
    print('6J{10  15  10;  7   7   9}: %.12g' % val6j)

    val9j = wig9jj( 1,  2,  3, \
                    4,  6,  8, \
                    3,  6,  9 )
    print('9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %.12g' % val9j)

    # Then the wrapped calls, taking an array as argument:
    val3j = wig3jj_array([2* 10 , 2* 15 , 2* 10 ,\
                          2*(-3), 2* 12 , 2*(-9)])
    print('3J(10  15  10; -3  12  -9): %.12g' % val3j)

    val6j = wig6jj_array([2* 10 , 2* 15 , 2* 10 ,\
                          2*  7,  2*  7 , 2*  9 ])
    print('6J{10  15  10;  7   7   9}: %.12g' % val6j)

    val9j = wig9jj_array( [1,  2,  3, \
                           4,  6,  8, \
                           3,  6,  9] )
    print('9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %.12g' % val9j)

    # Free memory space
    wig_temp_free()
    wig_table_free()
