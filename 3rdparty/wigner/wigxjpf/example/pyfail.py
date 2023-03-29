#!/usr/bin/python
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

import pywigxjpf as wig

import sys

failure = int(sys.argv[1])

print('WIGXJPF python error handling test (%d)' % failure)
sys.stdout.flush()

# Initialize
if (failure != 1):
    wig.wig_table_init(2*100 if failure != 2 else 2*1,
                       9     if failure != 3 else -1)
if (failure != 4):
    wig.wig_temp_init(2*100 if failure != 5 else 1)

# Note that arguments are in two_j = 2*j.
val3j = wig.wig3jj(2* 10 , 2* 15 , 2* 10 ,\
                   2*(-3), 2* 12 , 2*(-9))
if failure == 6:
    val3j = wig.wig3jj(2* 1000 , 2* 1500 , 2* 1000 ,\
                       2*(-300), 2* 1200 , 2*(-900))
print('3J(10  15  10; -3  12  -9): %.12g' % val3j)
sys.stdout.flush()

val6j = wig.wig6jj(2* 10 , 2* 15 , 2* 10 ,\
                   2*  7,  2*  7 , 2*  9)
if failure == 7:
    val6j = wig.wig6jj(2* 1000 , 2* 1500 , 2* 1000 ,\
                       2*  700,  2*  700 , 2*  900)
print('6J{10  15  10;  7   7   9}: %.12g' % val6j)
sys.stdout.flush()

val9j = wig.wig9jj( 1,  2,  3,\
                    4,  6,  8, \
                    3,  6,  9 )
if failure == 8:
    val9j = wig.wig9jj( 100,  200,  300,\
                        400,  600,  800, \
                        300,  600,  900 )
print('9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %.12g' % val9j)
sys.stdout.flush()

print('Using the (numpy) array interface')
sys.stdout.flush()

val3j = wig.wig3jj_array([2* 10 , 2* 15 , 2* 10 ,\
                          2*(-3), 2* 12 , 2*(-9)])
if failure == 9:
    val3j = wig.wig3jj_array([2* 1000 , 2* 1500 , 2* 1000 ,\
                              2*(-300), 2* 1200 , 2*(-900)])
print('3J(10  15  10; -3  12  -9): %.12g' % val3j)
sys.stdout.flush()

val6j = wig.wig6jj_array([2* 10 , 2* 15 , 2* 10 ,\
                          2*  7,  2*  7 , 2*  9])
if failure == 10:
    val6j = wig.wig6jj_array([2* 1000 , 2* 1500 , 2* 1000 ,\
                              2*  700,  2*  700 , 2*  900])
print('6J{10  15  10;  7   7   9}: %.12g' % val6j)
sys.stdout.flush()

val9j = wig.wig9jj_array([1,  2,  3,\
                          4,  6,  8, \
                          3,  6,  9])
if failure == 11:
    val9j = wig.wig9jj_array([100,  200,  300,\
                              400,  600,  800, \
                              300,  600,  900])
print('9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %.12g' % val9j)
sys.stdout.flush()

# Free memory space
wig.wig_temp_free()
wig.wig_table_free()

