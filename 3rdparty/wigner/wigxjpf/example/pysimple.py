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

print('WIGXJPF python test program')

# Initialize
wig.wig_table_init(2*100,9)
wig.wig_temp_init(2*100)

# Note that arguments are in two_j = 2*j.
val3j = wig.wig3jj(2* 10 , 2* 15 , 2* 10 ,\
                   2*(-3), 2* 12 , 2*(-9))
print('3J(10  15  10; -3  12  -9): %.12g' % val3j)

val6j = wig.wig6jj(2* 10 , 2* 15 , 2* 10 ,\
                   2*  7,  2*  7 , 2*  9)
print('6J{10  15  10;  7   7   9}: %.12g' % val6j)

val9j = wig.wig9jj( 1,  2,  3,\
                    4,  6,  8, \
                    3,  6,  9 )
print('9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %.12g' % val9j)

print('Using the (numpy) array interface')

val3j = wig.wig3jj_array([2* 10 , 2* 15 , 2* 10 ,\
                          2*(-3), 2* 12 , 2*(-9)])
print('3J(10  15  10; -3  12  -9): %.12g' % val3j)

val6j = wig.wig6jj_array([2* 10 , 2* 15 , 2* 10 ,\
                          2*  7,  2*  7 , 2*  9])
print('6J{10  15  10;  7   7   9}: %.12g' % val6j)

val9j = wig.wig9jj_array([1,  2,  3,\
                          4,  6,  8, \
                          3,  6,  9])
print('9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %.12g' % val9j)

# Free memory space
wig.wig_temp_free()
wig.wig_table_free()

