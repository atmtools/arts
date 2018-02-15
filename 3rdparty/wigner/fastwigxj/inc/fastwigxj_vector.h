
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of FASTWIGXJ.
 *
 *  FASTWIGXJ is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  FASTWIGXJ is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with FASTWIGXJ.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#ifndef __FASTWIGXJ_VECTOR_H__
#define __FASTWIGXJ_VECTOR_H__

#include <stdint.h>

typedef long long int  v2di __attribute__ ((vector_size (16)));
typedef int            v4si __attribute__ ((vector_size (16)));

typedef long long int  v4di __attribute__ ((vector_size (32)));

typedef v4si v4su32;
typedef v4di v4su64;

#endif/*__FASTWIGXJ_VECTOR_H__*/
