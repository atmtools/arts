
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of WIGXJPF.
 *
 *  WIGXJPF is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  WIGXJPF is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with WIGXJPF.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#include "c_wrap.h"

#include "wigxjpf_quadmath.h"

void wig3jj_float128(__float128 *result,
		     int two_j1, int two_j2, int two_j3,
		     int two_m1, int two_m2, int two_m3)
{
  calc_3j_float128(result,
		 two_j1, two_j2, two_j3,
		 two_m1, two_m2, two_m3,
		 wigxjpf_global_temp);
}

void wig6jj_float128(__float128 *result,
		     int two_j1, int two_j2, int two_j3,
		     int two_j4, int two_j5, int two_j6)
{
  calc_6j_float128(result,
		 two_j1, two_j2, two_j3,
		 two_j4, two_j5, two_j6,
		 wigxjpf_global_temp);
}

void wig9jj_float128(__float128 *result,
		     int two_j1, int two_j2, int two_j3,
		     int two_j4, int two_j5, int two_j6,
		     int two_j7, int two_j8, int two_j9)
{
  calc_9j_float128(result,
		 two_j1, two_j2, two_j3,
		 two_j4, two_j5, two_j6,
		 two_j7, two_j8, two_j9,
		 wigxjpf_global_temp);
}
