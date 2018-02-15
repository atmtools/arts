
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

/* FORTRAN sends pointers, so we expect those.
 *
 * We also take 32-bit (4 byte) parameters, which matches integer (or
 * better integer*4).  When used with interface declarations in the
 * fortran source, gfortran will detect arguments of the wrong type
 * even if the user used -fdefault-integer-8, but not
 * -finteger-4-integer-8 (as that modifies the interface declaration).
 *
 * An underscore is added to the function names, which seems to be the
 * standard name mangling.
 */

void fwig_table_init_(int32_t *max_two_j, int32_t *wigner_type)
{
  wig_table_init(*max_two_j, *wigner_type);
}

void fwig_table_free_()
{
  wig_table_free();
}

void fwig_temp_init_(int32_t *max_two_j)
{
  wig_temp_init(*max_two_j);
}

#if WIGXJPF_HAVE_THREAD
void fwig_thread_temp_init_(int32_t *max_two_j)
{
  wig_temp_init(*max_two_j);
}
#endif

void fwig_temp_free_()
{
  wig_temp_free();
}

double fwig3jj_(int32_t *two_j1, int32_t *two_j2, int32_t *two_j3,
		int32_t *two_m1, int32_t *two_m2, int32_t *two_m3)
{
  double result;

  calc_3j_double(&result,
		 *two_j1, *two_j2, *two_j3,
		 *two_m1, *two_m2, *two_m3,
		 wigxjpf_global_temp);

  return result;
}

double fwig6jj_(int32_t *two_j1, int32_t *two_j2, int32_t *two_j3,
		int32_t *two_j4, int32_t *two_j5, int32_t *two_j6)
{
  double result;

  calc_6j_double(&result,
		 *two_j1, *two_j2, *two_j3,
		 *two_j4, *two_j5, *two_j6,
		 wigxjpf_global_temp);

  return result;
}

double fwig9jj_(int32_t *two_j1, int32_t *two_j2, int32_t *two_j3,
		int32_t *two_j4, int32_t *two_j5, int32_t *two_j6,
		int32_t *two_j7, int32_t *two_j8, int32_t *two_j9)
{
  double result;

  calc_9j_double(&result,
		 *two_j1, *two_j2, *two_j3,
		 *two_j4, *two_j5, *two_j6,
		 *two_j7, *two_j8, *two_j9,
		 wigxjpf_global_temp);

  return result;
}
