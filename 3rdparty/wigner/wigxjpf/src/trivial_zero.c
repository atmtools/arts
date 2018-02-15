
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

#include <stdio.h>

#include "wigxjpf.h"

/* If a > b + c => (b + c) - a < 0, i.e. sign bit in result.
 * If a < |b - c| => a - (b - c) < 0 or a - (c - b) < 0, i.e.
 * sign bit in one of those two results.
 */

#define COLLECT_NEGATIVE(two_j1,two_j2,two_j3) do {	\
    collect_sign |= (two_j1) | (two_j2) | (two_j3);	\
  } while (0)

#define COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1,two_j2,two_j3) do {	\
    collect_odd |= (two_j1) + (two_j2) + (two_j3);			\
    collect_sign |= (((two_j2) + (two_j3)) - (two_j1));			\
    collect_sign |= ((two_j1) - ((two_j2) - (two_j3)));			\
    collect_sign |= ((two_j1) - ((two_j3) - (two_j2)));			\
  } while (0)

#define COLLECT_ABS_M_WITHIN_J(two_m,two_j) do {	\
    collect_odd |= ((two_m) + (two_j));			\
    collect_sign |= ((two_j) - (two_m));		\
    collect_sign |= ((two_j) + (two_m));		\
  } while (0)

int trivial_zero_3j(int two_j1, int two_j2, int two_j3,
		    int two_m1, int two_m2, int two_m3)
{
  int collect_sign = 0;
  int collect_odd = 0;

  COLLECT_NEGATIVE(two_j1, two_j2, two_j3);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j2, two_j3);

  COLLECT_ABS_M_WITHIN_J(two_m1, two_j1);
  COLLECT_ABS_M_WITHIN_J(two_m2, two_j2);
  COLLECT_ABS_M_WITHIN_J(two_m3, two_j3);

  return ((two_m1 + two_m2 + two_m3) |
	  (collect_sign & (1 << (sizeof (int) * 8 - 1))) |
	  (collect_odd & 1));
}

int trivial_zero_6j(int two_j1, int two_j2, int two_j3,
		    int two_j4, int two_j5, int two_j6)
{
  int collect_sign = 0;
  int collect_odd = 0;

  COLLECT_NEGATIVE(two_j1, two_j2, two_j3);
  COLLECT_NEGATIVE(two_j4, two_j5, two_j6);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j2, two_j3);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j5, two_j6);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j4, two_j2, two_j6);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j4, two_j5, two_j3);

  return ((collect_sign & (1 << (sizeof (int) * 8 - 1))) |
	  (collect_odd & 1));
}

int trivial_zero_9j(int two_j1, int two_j2, int two_j3,
		    int two_j4, int two_j5, int two_j6,
		    int two_j7, int two_j8, int two_j9)
{
  int collect_sign = 0;
  int collect_odd = 0;

  COLLECT_NEGATIVE(two_j1, two_j2, two_j3);
  COLLECT_NEGATIVE(two_j4, two_j5, two_j6);
  COLLECT_NEGATIVE(two_j7, two_j8, two_j9);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j2, two_j3);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j4, two_j5, two_j6);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j7, two_j8, two_j9);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j1, two_j4, two_j7);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j2, two_j5, two_j8);
  COLLECT_TRIANGLE_TRIVIAL_ZERO(two_j3, two_j6, two_j9);

  return ((collect_sign & (1 << (sizeof (int) * 8 - 1))) |
	  (collect_odd & 1));
}
