
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

/* This small program is for compile-time test of availability of long
 * double and __float128, threading support,
 */

#include <stdio.h>

#if TEST_LONG_DOUBLE
#include <math.h>
#endif
#if TEST_FLOAT128
#include "quadmath.h"
#endif
#if TEST_THREAD
__thread int global = 0;
#endif

int main()
{
#if TEST_LONG_DOUBLE
  long double a;
  long double b = 1.14, c = 2.00159;
  long double d;

  a = b + c;
  d = ldexpl(c, 5);

  printf ("%Lf\n", a);
  printf ("#define WIGXJPF_IMPL_LONG_DOUBLE 1\n");
#endif
#if TEST_FLOAT128
  __float128 a;
  __float128 b = 1.14, c = 2.00159;
  char s[64];

  a = b + c;

  quadmath_snprintf(s, sizeof(s), "%Qf", a);

  printf ("%s\n", s);
  printf ("#define WIGXJPF_IMPL_FLOAT128 1\n");
#endif
#if TEST_THREAD
  global = 2;
  printf ("#define WIGXJPF_HAVE_THREAD 1\n");
#endif
#if TEST_UINT128
  __int128    a;
  __uint128_t b;
  printf ("#define MULTI_WORD_INT_SIZEOF_ITEM 8\n");
#endif


  return 0;
}
