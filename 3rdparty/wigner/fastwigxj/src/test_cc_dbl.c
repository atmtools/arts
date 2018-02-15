
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

/* This small program is for compile-time test of availability of long
 * double, threading support,
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
#if TEST_AVX2
#include <x86intrin.h>
typedef long long v4di __attribute__ ((vector_size (32)));
#endif
#if TEST_SSE4_1
#include <x86intrin.h>
typedef double v2df __attribute__ ((vector_size (16)));
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
  printf ("#define FASTWIGXJ_USE_LONG_DOUBLE 1\n");
#endif  
#if TEST_FLOAT128
  __float128 a;
  __float128 b = 1.14, c = 2.00159;
  char s[64];

  a = b + c;

  quadmath_snprintf(s, sizeof(s), "%Qf", a);

  printf ("%s\n", s);
  printf ("#define FASTWIGXJ_USE_FLOAT128 1\n");
#endif
#if TEST_THREAD
  global = 2;
  printf ("#define FASTWIGXJ_HAVE_THREAD 1\n");
#endif
#if TEST_AVX2
  // It makes no sense to try to compile if the compiler
  // does not support.  And there was a machine crash on seeing
  // instructions it could not perform...
#if defined(__AVX2__)
  v4di a = { 1, 0., 0., 0. }, b = { 2, 0., 0., 0. };
  v4di c;
  c = _mm256_add_epi64(a,b);
  printf ("%lld\n",c[0]);
  printf ("#define FASTWIGXJ_HAVE_AVX2 1\n");
#endif
#endif
#if TEST_SSE4_1
#if defined(__SSE4_1__)
  v2df x = { 1.2, 0.}, y = { 1.3, 0.}, c = { 0., 0.};
  v2df z;
  z = __builtin_ia32_blendvpd(y,x,c);
  printf ("%f\n",z[0]);
  printf ("#define FASTWIGXJ_HAVE_SSE4_1 1\n");
#endif
#endif

  return 0;
}
