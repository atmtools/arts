
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

#ifndef __FASTWIGXJ_STRUCT_H__
#define __FASTWIGXJ_STRUCT_H__

#include "fastwigxj_vector.h"

/* When putting the 9j together, it makes more sense to put the
 * elements in a vector (struct) in memory, than to keep them in
 * registers.  This also fits the canonicaliser much better, as it
 * then can do just a few large loads from memory.
 */

typedef union wigner9j_symbol_t
{
  int v[9] __attribute__ ((aligned (32)));

  struct
  {
    int two_j1, two_j2, two_j3;
    int two_j4, two_j5, two_j6;
    int two_j7, two_j8, two_j9;
    int _dummy[7];
  } j;

} wigner9j_symbol;

typedef union wigner6j_symbol_t
{
  int v[8] __attribute__ ((aligned (32)));

  struct
  {
    int two_j1, two_j2, two_j3;
    int two_j4, two_j5, two_j6;
    int _dummy[2];
  } j;

} wigner6j_symbol;

typedef union wigner6j_4_symbol_t
{
  int v[4*6] __attribute__ ((aligned (32)));

  struct
  {
    int two_j1[4], two_j2[4], two_j3[4];
    int two_j4[4], two_j5[4], two_j6[4];
  } j;

  struct
  {
    v4su32 two_j1, two_j2, two_j3;
    v4su32 two_j4, two_j5, two_j6;
  } jv;

} wigner6j_4_symbol;

typedef union wigner3j_symbol_t
{
  int v[8] __attribute__ ((aligned (32)));

  struct
  {
    int two_j1, two_j2, two_j3;
    int two_m1, two_m2;
    int two_m3; /* not used */
    int _dummy[2];
  } j;

} wigner3j_symbol;

#endif/*__FASTWIGXJ_STRUCT_H__*/
