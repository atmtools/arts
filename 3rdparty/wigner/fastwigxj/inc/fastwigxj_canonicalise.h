
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

#ifndef __CANONICALISE_H__
#define __CANONICALISE_H__

/* All interfaces in this header are internal.
 * Do not use directly.
 */

#include "fastwigxj_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

union wigner6j_type_pun_v4di_uint64_t
{
  v4di     _v4di;
  uint64_t _u64[4];
};

uint64_t wigner3j_regge_canonicalise(const int *two_jv);

uint64_t wigner6j_regge_canonicalise_index(const int *two_jv);

void wigner6j_struct_array4_regge_canonicalise_index_no_t0_chk(v4su64 *key, const wigner6j_4_symbol *js);

uint64_t wigner9j_canonicalise(const int *two_jv);

#ifdef __cplusplus
}
#endif

#endif/*__CANONICALISE_H__*/
