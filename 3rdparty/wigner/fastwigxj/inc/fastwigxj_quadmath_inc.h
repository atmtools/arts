
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

#ifndef __FASTWIGXJ_QUADMATH_INC_H__
#define __FASTWIGXJ_QUADMATH_INC_H__

#ifdef __cplusplus
extern "C" {
#endif

void fastwig6jj_float128_fallback(__float128 *result, const int *two_jv);
void fastwig6jj_float128_fallback_stride4(__float128 *result,const int *two_jv);

// #if FASTWIGXJ_USE_FLOAT128
static inline void fw6jj_prefetch_float128(uint64_t x)
{
  if (x < TABLE_6J_128->_table_entries)
    {
      __builtin_prefetch(&((__float128 *) TABLE_6J_128->_table)[x], 0, 0);
    }
}
// #endif

// #if FASTWIGXJ_USE_FLOAT128
static inline void fw6jjs4_get_float128(__float128 *result,
					const int *two_jv, uint64_t x)
{
  if (x < TABLE_6J_128->_table_entries)
    {
      __float128 value = ((__float128 *) TABLE_6J_128->_table)[x];
      /*
      printf ("%016" PRIu64 " -> %.10f\n",
	      x, value);
      */
      *result = value;
      STATS_6J_128->_hits++;
      return;
    }
  STATS_6J_128->_calc++;
  fastwig6jj_float128_fallback_stride4(result, two_jv);
}
// #endif

static inline void fw6jj_get_float128(__float128 *result,
				      const int *two_jv, uint64_t x)
{
  if (x < TABLE_6J_128->_table_entries)
    {
      *result = ((__float128 *) TABLE_6J_128->_table)[x];
      STATS_6J_128->_hits++;
      return;
    }
  STATS_6J_128->_calc++;
  fastwig6jj_float128_fallback(result, two_jv);
}

static inline __float128 fw6jjl_float128(const int *two_jv)
{
  uint64_t x;

  fw6jj_canon(two_jv, &x);
  /* Prefetch only useful when unrelated code can be placed after it,
   * before the lookup_do.
   */
  /* fw6jj_fetch(x); */
  __float128 value;
  fw6jj_get_float128(&value, two_jv, x);
  return value;
}

#ifdef __cplusplus
}
#endif

#endif/*__FASTWIGXJ_QUADMATH_INC_H__*/
