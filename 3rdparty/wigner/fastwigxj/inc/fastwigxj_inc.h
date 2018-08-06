
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

#ifndef __FASTWIGXJ_INC_H__
#define __FASTWIGXJ_INC_H__

/* This file includes functions that for performance reasons are
 * inline, and therefore have to be compiled into the using program,
 * and not loaded from a library.
 *
 * It is included by fastwigxj.h.
 *
 * Note that it also defines some additional interfaces that should be
 * considered internal.  I.e., using programs should only use
 * functions and structure layouts defined in fastwigxj.h.
 */

#include "fastwigxj_canonicalise.h"
#include "fastwigxj_hash_fcn.h"

#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

struct fastwigxj_entry
{
  uint64_t _key;
  double   _value;
};

struct wigner369j_dyn_table_mutex;

struct wigner369j_dyn_table
{
  struct fastwigxj_entry *_table;
  uint8_t                *_last_use_table;

  size_t  _table_entries;
  size_t  _table_mask;

  volatile int     _table_clear_count; /* Clear counter. */
  volatile uint8_t _table_recent;
  uint8_t          _table_thin_index;

  struct wigner369j_dyn_table_mutex *_table_mutex;

  size_t  _table_used;
  size_t  _table_used_next_check;

  uint64_t _key_mask;
};

#define WIGNER369J_DYN_TABLE_INIT \
  { NULL, NULL,  0, 0,  0, 0, 0,  NULL,  0, 0, 0 }

struct wigner369j_table
{
  struct fastwigxj_entry *_table;
  
  size_t  _table_entries;
  size_t  _table_mask;

  struct wigner369j_dyn_table *_dyn_table;

  size_t  _map_size;
  int     _fd;
};

#define WIGNER369J_TABLE_INIT \
  { NULL,  0, 0,  NULL,  0, -1 }

struct wigner369j_stats
{
  uint64_t  _hits;
  uint64_t  _dyn_hits;
  uint64_t  _dyn_trip; /* Do not include in sum of lookups. */
  uint64_t  _calc;
  uint64_t  _trivial0; /* Do not include in sum of lookups. */
  uint64_t  _9j_by_6j; /* Do not include in sum of lookups. */
  uint64_t  _dyn_table_reduce;
};

#define WIGNER369J_STATS_INIT \
  { 0, 0, 0, 0, 0, 0, 0 }

union fastwigxj_double_uint64_t_type_pun
{
  double   _f;
  uint64_t _i;
};

#ifdef __cplusplus
extern "C" {
#endif

#define FASTWIGXJPF_TABLE_3J           0
#define FASTWIGXJPF_TABLE_6J           1
#define FASTWIGXJPF_TABLE_9J           2
#define FASTWIGXJPF_TABLE_6J_FLOAT128  3
#define FASTWIGXJPF_NUM_TABLES         4

extern struct wigner369j_table fastwigxj_tables[FASTWIGXJPF_NUM_TABLES];

#define TABLE_3J (&fastwigxj_tables[FASTWIGXJPF_TABLE_3J])
#define TABLE_6J (&fastwigxj_tables[FASTWIGXJPF_TABLE_6J])
#define TABLE_9J (&fastwigxj_tables[FASTWIGXJPF_TABLE_9J])

#define TABLE_6J_128 (&fastwigxj_tables[FASTWIGXJPF_TABLE_6J_FLOAT128])

extern struct wigner369j_stats fastwigxj_stats[FASTWIGXJPF_NUM_TABLES];

#define STATS_3J (&fastwigxj_stats[FASTWIGXJPF_TABLE_3J])
#define STATS_6J (&fastwigxj_stats[FASTWIGXJPF_TABLE_6J])
#define STATS_9J (&fastwigxj_stats[FASTWIGXJPF_TABLE_9J])

#define STATS_6J_128 (&fastwigxj_stats[FASTWIGXJPF_TABLE_6J_FLOAT128])

void wig369j_ht_clear(struct wigner369j_table *table);
void wig369j_ht_deinit(struct wigner369j_table *table);
size_t wig369j_ht_init(struct wigner369j_table *table,
		       const char *filename, int type, int c14n,
		       struct fastwigxj_header *header);

void wig369j_ht_dyn_clear(struct wigner369j_dyn_table *table);
void wig369j_ht_dyn_deinit(struct wigner369j_dyn_table *table);
size_t wig369j_ht_dyn_init(struct wigner369j_dyn_table *table,
			   size_t entries, uint64_t key_mask);

int wig369j_ht_dyn_lookup(struct wigner369j_dyn_table *table,
			  struct wigner369j_stats *stats,
			  uint64_t key, double *rvalue);
int wig369j_ht_dyn_lookup_mask1(struct wigner369j_dyn_table *table,
				struct wigner369j_stats *stats,
				uint64_t sign, uint64_t key, double *rvalue);
void wig369j_ht_dyn_insert(struct wigner369j_dyn_table *dyn_table,
			   struct wigner369j_stats *stats,
			   uint64_t sign, uint64_t key, double value); 

double fastwig3jj_fallback(uint64_t x, const int *two_jv);
double fastwig6jj_fallback(uint64_t x, const int *two_jv);
double fastwig6jj_fallback_stride4(const int *two_jv);

/*double fastwig9jj_calc_by_6j(const int *two_jv);*/
double fastwig9jj_fallback(uint64_t sign, uint64_t key, const int *two_jv);

static inline void fw3jj_canon(const int *two_jv, uint64_t *rx)
{
  uint64_t x =
    wigner3j_regge_canonicalise(two_jv);
  *rx = x;
}

static inline void fw3jj_prefetch(uint64_t x)
{
  uint64_t index = x >> 1;

  if (index < TABLE_3J->_table_entries)
    {
      __builtin_prefetch(&((double *) TABLE_3J->_table)[index], 0, 0);
    }
}

static inline double fw3jj_get(const int *two_jv, uint64_t x)
{
  uint64_t sign = x & 1;
  uint64_t index = x >> 1;

  if (index < TABLE_3J->_table_entries)
    {
      double value = ((double *) TABLE_3J->_table)[index];
      STATS_3J->_hits++;
      return sign ? -value : value;
    }
  if (TABLE_3J->_dyn_table)
    {
      double value;

      if (wig369j_ht_dyn_lookup(TABLE_3J->_dyn_table, STATS_3J, index, &value))
	return sign ? -value : value;
    }
  STATS_3J->_calc++;
  return fastwig3jj_fallback(x, two_jv);
}

static inline double fw3jjl(const int *two_jv)
{
  uint64_t x;

  fw3jj_canon(two_jv, &x);
  /* Prefetch only useful when unrelated code can be placed after it,
   * before the lookup_do.
   */
  /* fw3jj_fetch(x); */
  return fw3jj_get(two_jv, x);
}
  
static inline double fw3jja(int two_j1, int two_j2, int two_j3,
			    int two_m1, int two_m2)
{
  wigner3j_symbol js;

  js.j.two_j1 = two_j1;
  js.j.two_j2 = two_j2;
  js.j.two_j3 = two_j3;
  js.j.two_m1 = two_m1;
  js.j.two_m2 = two_m2;

  return fw3jjl(js.v);
}

static inline double fw3jja6(int two_j1, int two_j2, int two_j3,
			     int two_m1, int two_m2, int two_m3)
{
  wigner3j_symbol js;

  js.j.two_j1 = two_j1;
  js.j.two_j2 = two_j2;
  js.j.two_j3 = two_j3;
  js.j.two_m1 = two_m1;
  js.j.two_m2 = two_m2;

  if (two_m1 + two_m2 + two_m3 != 0)
    return 0.0;

  return fw3jjl(js.v);
}

static inline void fw6jj_canon(const int *two_jv, uint64_t *rx)
{
  uint64_t x =
    wigner6j_regge_canonicalise_index(two_jv);
  *rx = x;
}

static inline void fw6jjs4_canon_nz(const wigner6j_4_symbol *js, uint64_t *rx)
{
  union wigner6j_type_pun_v4di_uint64_t x;

  wigner6j_struct_array4_regge_canonicalise_index_no_t0_chk(&x._v4di, js);


  int i;

  for (i = 0; i < 4; i++)
    rx[i] = x._u64[i];
}

static inline void fw6jj_prefetch(uint64_t x)
{
  if (x < TABLE_6J->_table_entries)
    {
      __builtin_prefetch(&((double *) TABLE_6J->_table)[x], 0, 0);
    }
}

static inline double fw6jj_get(const int *two_jv, uint64_t x)
{
  if (x < TABLE_6J->_table_entries)
    {
      double value = ((double *) TABLE_6J->_table)[x];
      STATS_6J->_hits++;
      return value;
    }
  if (TABLE_6J->_dyn_table)
    {
      double value;

      if (wig369j_ht_dyn_lookup(TABLE_6J->_dyn_table, STATS_6J, x, &value))
	return value;
    }
  STATS_6J->_calc++;
  return fastwig6jj_fallback(x, two_jv);
}

static inline double fw6jjs4_get(const int *two_jv, uint64_t x)
{
  if (x < TABLE_6J->_table_entries)
    {
      double value = ((double *) TABLE_6J->_table)[x];
      /*
      printf ("%016" PRIu64 " -> %.10f\n",
	      x, value);
      */
      return value;
    }
  return fastwig6jj_fallback_stride4(two_jv);
}

static inline double fw6jjl(const int *two_jv)
{
  uint64_t x;

  fw6jj_canon(two_jv, &x);
  /* Prefetch only useful when unrelated code can be placed after it,
   * before the lookup_do.
   */
  /* fw6jj_fetch(x); */
  return fw6jj_get(two_jv, x);
}

static inline double fw6jja(int two_j1, int two_j2, int two_j3,
			    int two_j4, int two_j5, int two_j6)
{
  wigner6j_symbol js;

  js.j.two_j1 = two_j1;
  js.j.two_j2 = two_j2;
  js.j.two_j3 = two_j3;
  js.j.two_j4 = two_j4;
  js.j.two_j5 = two_j5;
  js.j.two_j6 = two_j6;

  return fw6jjl(js.v);
}

static inline void fw9jj_canon(const int *two_jv, uint64_t *rkey, uint64_t *rx)
{
  uint64_t key =
    wigner9j_canonicalise(two_jv);

  // and do the search in the hash-table

  uint64_t x = key & ~((uint64_t) 1);

  x = WIGNER369_HASH_FCN(x);

  x &= TABLE_9J->_table_mask;

  // We now have the address!

  *rkey = key;
  *rx = x;

}

static inline void fw9jj_prefetch(uint64_t x)
{
  __builtin_prefetch(&TABLE_9J->_table[x]._key, 0, 0);
}

static inline double fw9jj_get(const int *two_jv, uint64_t key, uint64_t x)
{
  uint64_t sign = key & 1;
  key &= ~((uint64_t) 1);

  for ( ; ; )
    {
      uint64_t table_key = TABLE_9J->_table[x]._key;

      if ((table_key & ~((uint64_t) 1)) == key)
	{
	  double value = TABLE_9J->_table[x]._value;
	  STATS_9J->_hits++;
	  return (sign & table_key) ? -value : value;
	}

      if (table_key == (uint64_t) -1)
	{
	  if (TABLE_9J->_dyn_table)
	    {
	      double value;

	      if (wig369j_ht_dyn_lookup_mask1(TABLE_9J->_dyn_table, STATS_9J,
					      sign, key, &value))
		return value;
	    }
	  STATS_9J->_calc++;
	  return fastwig9jj_fallback(sign, key, two_jv);
	}

      x++;
      x &= TABLE_9J->_table_mask;
    }
}

static inline double fw9jjl(const int *two_jv)
{
  uint64_t key, x;

  fw9jj_canon(two_jv, &key, &x);
  /* Prefetch only useful when unrelated code can be placed after it,
   * before the lookup_do.
   */
  /* fw9jj_fetch(x); */
  return fw9jj_get(two_jv, key, x);
}

static inline double fw9jja(int two_j1, int two_j2, int two_j3,
			    int two_j4, int two_j5, int two_j6,
			    int two_j7, int two_j8, int two_j9)
{
  wigner9j_symbol js;

  js.j.two_j1 = two_j1;
  js.j.two_j2 = two_j2;
  js.j.two_j3 = two_j3;
  js.j.two_j4 = two_j4;
  js.j.two_j5 = two_j5;
  js.j.two_j6 = two_j6;
  js.j.two_j7 = two_j7;
  js.j.two_j8 = two_j8;
  js.j.two_j9 = two_j9;

  return fw9jjl(js.v);
}

#ifdef __cplusplus
}
#endif

#endif/*__FASTWIGXJ_INC_H__*/
