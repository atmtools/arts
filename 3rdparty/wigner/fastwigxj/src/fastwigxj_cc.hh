
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

/* TO BE REMOVED FROM TAR-BALL. */

#ifndef __FASTWIGXJ_HH__
#define __FASTWIGXJ_HH__

#include "fastwigxj.h"
#include "fastwigxj_config.h"
#if FASTWIGXJ_USE_FLOAT128
#include "fastwigxj_quadmath_inc.h"
#endif

#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

class wigner369j_hash
{
public:
  wigner369j_hash();
  ~wigner369j_hash();

public:
  void init(const char *filename, int type,
	    fastwigxj_header *header = NULL);
  void deinit();
};

class wigner3j_table :
  public wigner369j_hash
{
public:
  void init(const char *filename, fastwigxj_header *header = NULL)
  {
    wigner369j_hash::init(filename, 3, header);
  }

public:
  double lookup(const int *two_jv) const;

  double lookup(int two_ja, int two_jb, int two_jc,
		int two_ma, int two_mb) const;
};

inline double wigner3j_table::lookup(const int *two_jv) const
{
  return fw3jjl(two_jv);
}

inline double wigner3j_table::lookup(int two_ja, int two_jb, int two_jc,
				     int two_ma, int two_mb) const
{
  return fw3jja(two_ja, two_jb, two_jc, 
		two_ma, two_mb);
}

class wigner6j_table :
  public wigner369j_hash
{
public:
  void init(const char *filename, fastwigxj_header *header = NULL)
  {
    wigner369j_hash::init(filename, 6, header);
  }

public:
  double lookup(const int *two_jv) const
  {
    return fw6jjl(two_jv);
  }

  double lookup(int two_ja, int two_jb, int two_jc,
		int two_jd, int two_je, int two_jf) const
  {
    int two_jv[6] = { two_ja, two_jb, two_jc,
			   two_jd, two_je, two_jf};

    return fw6jjl(two_jv);
  }

  void lookup_pre(const int *two_jv,
		  uint64_t &rkey_dummy, uint64_t &rx) const
  {
    (void) rkey_dummy;
    fw6jj_canon(two_jv, &rx);
  }

  void lookup_prefetch(uint64_t x) const
  {
    fw6jj_prefetch(x);
  }

  double lookup_do(const int *two_jv,
		   uint64_t key_dummy, uint64_t x) const
  {
    (void) key_dummy;
    return fw6jj_get(two_jv, x);
  }
};

#if FASTWIGXJ_USE_FLOAT128
class wigner6j_table_float128 :
  public wigner369j_hash
{
public:
  void init(const char *filename, fastwigxj_header *header = NULL)
  {
    wigner369j_hash::init(filename, 7, header);
  }

public:
  double lookup(const int *two_jv) const
  {
    return (double) fw6jjl_float128(two_jv);
  }

  double lookup(int two_ja, int two_jb, int two_jc,
		int two_jd, int two_je, int two_jf) const
  {
    int two_jv[6] = { two_ja, two_jb, two_jc,
		      two_jd, two_je, two_jf};

    return (double) fw6jjl_float128(two_jv);
  }
};
#endif

class wigner9j_hash :
  public wigner369j_hash
{
public:
  void init(const char *filename, fastwigxj_header *header = NULL)
  {
    wigner369j_hash::init(filename, 9, header);
  }

public:
  void set_table_6j(wigner6j_table *table_6j);

public:
  double lookup(const int *two_jv) const
  {
    return fw9jjl(two_jv);
  }

  void lookup_pre(const int *two_jv,
		  uint64_t &rkey, uint64_t &rx) const
  {
    fw9jj_canon(two_jv, &rkey, &rx);
  }

  void lookup_prefetch(uint64_t x) const
  {
    fw9jj_prefetch(x);
  }

  double lookup_do(const int *two_jv,
		   uint64_t key, uint64_t x) const
  {
    return fw9jj_get(two_jv, key, x);
  }
};

#endif//__FASTWIGXJ_HH__
