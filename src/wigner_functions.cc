/* Copyright (C) 2012
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/**
 * @file wigner_functions.cc
 * @author Richard Larsson
 * @date 2013-06-19
 * 
 * @brief Wigner symbol interactions
 */

#include "wigner_functions.h"

#include <sys/errno.h>

#include <algorithm>

#include "arts_omp.h"
#include "constants.h"
#include "debug.h"

#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

Numeric wigner3j(const Rational j1,
                 const Rational j2,
                 const Rational j3,
                 const Rational m1,
                 const Rational m2,
                 const Rational m3) {
  errno = 0;

  const int a = (2 * j1).toInt(), b = (2 * j2).toInt(), c = (2 * j3).toInt(),
            d = (2 * m1).toInt(), e = (2 * m2).toInt(), f = (2 * m3).toInt();
  double g;
  const int j = std::max({std::abs(a),
                          std::abs(b),
                          std::abs(c),
                          std::abs(d),
                          std::abs(e),
                          std::abs(f)}) *
                    3 / 2 +
                1;

  wig_thread_temp_init(j);
  g = WIGNER3(a, b, c, d, e, f);
  wig_temp_free();

  if (errno == EDOM) {
    errno = 0;
    ARTS_USER_ERROR("Bad state, perhaps you need to call Wigner3Init?")
  }

  return Numeric(g);
}

Numeric wigner6j(const Rational j1,
                 const Rational j2,
                 const Rational j3,
                 const Rational l1,
                 const Rational l2,
                 const Rational l3) {
  errno = 0;

  const int a = (2 * j1).toInt(), b = (2 * j2).toInt(), c = (2 * j3).toInt(),
            d = (2 * l1).toInt(), e = (2 * l2).toInt(), f = (2 * l3).toInt();
  double g;
  const int j = std::max({std::abs(a),
                          std::abs(b),
                          std::abs(c),
                          std::abs(d),
                          std::abs(e),
                          std::abs(f)});

  wig_thread_temp_init(j);
  g = WIGNER6(a, b, c, d, e, f);
  wig_temp_free();

  if (errno == EDOM) {
    errno = 0;
    ARTS_USER_ERROR("Bad state, perhaps you need to call Wigner6Init?")
  }

  return Numeric(g);
}

std::pair<Rational, Rational> wigner_limits(std::pair<Rational, Rational> a,
                                            std::pair<Rational, Rational> b) {
  const bool invalid = a.first.isUndefined() or b.first.isUndefined();
  if (invalid) {
    return {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED};
  }

  std::pair<Rational, Rational> out{max(a.first, b.first),
                                    min(a.second, b.second)};
  if (out.first > out.second) out = {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED};
  return out;
}

Index make_wigner_ready(int largest, [[maybe_unused]] int fastest, int size) {
  if (size == 3) {
#if DO_FAST_WIGNER
    fastwigxj_load(FAST_WIGNER_PATH_3J, 3, NULL);
#ifdef _OPENMP
    fastwigxj_thread_dyn_init(3, fastest);
#else
    fastwigxj_dyn_init(3, fastest);
#endif
#endif
    wig_table_init(largest, 3);

    return largest;
  }

  if (size == 6) {
#if DO_FAST_WIGNER
    fastwigxj_load(FAST_WIGNER_PATH_3J, 3, NULL);
    fastwigxj_load(FAST_WIGNER_PATH_6J, 6, NULL);
#ifdef _OPENMP
    fastwigxj_thread_dyn_init(3, fastest);
    fastwigxj_thread_dyn_init(6, fastest);
#else
    fastwigxj_dyn_init(3, fastest);
    fastwigxj_dyn_init(6, fastest);
#endif
#endif
    wig_table_init(largest * 2, 6);

    return largest;
  }

  return 0;
}

bool is_wigner_ready(int j) {
  extern int wigxjpf_max_prime_decomp;
  return not(j > wigxjpf_max_prime_decomp);
}

bool is_wigner3_ready(const Rational& J) {
  const int test = J.toInt(6) / 2 + 1;  // nb. J can be half-valued
  return is_wigner_ready(test);
}

bool is_wigner6_ready(const Rational& J) {
  const int test = J.toInt(4) + 1;  // nb. J can be half-valued
  return is_wigner_ready(test);
}
