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

Numeric dwigner3j(Index M, Index J1, Index J2, Index J) {
  auto CJM = [](Index J, Index M) {
    Numeric CJM = 1.;
    for (Index I = 1; I <= J; I++) CJM *= (1. - .5 / static_cast<Numeric>(I));
    for (Index K = 1; K <= M; K++)
      CJM *= static_cast<Numeric>(J + 1 - K) / static_cast<Numeric>(J + K);
    return CJM;
  };

  Numeric GCM = 0.;
  if (J1 < 0 or J2 < 0 or J < 0) return GCM;
  Index JM = J1 + J2;
  Index J0 = std::abs(J1 - J2);
  if (J > JM or J < J0) return GCM;
  Index JS = std::max(J1, J2);
  Index JI = std::min(J1, J2);
  Index MA = std::abs(M);
  if (MA > JI) return GCM;
  Index UN = 1 - 2 * (JS % 2);
  Index QM = M + M;
  Numeric CG0 = 0.;
  GCM = static_cast<Numeric>(UN) *
        std::sqrt(CJM(JI, MA) / CJM(JS, MA) * CJM(J0, 0) /
                  static_cast<Numeric>(JS + JS + 1));
  Index AJ0 = J0;
  Index AJM = JM + 1;
  Index AJ02 = AJ0 * AJ0;
  Index AJM2 = AJM * AJM;
  Numeric ACG0 = 0.;
  for (Index I = J0 + 1; I <= J; I++) {
    Index AI = I;
    Index AI2 = AI * AI;
    Numeric ACG = std::sqrt((AJM2 - AI2) * (AI2 - AJ02));
    Numeric CG1 =
        (static_cast<Numeric>(QM) * static_cast<Numeric>(I + I - 1) * GCM -
         ACG0 * CG0) /
        ACG;
    CG0 = GCM;
    GCM = CG1;
    ACG0 = ACG;
  }
  return GCM;
}

Numeric dwigner6j(Index A, Index B, Index C, Index D, Index F) {
  Numeric SIXJ, TERM;
  if (std::abs(A - C) > F or std::abs(B - D) > F or (A + C) < F or (B + D < F))
    goto x1000;
  switch (C - D + 2) {
    case 2:
      goto x2;
    case 3:
      goto x3;
    case 4:
      goto x1000;
    default:
      goto x1;
  }

x1:
  switch (A - B + 2) {
    case 2:
      goto x11;
    case 3:
      goto x12;
    default:
      goto x10;
  }

x10:
  TERM =
      (static_cast<Numeric>(F + B + D + 1) * static_cast<Numeric>(F + B + D) *
       static_cast<Numeric>(B + D - F) * static_cast<Numeric>(B + D - (1 + F)));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B) *
           static_cast<Numeric>(2 * B - 1) * static_cast<Numeric>(D) *
           static_cast<Numeric>(2 * D - 1) * static_cast<Numeric>(2 * D + 1));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C + F)) * std::sqrt(TERM);
  return SIXJ;

x11:
  TERM = (static_cast<Numeric>(F + B + D + 1) *
          static_cast<Numeric>(F + B - D + 1) *
          static_cast<Numeric>(F + D - B) * static_cast<Numeric>(B + D - F));
  TERM /= (4. * static_cast<Numeric>(B) * static_cast<Numeric>(2 * B + 1) *
           static_cast<Numeric>(B + 1) * static_cast<Numeric>(D) *
           static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(2 * D - 1));
  SIXJ =
      static_cast<Numeric>(pow_negative_one(A + C - F - 1)) * std::sqrt(TERM);
  return SIXJ;

x12:
  TERM =
      (static_cast<Numeric>(F + D - B) * static_cast<Numeric>(F + D - B - 1) *
       static_cast<Numeric>(F + B - D + 2) *
       static_cast<Numeric>(F + B - D + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B + 1) *
           static_cast<Numeric>(2 * B + 3) * static_cast<Numeric>(2 * D - 1) *
           static_cast<Numeric>(D) * static_cast<Numeric>(2 * D + 1));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x2:
  switch (A - B + 2) {
    case 2:
      goto x21;
    case 3:
      goto x22;
    default:
      goto x20;
  }

x20:
  TERM =
      (static_cast<Numeric>(F + B + D + 1) * static_cast<Numeric>(D + B - F) *
       static_cast<Numeric>(F + B - D) * static_cast<Numeric>(F + D - B + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B) *
           static_cast<Numeric>(2 * B - 1) * static_cast<Numeric>(D) *
           static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(D + 1));
  SIXJ =
      static_cast<Numeric>(pow_negative_one(A + C - F - 1)) * std::sqrt(TERM);
  return SIXJ;

x21:
  TERM = (static_cast<Numeric>(B) * static_cast<Numeric>(B + 1) +
          static_cast<Numeric>(D) * static_cast<Numeric>(D + 1) -
          static_cast<Numeric>(F) * static_cast<Numeric>(F + 1));
  TERM /=
      std::sqrt(4. * static_cast<Numeric>(B) * static_cast<Numeric>(B + 1) *
                static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(D) *
                static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(D + 1));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C - F - 1)) * TERM;
  return SIXJ;

x22:
  TERM =
      (static_cast<Numeric>(F + D + B + 2) *
       static_cast<Numeric>(F + B - D + 1) *
       static_cast<Numeric>(B + D - F + 1) * static_cast<Numeric>(F + D - B));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B + 1) *
           static_cast<Numeric>(2 * B + 3) * static_cast<Numeric>(D) *
           static_cast<Numeric>(D + 1) * static_cast<Numeric>(2 * D + 1));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x3:
  switch (A - B + 2) {
    case 2:
      goto x31;
    case 3:
      goto x32;
    default:
      goto x30;
  }

x30:
  TERM =
      (static_cast<Numeric>(F + B - D) * static_cast<Numeric>(F + B - D - 1) *
       static_cast<Numeric>(F + D - B + 2) *
       static_cast<Numeric>(F + D - B + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(B) *
           static_cast<Numeric>(2 * B - 1) * static_cast<Numeric>(D + 1) *
           static_cast<Numeric>(2 * D + 1) * static_cast<Numeric>(2 * D + 3));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x31:
  TERM =
      (static_cast<Numeric>(F + D + B + 2) *
       static_cast<Numeric>(B + D - F + 1) *
       static_cast<Numeric>(F + D - B + 1) * static_cast<Numeric>(F + B - D));
  TERM /= (4. * static_cast<Numeric>(B) * static_cast<Numeric>(2 * B + 1) *
           static_cast<Numeric>(B + 1) * static_cast<Numeric>(2 * D + 1) *
           static_cast<Numeric>(D + 1) * static_cast<Numeric>(2 * D + 3));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x32:
  TERM = (static_cast<Numeric>(F + D + B + 3) *
          static_cast<Numeric>(F + B + D + 2) *
          static_cast<Numeric>(B + D - F + 2) *
          static_cast<Numeric>(B + D - F + 1));
  TERM /= (4. * static_cast<Numeric>(2 * B + 3) * static_cast<Numeric>(B + 1) *
           static_cast<Numeric>(2 * B + 1) * static_cast<Numeric>(2 * D + 3) *
           static_cast<Numeric>(D + 1) * static_cast<Numeric>(2 * D + 1));
  SIXJ = static_cast<Numeric>(pow_negative_one(A + C - F)) * std::sqrt(TERM);
  return SIXJ;

x1000:
  SIXJ = 0.;
  return SIXJ;
}
