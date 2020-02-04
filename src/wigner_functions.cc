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
#include <algorithm>
#include "arts_omp.h"

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

  wig_temp_init(j);
  g = WIGNER3(a, b, c, d, e, f);
  wig_temp_free();

  return Numeric(g);
}

Numeric wigner6j(const Rational j1,
                 const Rational j2,
                 const Rational j3,
                 const Rational l1,
                 const Rational l2,
                 const Rational l3) {
  const int a = (2 * j1).toInt(), b = (2 * j2).toInt(), c = (2 * j3).toInt(),
            d = (2 * l1).toInt(), e = (2 * l2).toInt(), f = (2 * l3).toInt();
  double g;
  const int j = std::max({std::abs(a),
                          std::abs(b),
                          std::abs(c),
                          std::abs(d),
                          std::abs(e),
                          std::abs(f)});

  wig_temp_init(j);
  g = WIGNER6(a, b, c, d, e, f);
  wig_temp_free();

  return Numeric(g);
}

Numeric co2_ecs_wigner_symbol(
    int Ji, int Jf, int Ji_p, int Jf_p, int L, int li, int lf) {
  return WIGNER3(Ji_p, L, Ji, li, 0, -li) * WIGNER3(Jf_p, L, Jf, -lf, 0, lf) *
         WIGNER6(Ji, Jf, 2, Jf_p, Ji_p, L) * Numeric(L + 1);
}

Numeric o2_ecs_wigner_symbol(
    int Nl, int Nk, int Jl, int Jk, int Jl_p, int Jk_p, int L) {
  return WIGNER3(Nl, Nk, L, 0, 0, 0) * WIGNER6(L, Jk, Jl, 2, Nl, Nk) *
         WIGNER6(L, Jk_p, Jl_p, 2, Nl, Nk) * WIGNER6(L, Jk, Jl, 2, Jl_p, Jk_p);
}

Numeric o2_ecs_wigner_symbol_tran(
  int Ni,  int Nf,  int Ji,  int Jf,
  int Nip, int Nfp, int Jip, int Jfp,
  int L, int Si, int Sf, int n)
{
  auto a = WIGNER3(2*Nip, 2*Ni, 2*L, 0, 0, 0);
  auto b = WIGNER3(2*Nfp, 2*Nf, 2*L, 0, 0, 0);
  auto c = WIGNER6(2*L, 2*Ji, 2*Jip, 2*Si, 2*Nip, 2*Ni);
  auto d = WIGNER6(2*L, 2*Jf, 2*Jfp, 2*Sf, 2*Nfp, 2*Nf);
  auto e = WIGNER6(2*L, 2*Ji, 2*Jip, 2*n, 2*Jfp, 2*Jf);
  return (a * b * c * d * e) * (2*L + 1);
}

bool is_wigner_ready(int j) {
  extern int wigxjpf_max_prime_decomp;
  return not(j > wigxjpf_max_prime_decomp);
}

bool is_wigner3_ready(const Rational& J) {
  const int test = (J * 6).toInt() / 2 + 1;  // nb. J can be half-valued
  return is_wigner_ready(test);
}

bool is_wigner6_ready(const Rational& J) {
  const int test = (J * 4).toInt() + 1;  // nb. J can be half-valued
  return is_wigner_ready(test);
}
