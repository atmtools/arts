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

#include "constants.h"
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


struct Wigner3JTriangleLimit {
  int lower;
  int upper;
};



/** Finds the upper and lower limits of L3 given a Wigner-3J symbol:
 * 
 * /            \
 * |  L1 L2 L3  |
 * |            |
 * |  M1 M2 M3  |
 * \            /
 * 
 * Based on abs(L1-L2) <= L3 <= L1+L2, where L1 and L2 are positives
 * 
 * @param[in] L1 As in Wigner symbol
 * @param[in] L2 As in Wigner symbol
 * @return The limits
 */
constexpr Wigner3JTriangleLimit find_wigner3j_limits(const Rational L1,
                                                     const Rational L2)
{
  return {abs(L1-L2).toInt(), (L1+L2).toInt()};
}


/** Combines two limits and return the highest even-numbered low value and the
 * lowest numbered even high value
 * 
 * @param[in] lim1 Limit 1
 * @param[in] lim2 Limit 2
 * @return More limited limit
 */
constexpr Wigner3JTriangleLimit find_even_limits(const Wigner3JTriangleLimit lim1,
                                                 const Wigner3JTriangleLimit lim2)
{
  return 
  {
  lim1.lower > lim2.lower ?
    (lim1.lower % 2 ?
      lim1.lower + 1 :
      lim1.lower) :
    (lim2.lower % 2 ?
      lim2.lower + 1 :
      lim2.lower),
  lim1.upper < lim2.upper ?
    (lim1.upper % 2 ?
      lim1.upper - 1 :
      lim1.upper) :
    (lim2.upper % 2?
      lim2.upper - 1 :
      lim2.upper)
  };
}


Numeric o2_ecs_erot_jn_same(Rational J)
{
  using Constant::h;
  using Constant::pow2;
  using Constant::pow3;
  
  constexpr auto B = 43100.44276e6;
  constexpr auto D = 145.1271e3;
  constexpr auto H = 49e-3;
  constexpr auto lB = 59501.3438e6;
  constexpr auto lD = 58.3680e3;
  constexpr auto lH = 290.8e-3;
  constexpr auto gB = -252.58634e6;
  constexpr auto gD = -243.42;
  constexpr auto gH = -1.46e-3;
  
  const auto lJ = (J * (J + 1)).toNumeric();
  return h * (
  B*lJ - D*pow2(lJ) + H*pow3(lJ) -
  (gB + gD*lJ + gH*pow2(lJ)) +
  2.0/3.0 * (lB + lD*lJ + lH*pow2(lJ)));
}


Numeric o2_ecs_ql_makarov(Rational L, Numeric T)
{
  using Constant::k;
  
  auto lambda =  0.39;
  auto beta = 0.567;
  auto el = o2_ecs_erot_jn_same(L);
  
  return
  (2*L + 1).toNumeric() / pow(L * (L+1), lambda) * std::exp(-beta * el / (k*T));
}


Numeric o2_ecs_adiabatic_factor_makarov(Rational N, Numeric T)
{
  using Conversion::angstrom2meter;
  using Constant::h_bar;
  using Constant::pow2;
  using Constant::m_u;
  using Constant::pi;
  using Constant::k;
  
  auto dc =  angstrom2meter(0.61);
  
  auto en = o2_ecs_erot_jn_same(N);
  auto enm2 = o2_ecs_erot_jn_same(N-2);
  auto wnnm2 = (en - enm2) / h_bar;
  
  auto v_bar = std::sqrt(8*k*T/(pi * 31.989830 * m_u));
  auto tauc = dc / v_bar;
  
  return 1.0 / pow2(1 + 1.0/24.0 * pow2(wnnm2 * tauc));
}


Numeric o2_ecs_wigner_symbol_tran(
  const Rational& Ji,
  const Rational& Jf,
  const Rational& Ni,
  const Rational& Nf,
  const Rational& Si,
  const Rational& Sf,
  const Rational& Ji_p,
  const Rational& Jf_p,
  const Rational& Ni_p,
  const Rational& Nf_p,
  const Rational& n,
  const Numeric& T)
{
  Numeric o2_ecs_wigner_symbol_tran = 0;
  
  // Limits above 2 and below the largest index there is
  auto lims = find_even_limits({2, std::numeric_limits<int>::max()},
                               find_even_limits(find_wigner3j_limits(Ni_p, Ni),
                                                find_wigner3j_limits(Nf_p, Nf)));
  
  auto f = sqrt(2*Ni+1) * sqrt(2*Ni_p+1) * sqrt(2*Jf+1) * sqrt(2*Jf_p+1) * sqrt(2*Nf+1) * sqrt(2*Nf_p+1) * (2*Ji_p+1).toNumeric();
  auto g = even(Ji_p + Ji + n) ? 1 : -1;  // -1^(Ji_p + Ji + n)
  auto OmegaNi = o2_ecs_ql_makarov(Ni, T);
  
  for (int L=lims.lower; L<=lims.upper; L+=2) {
    auto OmegaL = o2_ecs_adiabatic_factor_makarov(L, T);
    auto QL = o2_ecs_ql_makarov(L, T);
    auto a = WIGNER3(Ni_p.toInt(2), Ni.toInt(2), 2*L, 0, 0, 0);
    auto b = WIGNER3(Nf_p.toInt(2), Nf.toInt(2), 2*L, 0, 0, 0);
    auto c = WIGNER6(2*L, Ji.toInt(2), Ji_p.toInt(2), Si.toInt(2), Ni_p.toInt(2), Ni.toInt(2));
    auto d = WIGNER6(2*L, Jf.toInt(2), Jf_p.toInt(2), Sf.toInt(2), Nf_p.toInt(2), Nf.toInt(2));
    auto e = WIGNER6(2*L, Ji.toInt(2), Ji_p.toInt(2), n.toInt(2), Jf_p.toInt(2), Jf.toInt(2));
    
    o2_ecs_wigner_symbol_tran += (a * b * c * d * e * f * g) * (2*L + 1) * QL * (OmegaNi / OmegaL);
  }
  
  return o2_ecs_wigner_symbol_tran;
}


Numeric o2_makarov2013_reduced_dipole(const Rational& Jup, const Rational& Jlo, const Rational& N)
{
  return 
  (even(Jlo + N) ? 1 : -1) *
  sqrt(6 * (2*Jlo + 1) * (2*Jup + 1)) *
  WIGNER6(2, 2, 2, Jup.toInt(2), Jlo.toInt(2), N.toInt(2)) / (2*N + 1).toNumeric();
}


Index make_wigner_ready(int largest, 
       [[maybe_unused]] int fastest,
                        int size)
{
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
  } else if (size == 6) {
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
  } else {
    return 0;
  }
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
