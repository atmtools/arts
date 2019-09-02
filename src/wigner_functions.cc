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

// #include "wigxjpf/inc/wigxjpf.h"
#include "wigner_functions.h"
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include "arts_omp.h"

#if DO_FAST_WIGNER
#define WIGNER3 fw3jja6
#define WIGNER6 fw6jja
#else
#define WIGNER3 wig3jj
#define WIGNER6 wig6jj
#endif

/*!
  Run wigxjpf wig3jj for Rational symbol
  
  /                \
  |  j1   j2   j3  |
  |                |
  |  m1   m2   m3  |
  \                /
  
  See for definition: http://dlmf.nist.gov/34.2
*/
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

/*!
 Run wigxjpf wig6jj for Rational symbol
 
 /                \
 |  j1   j2   j3  |
 <                >
 |  l1   l2   l3  |
 \                /
 
 See for definition: http://dlmf.nist.gov/34.4
*/
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

/*! Returns the wigner symbol used in Niro etal 2004

 Symbol:
 /              \  /              \  /               \
 | Ji_p  L   Ji |  |  Jf_p  L  Jf |  | Ji    Jf    1 |
 |              |  |              |  <               >  (2L + 1)
 | li    0  -li |  | -lf    0  lf |  | Jf_p  Ji_p  L |
 \              /  \              /  \               /
 
 Note: The wigner library takes two times the physical values
       so, e.g., the 1 must be 2.  This hold true for all user inputs as well!
       
 Reference: 
 Spectra calculations in central and wing regions of CO2 IR bands between 10 and 20 mcrons.
 I: model and laboratory measurements. F. Niro, C. Boulet, J.-M. Hartmann. 
 JQSRT 88 (2004) 483 – 498. Equation 4 page 488.
 
 Note: Ignore typos, the above is tested in relmat
*/
Numeric co2_ecs_wigner_symbol(
    int Ji, int Jf, int Ji_p, int Jf_p, int L, int li, int lf) {
  return WIGNER3(Ji_p, L, Ji, li, 0, -li) * WIGNER3(Jf_p, L, Jf, -lf, 0, lf) *
         WIGNER6(Ji, Jf, 2, Jf_p, Ji_p, L) * Numeric(L + 1);
}

/*! Returns the wigner symbol used in Makarov etal 2013
 
 Symbol:
 /             \  /             \  /                 \  /                 \
 |  Nl  Nk  L  |  |  L  Jk  Jl  |  |  L  Jk_p  Jl_p  |  |  L  Jk    Jl    |
 |             |  <             >  <                 >  <                 >
 |  0   0   0  |  |  1  Nl  Nk  |  |  1  Nl    Nk    |  |  1  Jl_p  Jk_p  |
 \             /  \             /  \                 /  \                 /
 
 
 Note: The wigner library takes two times the physical values
       so, e.g., the 1 must be 2.  This hold true for all user inputs as well!
  
  Reference:
  D.S. Makarov, M.Yu. Tretyakov, C. Boulet,
  Line mixing in the 60-GHz atmospheric oxygen band: Comparison of the MPM and ECS model,
  Journal of Quantitative Spectroscopy and Radiative Transfer,
  Volume 124,
  2013,
  Pages 1-10,
  ISSN 0022-4073,
  https://doi.org/10.1016/j.jqsrt.2013.02.019.
  (http://www.sciencedirect.com/science/article/pii/S0022407313000745)
  Abstract: Precise profiles of the 60-GHz molecular oxygen band recorded earlier in a wide temperature range by means of a resonator spectrometer at atmospheric pressure were treated. High signal-to-noise ratio allows careful study of the band shape taking into consideration the mixing effect. Comparative analysis of the band profile calculated by an extended MPM (Millimeter-wave Propagation Model) and by the ECS (Energy Corrected Sudden) approximation model is presented. Some limitations of the MPM approach are discussed on the basis of the comparison of the two models. Interbranch coupling is shown to make a noticeable contribution to absorption which means that it should be taken into account as it is expected to improve band profile modeling accuracy.
  Keywords: Molecular oxygen; Microwave spectroscopy; Profile shape modeling; Collisional coupling
  
  Note:  The ARTS implementation has not been tested in detail
  
  SIXJ(A,B,C,D,F) 
  C            |  A  B  1  |       A+B+C+D                                
  C            |  D  C  F  | = (-1)         W(ABCD;1F)                    
  C            |_         _|                           
 */
Numeric o2_ecs_wigner_symbol(
    int Nl, int Nk, int Jl, int Jk, int Jl_p, int Jk_p, int L) {
  return WIGNER3(Nl, Nk, L, 0, 0, 0) * WIGNER6(L, Jk, Jl, 2, Nl, Nk) *
         WIGNER6(L, Jk_p, Jl_p, 2, Nl, Nk) * WIGNER6(L, Jk, Jl, 2, Jl_p, Jk_p);
}

void ECS_wigner_CO2(Matrix& M,
                    const ArrayOfRational& Jl,
                    const ArrayOfRational& Ju,
                    const Rational& ll,
                    const Rational& lu,
                    ConstVectorView G0,
                    ConstVectorView population) {
  // Size of the problem
  const Index nj = Jl.nelem();
  M.resize(nj, nj);

  // wig3jj and wig6jj operates on integers
  int li = int((2 * ll).toIndex());
  int lf = int((2 * lu).toIndex());

  // Size of the wigner calculators memory requirements
  Index size = 0;
  for (Index i = 0; i < nj; i++) {
    Index tmp;
    tmp = (2 * Ju[i]).toIndex();
    if (tmp > size) size = tmp;
    tmp = (2 * Jl[i]).toIndex();
    if (tmp > size) size = tmp;
  }

  // Main loop over all the lines
  for (Index j = 0; j < nj; j++)  // For all lines
  {
    for (Index k = 0; k <= j; k++)  // For all lines up til now
    {
      if (j == k)
        M(j, k) = 2.0 * G0[j];
      else {
        // Only perform this step for downwards transitions
        const bool jbig = Jl[j] >= Jl[k];
        const Index big = jbig ? j : k, small = jbig ? k : j;

        // Conversion of type to fit with wigxjpf
        const int Ji = int((2 * Ju[big]).toIndex()),
                  Jf = int((2 * Jl[big]).toIndex()),
                  Ji_p = int((2 * Ju[small]).toIndex()),
                  Jf_p = int((2 * Jl[small]).toIndex());

        // Find potential start and end point of loop by relevance
        const int st = std::max(Ji - Ji_p, Jf - Jf_p),
                  en = std::min(Ji + Ji_p, Jf + Jf_p);

        // FIXME:  Adiabatic factor 1 and K1 goes here in relmat code...

        Numeric x = 0.0;

        // Loop over all relevant L: all the even numbers but not zero
        for (int L = st ? st : 4; L <= en; L += 4) {
          // FIXME:  Compute the basis rate (using coefficients or linearization)

          // FIXME:  Adiabatic factor 2 shuold be computed here if coefficients method
          const Numeric AF2 = 1.0;

          // This levels term... (n.b., there is 2L + 1, but sinc wigxjpf works on double digits, no need to write 2L)
          const Numeric y =
              co2_ecs_wigner_symbol(Ji, Jf, Ji_p, Jf_p, L, li, lf) / AF2;

          // Sum to the total
          x += y;
        }

        const Numeric rk = population[big] / population[small];
        M(big, small) = x;
        M(small, big) = rk * M(big, small);
      }
    }
  }
}

bool is_wigner_ready(int j) {
  extern int wigxjpf_max_prime_decomp;
  return not(j > wigxjpf_max_prime_decomp);
}

bool is_Wigner3_ready(const Rational& J) {
  const int test = (J * 6).toInt() / 2 + 1;  // nb. J can be half-valued
  return is_wigner_ready(test);
}

bool is_Wigner6_ready(const Rational& J) {
  const int test = (J * 4).toInt() + 1;  // nb. J can be half-valued
  return is_wigner_ready(test);
}
