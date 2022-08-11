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
 * @file wigner_functions.h
 * @author Richard Larsson
 * @date 2013-06-19
 * 
 * @brief Wigner symbol interactions
 */

#ifndef wigner_functions_h
#define wigner_functions_h

#include "rational.h"

#include <wigner/wigxjpf/inc/wigxjpf.h>

#include <algorithm>
#include <array>

#ifdef FAST_WIGNER_PATH_3J
#define DO_FAST_WIGNER 1
#include <wigner/fastwigxj/inc/fastwigxj.h>
#else
#define DO_FAST_WIGNER 0
#endif


/*! Refine limits from multiple inputs
 * 
 * @param[in] a A limit
 * @param[in] a Another limit
 * @return [low, high] for valid ranges of numbers or two undefined rationals
 */
std::pair<Rational, Rational> wigner_limits(std::pair<Rational, Rational> a, std::pair<Rational, Rational> b);

/*! Return the limits where a wigner3j symbol can be non-zero
 * 
 * The values a-e are as in a call to wigner3j with pos
 * determining the offset.  The output is the range of
 * valid values at the offset position given the relevant
 * triangle equaility |x - y| <= z <= x + y and that the
 * m1 + m2 = -m3 condition for the lower row values
 * 
 * Positional information:
 *    pos == 1:
 *        wigner3j(X, a, b, c, d, e):
 *            b <= X + a -> X <= b - a, and
 *            |X - a| <= b -> a - b <= X <= a + b
 *    pos == 2:
 *        wigner3j(a, X, b, c, d, e):
 *            b <= a + X -> X <= b - a, and
 *            |a - X| <= b -> a + b <= X <= a - b
 *    pos == 3:
 *        wigner3j(a, b, X, c, d, e):
 *            |a - b| <= X <= a + b -> X <= b - a
 *    pos == 4:
 *        wigner3j(a, b, c, X, d, e):
 *            -|a| <= X <= |a|, and
 *            X + d = - e -> X = - e - d
 *    pos == 5:
 *        wigner3j(a, b, c, d, X, e):
 *            -|b| <= X <= |b|, and
 *            d + X = - e -> X = - e - d
 *    pos == 6:
 *        wigner3j(a, b, c, d, e, X):
 *            -|c| <= X <= |c|, and
 *            d + e = - X -> X = - e - d
 * 
 * If there is no valid range, the function returns
 * {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED}
 * 
 * @param[in] pos Position of value
 * @param[in] a An input
 * @param[in] b An input
 * @param[in] c An input
 * @param[in] d An input
 * @param[in] e An input
 * @return A valid range where both start and end are valid, or invalid numbers
 */
template<Index pos> constexpr
std::pair<Rational, Rational> wigner3j_limits([[maybe_unused]] const Rational a=0,
                                              [[maybe_unused]] const Rational b=0,
                                              [[maybe_unused]] const Rational c=0,
                                              [[maybe_unused]] const Rational d=0,
                                              [[maybe_unused]] const Rational e=0) {
  static_assert(pos < 7 and pos > 0, "Only valid for pos := 1, 2, 3, 4, 5, and 6");
  
  if constexpr (pos == 1 or pos == 2) {
    const Rational maxX = b - a;
    std::pair<Rational, Rational> out {-maxX, a + b};
    if (out.first > out.second) swap(out.first, out.second);
    if (out.second > maxX) out.second = maxX;
    if (out.first > maxX) out = {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED};
    return out;
  } else if constexpr (pos == 3) {
    const Rational maxX = a + b;
    const Rational minX = abs(a - b);
    if (maxX >= minX) return {minX, maxX};
    return {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED};
  } else {
    const Rational lim = pos == 4 ? abs(a) :
                         pos == 5 ? abs(b) :
                       /*pos == 6*/ abs(c);
    const Rational val = - e - d;
    if (-lim <= val and val <= lim) return {val, val};
    return {RATIONAL_UNDEFINED, RATIONAL_UNDEFINED};
  }
}

/** Wigner 3J symbol
  * 
  * Run wigxjpf wig3jj for Rational symbol
  * 
  * /                \
  * |  j1   j2   j3  |
  * |                |
  * |  m1   m2   m3  |
  * \                /
  * 
  * See for definition: http://dlmf.nist.gov/34.2
  * 
  * @param[in] j1 as above
  * @param[in] j2 as above
  * @param[in] j3 as above
  * @param[in] m1 as above
  * @param[in] m2 as above
  * @param[in] m3 as above
  * @return Numeric Symbol value
  */
Numeric wigner3j(const Rational j1,
                 const Rational j2,
                 const Rational j3,
                 const Rational m1,
                 const Rational m2,
                 const Rational m3);

/** Wigner 6J symbol
  * 
  * Run wigxjpf wig6jj for Rational symbol
  *
  * /                \
  * |  j1   j2   j3  |
  * <                >
  * |  l1   l2   l3  |
  * \                /
  *
  * See for definition: http://dlmf.nist.gov/34.4
  * 
  * @param[in] j1 as above
  * @param[in] j2 as above
  * @param[in] j3 as above
  * @param[in] l1 as above
  * @param[in] l2 as above
  * @param[in] l3 as above
  * @return Numeric Symbol value
  */
Numeric wigner6j(const Rational j1,
                 const Rational j2,
                 const Rational j3,
                 const Rational l1,
                 const Rational l2,
                 const Rational l3);

/** Ready Wigner
 * 
 * @param[in] largest
 * @param[in] fastest
 * @param[in] size [3 or 6]
 * @return largest if successful
 */
Index make_wigner_ready(int largest, int fastest, int size);

/** Tells if the function can deal with the input integer
 * 
 * @param[in] j 
 * @return true If j is less than max allowed j
 * @return false Otherwise
 */
bool is_wigner_ready(int j);

/** Tells if the function is ready for Wigner 3J calculations
 * 
 * @param[in] J Largest input into a Wigner 3J function call
 * @return true If is_wigner_ready(3J + 1) does
 * @return false Otherwise
 */
bool is_wigner3_ready(const Rational& J);

/** Tells if the function is ready for Wigner 6J calculations
 * 
 * @param[in] J Largest input into a Wigner 6J function call
 * @return true If is_wigner_ready(4J + 1)
 * @return false Otherwise
 */
bool is_wigner6_ready(const Rational& J);

template <class ... Integer>
int temp_init_size(Integer... vals) {
  constexpr auto N = sizeof...(Integer);
  static_assert(N > 0);
  const std::array<int, N> v{int(vals)...};
  return 1 + 3 * (*std::max_element(v.begin(), v.end()));
}

/** Computes the wigner 3J symbol with floating point precision

         /               \
         |  J1   J2   J  |
output = |               |
         |  M   -M    0  |
         \               /

 @param M Input as above
 @param J1 Input as above
 @param J2 Input as above
 @param J Input as above
 @return Numeric 
 */
Numeric dwigner3j(Index M, Index J1, Index J2, Index J);

/** Computes the wigner 6J symbol with floating point precision

         /             \
         |  A   B   1  |
output = <             >
         |  D   C   F  |
         \             /

  @param A Input as above
  @param B Input as above
  @param C Input as above
  @param D Input as above
  @param F Input as above
  @return Numeric 
 */
Numeric dwigner6j(Index A, Index B, Index C, Index D, Index F);

#endif  // wigner_functions_h
