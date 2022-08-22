/* Copyright (C) 2019
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

/*!
 * \file   arts_constexpr_math.h
 * \brief  Simple constexpr math that we can make use of in Arts
 *
 * The main advantage of maning things like here is to avoid repetition
 * later x*x*x*x is more difficult than pow4(x).
 * 
 * The main point of these functions is to avoid repetition as above in
 * context where we also want to support compile time computations
 * 
 * \author Richard Larsson
 * \date   2019-04-01
 */

#ifndef CONSTANT_MATHS_IN_ARTS_H
#define CONSTANT_MATHS_IN_ARTS_H

/** Namespace containing several constants, physical and mathematical **/
namespace Math {
/** power of two */
template <class T>
constexpr auto pow2(T x) noexcept -> decltype(x * x) {
  return x * x;
}

/** power of three */
template <class T>
constexpr auto pow3(T x) noexcept -> decltype(pow2(x) * x) {
  return pow2(x) * x;
}

/** power of four */
template <class T>
constexpr auto pow4(T x) noexcept -> decltype(pow2(pow2(x))) {
  return pow2(pow2(x));
}
}  // namespace Math

#endif 
