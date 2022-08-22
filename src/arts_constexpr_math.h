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
 * \file   arts_constants.h
 * \brief  Constants of physical expressions as constexpr
 * 
 * Following May 2019 reported SI unit updates, several previously
 * important constant are proper constants and not derived constants.
 * 
 * This file contains the intended May 2019 constants update meaning
 * that the Planck constant, Boltzmann constant, elementary charge,
 * Avogadro number, the definition of a luminosity, the frequency of
 * Cesium, and the speed of light have been given constant values.
 * 
 * All other values are derived.  To be internally consistent in ARTS,
 * the adaptation of these rules below take the approach that only the
 * fine structure constant and the Rydberg constant are 'measurable'
 * quantities.  All other quantities are derived by their formal
 * relationships.  As an example, the resting mass of an electron is
 * defined as 2 h R_inf / c alpha^2.  This means that they should be
 * constant in their internal relationships to the level permitted by
 * the select floating point accuracy.
 * 
 * Note 1: The constants below contains both named and a few convenience
 * variables
 * 
 * Note 2: Derived constants of convenience, such as the Doppler broadening
 * constant, should also go into here if they have a clear name.
 * 
 * Note 3: The PrintPhysicalConstants convenience function is part of the 
 * global ARTS namespace.  Please update this if and when new constants
 * are added to this namespace.
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
