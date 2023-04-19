/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
 * @file   physics_funcs.h
 * @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
 * @date   2002-05-08
 *
 * @brief This file contains declerations of functions of physical character.
*/

#ifndef physics_h
#define physics_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts_conversions.h"
#include "arts.h"
#include "matpack_data.h"
#include "matpack_complex.h"

/*===========================================================================
  === Functions in physics_funcs.h
  ===========================================================================*/

Numeric barometric_heightformula(const Numeric& p, const Numeric& dh);

Numeric dinvplanckdI(const Numeric& i, const Numeric& f);

void fresnel(Complex& Rv,
             Complex& Rh,
             const Complex& n1,
             const Complex& n2,
             const Numeric& theta);

Numeric invplanck(const Numeric& i, const Numeric& f);

Numeric invrayjean(const Numeric& i, const Numeric& f);

/** number_density
 * 
 * Calculates the atmospheric number density.
 *
 * @param[in]  p  Pressure.
 * @param[in]  t  Temperature.
 *
 * @return     Number density.
 *
 * @author Patrick Eriksson
 * @date   2000-04-08
 */
constexpr Numeric number_density(Numeric p, Numeric t) noexcept {return p / (Constant::k * t);}

/** dnumber_density_dT
 * 
 * Calculates the atmospheric number density derivative with temperature.
 *
 * @param[in]  p  Pressure.
 * @param[in]  t  Temperature.
 *
 * @return     Number density.
 *
 * @author Richard Larsson
 * @date   2015-09-22
 */
constexpr Numeric dnumber_density_dt(Numeric p, Numeric t) noexcept {return - p / (Constant::k * Math::pow2(t));}

Numeric planck(const Numeric& f, const Numeric& t);

void planck(VectorView b, const ConstVectorView& f, const Numeric& t);

Vector planck(const ConstVectorView& f, const Numeric& t);

Numeric dplanck_dt(const Numeric& f, const Numeric& t);

void dplanck_dt(VectorView b, const ConstVectorView& f, const Numeric& t);

Vector dplanck_dt(const ConstVectorView& f, const Numeric& t);

Numeric dplanck_df(const Numeric& f, const Numeric& t);

Vector dplanck_df(const ConstVectorView& f, const Numeric& t);

Numeric rayjean(const Numeric& f, const Numeric& t);

#endif  // physics_h
