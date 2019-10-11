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

#include "arts.h"
#include "complex.h"
#include "matpackI.h"

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

Numeric number_density(const Numeric& p, const Numeric& t);

Numeric dnumber_density_dt(const Numeric& p, const Numeric& t);

Numeric planck(const Numeric& f, const Numeric& t);

void planck(VectorView b, ConstVectorView f, const Numeric& t);

Numeric dplanck_dt(const Numeric& f, const Numeric& t);

void dplanck_dt(VectorView b, ConstVectorView f, const Numeric& t);

Numeric dplanck_df(const Numeric& f, const Numeric& t);

Numeric rayjean(const Numeric& f, const Numeric& t);

#ifdef ENABLE_REFICE
extern "C" {
#endif
/** Calculates complex refractive index of Ice 1H
 *
 * Valid range for wavelengths is between 45 nm and 8.6 m.
 * For wavelengths above 167 microns, temperature dependence is
 * included for temperatures between 213 and 272K.
 * Mainly intended for applications in Earth ice
 * clouds and snow, not other planets or interstellar space;
 * the temperature dependence or crystalline form of ice may be
 * incorrect for these latter applications.
 *
 * Authors of Fortran function:<br>
 * Stephen Warren, Univ. of Washington (1983)<br>
 * Bo-Cai Gao, JCESS, Univ. of Maryland (1995)<br>
 * Warren Wiscombe, NASA Goddard (1995)
 *
 * References:<br>
 * Warren, S., 1984: Optical Constants of Ice from the Ultraviolet
 * to the Microwave, Appl. Opt. 23, 1206-1225
 *
 * Kou, L., D. Labrie, and P. Chylek, 1994: Refractive indices
 * of water and ice in the 0.65- to 2.5-micron spectral range,
 * Appl. Opt. 32, 3531-3540
 *
 * Perovich, D., and J. Govoni, 1991: Absorption Coefficients
 * of Ice from 250 to 400 nm, Geophys. Res. Lett. 18, 1233-1235
 *
 * @param[in] wavlen Wavelength in microns
 * @param[in] temp   Temperature in Kelvin
 *
 * returns Refractive index as complex number
 */
Complex refice_(const Numeric& wavlen, const Numeric& temp);
#ifdef ENABLE_REFICE
}
#endif

#endif  // physics_h
