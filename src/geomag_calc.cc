/* Copyright (C) 2003-2012 Nikolay Koulev <nkoulev@uni-bremen.de>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   geomag_calc.cc
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Mon Jul 28 11:38:22 2003
  
  \brief  Routine for calculating the geomagnetic field.
  
  
*/
#include <cmath>
#include <iostream>
#include "arts.h"
#include "legendre.h"
#include "matpackI.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric EARTH_RADIUS;

void magfield_nk(   // Output
    Numeric& B_r,   // radial component of the geomagnetic field in [nT].
    Numeric& B_th,  // latitudinal component of the geomagnetic field in [nT].
    Numeric& B_ph,  // longitudinal component of the geomagnetic field in [nT].

    // Input
    const Numeric r,      // radial distance to the point in [km]
    const Numeric theta,  // geocentric colatitude of the point in [deg]
    const Numeric phi,    // longitude of the point in [deg].
    // All coordinates - geocentric!

    const Index Ny,  // number of elapsed years after an epoch year, J - [0,4]
    const Verbosity& verbosity)

{
  Numeric Phi = phi * DEG2RAD;               // Longitude angle in radian.
  Numeric Theta = PI / 2 - theta * DEG2RAD;  // Colatitude angle in radian.

  // Initializing values of the magnetic field components.
  B_r = 0;
  B_th = 0;
  B_ph = 0;

  Matrix M;

  xml_read_from_file("geomag_coefficients.xml", M, verbosity);

  // M(i,0) and M(i,1) - the vectors with the values of the first and second coefficients
  // of the IGRF model.
  // M(i,2) and M(i,3) - the vectors with the values of the anual rate of change of the
  // first and second coefficient of of the IGRF model.

  // Loop over the degree number l of the Legendre polynomes.
  for (Index l = 1; l <= 10; l++) {
    // Loop over the order number m of the Legendre polynomes.
    for (Index m = 0; m <= l; m++) {
      // Relating the row index in M to the coresponding
      // degree number l and order number m.
      Index j = l * (l + 1) / 2 + m - 1;

      // Calculating the associated Schmidt quasi-normalized Legendre
      // polynomial for a degree number l and order number m.
      Numeric P_lm = g_legendre_poly_norm_schmidt(l, m, cos(Theta));

      // Calculating the derivative of the associated Schmidt quasi-normalized
      // Legendre polynomial for a degree number l and order number m.
      Numeric dP_lm = g_legendre_poly_norm_schmidt_deriv3(l, m, cos(Theta));

      // Calculating the radial (upward) component of the magnetic field.
      B_r += pow((Numeric)(l + 2), EARTH_RADIUS / r) * (Numeric)(l + 1) *
             ((M(j, 0) + (Numeric)Ny * M(j, 2)) * cos((Numeric)m * Phi) +
              (M(j, 1) + (Numeric)Ny * M(j, 3)) * sin((Numeric)m * Phi)) *
             P_lm;

      // Calculating the latitudinal (southward) component of the magnetic field.
      B_th += pow(l + 2, EARTH_RADIUS / r) *
              ((M(j, 0) + (Numeric)Ny * M(j, 2)) * cos((Numeric)m * Phi) +
               (M(j, 1) + (Numeric)Ny * M(j, 3)) * sin((Numeric)m * Phi)) *
              dP_lm * sin(Theta);

      // Calculating the longitudinal (eastward) component of the magnetic field.
      B_ph += pow(l + 2, EARTH_RADIUS / r) * (Numeric)m *
              ((M(j, 0) + (Numeric)Ny * M(j, 2)) * sin((Numeric)m * Phi) -
               (M(j, 1) + (Numeric)Ny * M(j, 3)) * cos((Numeric)m * Phi)) *
              P_lm / sin(Theta);
    }
  }
}
