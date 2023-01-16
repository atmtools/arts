/* Copyright (C) 2003-2012 Oliver Lemke  <olemke@core-dump.info>

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
   USA.
 */

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
  \file   legendre.h

  Contains the code to calculate Legendre polynomials.

  \author Oliver Lemke
  \date 2003-08-14
  */

#ifndef legendre_h
#define legendre_h

#include "arts.h"
#include "grids.h"
#include "matpack_data.h"

struct gsl_integration_glfixed_table {
  size_t n;        /* number of points */
  double* x;       /* Gauss abscissae/points */
  double* w;       /* Gauss weights for each abscissae */
  int precomputed; /* high precision abscissae/weights precomputed? */
};

bool gsl_integration_glfixed_table_alloc(Vector& x, Vector& w, Index n);


namespace Legendre {

//! Stores the up (U), south (S), east (E) values of a field relative to the sphere
struct SphericalField {
  Numeric U{0};
  Numeric S{0};
  Numeric E{0};
  
  //! Returns the total strength of the field
  Numeric total() const noexcept {return std::hypot(U, S, E);}
  
  //! Returns the total strength of the field
  Numeric total_horizontal() const noexcept {return std::hypot(S, E);}
  
  //! Always construct to zeroes explicitly
  constexpr SphericalField() noexcept = default;
};


//! Holds a SphericalField for multiple dimensions (radius times longitudes)
using MatrixOfSphericalField = Grid<SphericalField, 2>;


/** Computes the spherical field
 * 
 * If latitude is beyond a limit, the longitude is set to zero
 * and the latitude at the limit away from the pole is used. For these
 * position, the absolute of the horizontal component is kept around in
 * the south-facing component and the east-facing component is set to zero
 * 
 * The latitude limit is defined in the ColatitudeConversion struct.
 * 
 * The longitude is constrained to the range [-180, 180) to have consistent
 * behavior at the daytime border.
 * 
 * This function is taken from the implementation by Isabela de Oliveira Martins
 * at https://github.com/de-oliveira/IsabelaFunctions/blob/master/IsabelaFunctions/fieldmodel.py
 * (2021-05-06).
 * 
 * @param[in] g A N x N matrix of g-coefficients
 * @param[in] h A N x N matrix of h-coefficients
 * @param[in] r0 The reference radius (spherical)
 * @param[in] r The actual radius (spherical)
 * @param[in] lat The latitude (spherical)
 * @param[in] lon The longitude (spherical)
 * @return A spherical field
 */
SphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Numeric r, const Numeric lat, const Numeric lon) ARTS_NOEXCEPT;


/** Computes the spherical field for many radius and longitudes
 * 
 * See purely Numeric implementation for more information.
 * 
 * @param[in] g A N x N matrix of g-coefficients
 * @param[in] h A N x N matrix of h-coefficients
 * @param[in] r0 The reference radius (spherical)
 * @param[in] r The actual radius (spherical)
 * @param[in] lat The latitude (spherical)
 * @param[in] lon The longitude (spherical)
 * @return A spherical field
 */
MatrixOfSphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Vector& r, const Numeric lat, const Vector& lon) ARTS_NOEXCEPT;
} // namespace Legendre

#endif /* legendre_h */
