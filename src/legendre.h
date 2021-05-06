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
#include "matpackI.h"

typedef struct {
  size_t n;        /* number of points */
  double* x;       /* Gauss abscissae/points */
  double* w;       /* Gauss weights for each abscissae */
  int precomputed; /* high precision abscissae/weights precomputed? */
} gsl_integration_glfixed_table;

bool gsl_integration_glfixed_table_alloc(Vector& x, Vector& w, Index n);

namespace Legendre {
  struct SphericalField {
    Numeric U;
    Numeric S;
    Numeric E;
    Numeric total() const noexcept {return std::hypot(U, S, E);}
    constexpr SphericalField() noexcept : U(0), S(0), E(0) {}
  };
  
  using MatrixOfSphericalField = Grid<SphericalField, 2>;
  
  SphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Numeric r, const Numeric lat, const Numeric lon) ARTS_NOEXCEPT;
  
  MatrixOfSphericalField schmidt_fieldcalc(const Matrix& g, const Matrix& h, const Numeric r0, const Vector& r, const Numeric lat, const Vector& lon) ARTS_NOEXCEPT;
}

#endif /* legendre_h */
