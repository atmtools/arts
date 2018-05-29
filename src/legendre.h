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
#include "matpackI.h"

typedef struct
{
    size_t n;         /* number of points */
    double *x;        /* Gauss abscissae/points */
    double *w;        /* Gauss weights for each abscissae */
    int precomputed;  /* high precision abscissae/weights precomputed? */
}
        gsl_integration_glfixed_table;

Numeric
legendre_poly (Index l, Index m, Numeric x);

Numeric
legendre_poly_norm_schmidt (Index l, Index m, Numeric x);

Numeric
legendre_poly_deriv (Index l, Index m, Numeric x);

Numeric
legendre_poly_norm_schmidt_deriv (Index l, Index m, Numeric x);

Numeric
g_legendre_poly (Index l, Index m, Numeric x);

Numeric
g_legendre_poly_norm_schmidt (Index l, Index m, Numeric x);

Numeric
g_legendre_poly_deriv (Index l, Index m, Numeric x);

Numeric
g_legendre_poly_norm_schmidt_deriv (Index l, Index m, Numeric x);

Numeric
g_legendre_poly_norm_schmidt_deriv1 (Index l, Index m, Numeric x);

Numeric
g_legendre_poly_norm_schmidt_deriv2 (Index l, Index m, Numeric x);

Numeric
g_legendre_poly_norm_schmidt_deriv3 (Index l, Index m, Numeric x);

Numeric
g_legendre_poly_norm_schmidt_deriv4 (Index l, Index m, Numeric x);

bool
gsl_integration_glfixed_table_alloc(Vector& x, Vector& w, Index n);

#endif  /* legendre_h */

