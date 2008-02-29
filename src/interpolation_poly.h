/* Copyright (C) 2008 Stefan Buehler <sbuehler(at)ltu.se>

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

/*!
  \file   interpolation_poly.h
  \author Stefan Buehler <sbuehler(at)ltu.se>
  \date   2008-02-04
  
  \brief  Header file for interpolation_poly.cc.
*/

#ifndef interpolation_poly_h
#define interpolation_poly_h

#include "matpackI.h"

//! Structure to store a grid position for higher order interpolation.
/*! 
  This serves the same purpose as GridPos for linear
  interpolation. The main difference is that we store not fractional
  distances (fd), but interpolation weights (w).

  There is some confusion between the meaning of fractional distances
  and interpolation weights. In fact the two are almost the same in the 1D
  case! (But the sorting and signs can be different.) 

  The w here correspond exactly to the terms in front of the yi in
  Numerical Recipes, 2nd edition, section 3.1, eq. 3.1.1.

  Only for 2D or higher dimensional interpolation are further
  calculations necessary, to multiply the w for the individual
  dimensions. 

  The size of w must equal the number of points used in the
  interpolation m. (m=2 For linear, m=3 for quadratic, etc..) 
*/
struct GridPosPoly {
  //! Index to first point to use for interpolation in the original grid.
  Index  idx;
  /*! Interpolation weight for each grid point to use.
  (Dimension is the number of points in the interpolation, m.)  */
  Vector w;
};

//! An Array of grid positions.
/*! 
  See \ref GridPosPoly for details.
*/

typedef Array<GridPosPoly> ArrayOfGridPosPoly;

void gridpos_poly(ArrayOfGridPosPoly& gp,
                  ConstVectorView old_grid,
                  ConstVectorView new_grid,
                  const Index order,
                  const Numeric  extpolfac = 0.5);

void interpweights_poly( Vector& itw,
                         const GridPosPoly& tc );

Numeric interp_poly( ConstVectorView    itw,
                     ConstVectorView    a,    
                     const GridPosPoly& tc );

void interpweights_poly( Matrix& itw,
                         const ArrayOfGridPosPoly& cgp );

void interp_poly( VectorView            ia,
                  ConstMatrixView       itw,
                  ConstVectorView       a,    
                  const ArrayOfGridPosPoly& cgp);

#endif // interpolation_poly_h
