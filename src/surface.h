/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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

/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
   \file   surface.h
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date   2012-02-06 

   This file contains definitions of internal functions associated with 
   the surface.
*/

#ifndef surface_h
#define surface_h

#include "complex.h"
#include "matpackIV.h"
#include "mystring.h"
#include "ppath.h"

Numeric calc_incang(ConstVectorView rte_los, ConstVectorView specular_los);

void surface_calc(Matrix& iy,
                  ConstTensor3View I,
                  ConstMatrixView surface_los,
                  ConstTensor4View surface_rmatrix,
                  ConstMatrixView surface_emission);

void surface_specular_R_and_b(MatrixView surface_rmatrix,
                              VectorView surface_emission,
                              const Complex& Rv,
                              const Complex& Rh,
                              const Numeric& f,
                              const Index& stokes_dim,
                              const Numeric& surface_skin_t);

void surface_props_check(const Index& atmosphere_dim,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& surface_props_data,
                         const ArrayOfString& surface_props_names);

void surface_props_interp(Vector& v,
                          const String& vname,
                          const Index& atmosphere_dim,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const Matrix& itw,
                          const Tensor3& surface_props_data,
                          const ArrayOfString& surface_props_names);

void dsurface_check(const ArrayOfString& surface_props_names,
                    const ArrayOfString& dsurface_names,
                    const ArrayOfTensor4 dsurface_rmatrix_dx,
                    const ArrayOfMatrix& dsurface_emission_dx);

#endif  // surface_h
