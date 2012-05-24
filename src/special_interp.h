/* Copyright (C) 2002-2012 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
  \file   special_interp.h
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-11-14   
  
  \brief  Header file for special_interp.cc.
*/

#ifndef special_interp_h
#define special_interp_h

#include "interpolation.h"
#include "interpolation_poly.h"
#include "gridded_fields.h"




/*===========================================================================
  === Interpolation functions for atmospheric grids, fields and surfaces
  ===========================================================================*/

void interp_atmfield_gp2itw( 
              Matrix&           itw, 
        const Index&            atmosphere_dim,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon );

void interp_atmfield_by_itw( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstTensor3View        x_field,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon,
        ConstMatrixView         itw );

void interp_atmfield_by_gp( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstTensor3View        x_field,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon );

Numeric interp_atmfield_by_gp( 
        const Index&            atmosphere_dim,
        ConstTensor3View        x_field,
        const GridPos&          gp_p,
        const GridPos&          gp_lat,
        const GridPos&          gp_lon );

void interp_cloudfield_gp2itw( 
              Matrix&           itw, 
              ArrayOfGridPos&   gp_p_out,
              ArrayOfGridPos&   gp_lat_out,
              ArrayOfGridPos&   gp_lon_out,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon,
        const Index&            atmosphere_dim,
        const ArrayOfIndex&     cloudbox_limits );

void interp_atmsurface_gp2itw( 
              Matrix&           itw, 
        const Index&            atmosphere_dim,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon );

void interp_atmsurface_by_itw(
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstMatrixView         x_surface,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon,
        ConstMatrixView         itw );

void interp_atmsurface_by_gp( 
              VectorView        x, 
        const Index&            atmosphere_dim,
        ConstMatrixView         x_field,
        const ArrayOfGridPos&   gp_lat,
        const ArrayOfGridPos&   gp_lon );

Numeric interp_atmsurface_by_gp( 
        const Index&            atmosphere_dim,
        ConstMatrixView         x_field,
        const GridPos&          gp_lat,
        const GridPos&          gp_lon );

void interp_gfield3( 
                 Numeric&   value,
           const GriddedField3&   gfield3,
           const Index&     effective_dim,
           const Numeric&   x,
           const Numeric&   y,
           const Numeric&   z,
           const String&    dim0,
           const String&    dim1,
           const String&    dim2 );

void itw2p(
              VectorView       p_values,
        ConstVectorView        p_grid,
        const ArrayOfGridPos&  gp,
        ConstMatrixView        itw );

void p2gridpos(
             ArrayOfGridPos&   gp,
             ConstVectorView   old_pgrid,
             ConstVectorView   new_pgrid,   
             const Numeric&    extpolfac=0.5 );

void p2gridpos_poly(
               ArrayOfGridPosPoly&   gp,
               ConstVectorView       old_pgrid,
               ConstVectorView       new_pgrid,   
               const Index           order,
               const Numeric&        extpolfac=0.5 );

void z_at_lat_2d(
             VectorView   z,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        const GridPos&    gp_lat );

void z_at_latlon(
             VectorView    z,
        ConstVectorView    p_grid,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        ConstTensor3View   z_field,
        const GridPos&     gp_lat,
        const GridPos&     gp_lon );


#endif // special_interp_h
