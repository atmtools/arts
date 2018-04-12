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
              VectorView      itw, 
              GridPos&        gp_p_out,
              GridPos&        gp_lat_out,
              GridPos&        gp_lon_out,
        const GridPos&        gp_p_in,
        const GridPos&        gp_lat_in,
        const GridPos&        gp_lon_in,
        const Index&          atmosphere_dim,
        const ArrayOfIndex&   cloudbox_limits );

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

void regrid_atmfield_by_gp( 
         Tensor3&          field_new, 
   const Index&            atmosphere_dim, 
   ConstTensor3View        field_old, 
   const ArrayOfGridPos&   gp_p,
   const ArrayOfGridPos&   gp_lat,
   const ArrayOfGridPos&   gp_lon );

void regrid_atmsurf_by_gp( 
         Matrix&           field_new, 
   const Index&            atmosphere_dim, 
   ConstMatrixView         field_old, 
   const ArrayOfGridPos&   gp_lat,
   const ArrayOfGridPos&   gp_lon );

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

void rte_pos2gridpos(
         GridPos&     gp_p,
         GridPos&     gp_lat,
         GridPos&     gp_lon,
   const Index&       atmosphere_dim,
   ConstVectorView    p_grid,
   ConstVectorView    lat_grid,
   ConstVectorView    lon_grid,
   ConstTensor3View   z_field,
   ConstVectorView    rte_pos );

void rte_pos2gridpos(
         GridPos&     gp_lat,
         GridPos&     gp_lon,
   const Index&       atmosphere_dim,
   ConstVectorView    lat_grid,
   ConstVectorView    lon_grid,
   ConstVectorView    rte_pos );

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

void complex_n_interp(
         MatrixView       n_real,
         MatrixView       n_imag,
   const GriddedField3&   complex_n,
   const String&          varname,
   ConstVectorView        f_grid,
   ConstVectorView        t_grid );

#endif // special_interp_h
