/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>

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
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-11-14   
  
  \brief  Header file for special_interp.cc.
*/

#ifndef special_interp_h
#define special_interp_h

#include "arts.h"
#include "interpolation.h"



void ArrayOfGridPosPrint(
        const ArrayOfGridPos&   x,
        const String&           x_name );



/*===========================================================================
  === Interpolation functions for atmospheric grids, fields and surfaces
  ===========================================================================*/

void interp_atm_field_gp2itw( 
              Matrix&           itw, 
        const Index&          	atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
	const ArrayOfGridPos&   gp_lon );

void interp_atmfield_by_itw( 
              VectorView        x, 
        const Index&          	atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstTensor3View        x_field,
 	const String&           x_field_name,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
	const ArrayOfGridPos&   gp_lon,
	ConstMatrixView         itw );

void interp_atmfield_by_gp( 
              VectorView        x, 
        const Index&          	atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstTensor3View        x_field,
 	const String&           x_field_name,
        const ArrayOfGridPos&   gp_p,
        const ArrayOfGridPos&   gp_lat,
	const ArrayOfGridPos&   gp_lon );

Numeric interp_atmfield_by_gp( 
        const Index&          	atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstTensor3View        x_field,
 	const String&           x_field_name,
        const GridPos&   	gp_p,
        const GridPos&   	gp_lat,
	const GridPos&   	gp_lon );

void interp_atmsurface_gp2itw( 
              Matrix&           itw, 
        const Index&          	atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        const ArrayOfGridPos&   gp_lat,
	const ArrayOfGridPos&   gp_lon );

void interp_atmsurface_by_itw(
              VectorView        x, 
        const Index&          	atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstMatrixView         x_surface,
 	const String&           x_surface_name,
        const ArrayOfGridPos&   gp_lat,
	const ArrayOfGridPos&   gp_lon,
        ConstMatrixView         itw );

void interp_atmsurface_by_gp( 
              VectorView        x, 
        const Index&          	atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstMatrixView         x_field,
 	const String&           x_field_name,
        const ArrayOfGridPos&   gp_lat,
	const ArrayOfGridPos&   gp_lon );

Numeric interp_atmsurface_by_gp( 
        const Index&          	atmosphere_dim,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
	ConstMatrixView         x_field,
 	const String&           x_field_name,
        const GridPos&   	gp_lat,
	const GridPos&   	gp_lon );

void itw2p(
              VectorView       p_values,
        ConstVectorView        p_grid,
	const ArrayOfGridPos&  gp,
	ConstMatrixView        itw );

void z_at_lat_2d(
	     VectorView   	 z,
        ConstVectorView   	 p_grid,
        ConstVectorView   	 lat_grid,
        ConstMatrixView          z_field,
	const ArrayOfGridPos&    gp_lat );


#endif // special_interp_h
