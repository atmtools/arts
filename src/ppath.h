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



/*===========================================================================
  ===  File description 
  ===========================================================================*/

/*!
  \file   ppath.h
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-02
  
  \brief  Propagation path structure and functions.
  
   This file contains the definition of the Ppath structure and the
   functions in ppath.cc that are of interest elsewhere.
*/


#ifndef ppath_h
#define ppath_h


#include "array.h"
#include "arts.h"
#include "interpolation.h"
#include "matpackI.h"
#include "mystring.h"



/*===========================================================================
  === The Ppath structure
  ===========================================================================*/

//! The structure to describe a propagation path and releated quantities.
/*! 
   The fields of the structure are described in the ARTS user guide (AUG).
   It is listed as a sub-entry to "data structures".  
*/

struct Ppath {
  Index             dim;
  Index             np;
  Index             refraction;
  String            method;
  Numeric           constant;
  Matrix            pos;
  Vector            z;
  Vector            l_step;
  ArrayOfGridPos    gp_p;
  ArrayOfGridPos    gp_lat;
  ArrayOfGridPos    gp_lon;
  Matrix            los;
  String            background;
  Vector            tan_pos;
  Vector            geom_tan_pos;
};



/*===========================================================================
  === Functions from ppath.cc
  ===========================================================================*/

Numeric geometrical_ppc( const Numeric& r, const Numeric& za );

bool is_los_downwards( 
        const Numeric&   za,
        const Numeric&   tilt );

Numeric psurface_slope_2d(
        ConstVectorView   lat_grid,           
        ConstVectorView   r_geoid,
        ConstVectorView   z_surf,
        const GridPos&    gp,
        const Numeric&    za );

Numeric psurface_slope_3d(
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,  
        ConstMatrixView   r_geoid,
        ConstMatrixView   z_surf,
        const GridPos&    gp_lat,
        const GridPos&    gp_lon,
        const Numeric&    aa );

Numeric psurface_angletilt(
        const Numeric&   r,
        const Numeric&   c );

void do_gridcell_2d(
              Vector&    r_v,
              Vector&    lat_v,
              Vector&    za_v,
              Numeric&   lstep,
              Index&     endface,
        const Numeric&   r_start,
        const Numeric&   lat_start,
        const Numeric&   za_start,
        const Numeric&   ppc,
        const Numeric&   lmax,
        const Numeric&   r1,
        const Numeric&   r2,
        const Numeric&   r3,
        const Numeric&   r4,
        const Numeric&   lat1,
        const Numeric&   lat3,
        const Numeric&   rground1,
        const Numeric&   rground2 );

void do_gridcell_3d(
              Vector&    r_v,
              Vector&    lat_v,
              Vector&    lon_v,
              Vector&    za_v,
              Vector&    aa_v,
              Numeric&   lstep,
              Index&     endface,
        const Numeric&   r_start,
        const Numeric&   lat_start,
        const Numeric&   lon_start,
        const Numeric&   za_start,
        const Numeric&   aa_start,
        const Numeric&   ppc,
        const Numeric&   lmax,
        const Numeric&   r1a,
        const Numeric&   r2a,
        const Numeric&   r3a,
        const Numeric&   r4a,
        const Numeric&   r1b,
        const Numeric&   r2b,
        const Numeric&   r3b,
        const Numeric&   r4b,
        const Numeric&   lat1,
        const Numeric&   lat3,
        const Numeric&   lon5,
        const Numeric&   lon6,
        const Numeric&   rground1a,
        const Numeric&   rground2a,
        const Numeric&   rground1b,
	const Numeric&   rground2b );

void ppath_init_structure( 
              Ppath&      ppath,
        const Index&      atmosphere_dim,
        const Index&      np );

void ppath_set_background( 
              Ppath&      ppath,
        const Index&      case_nr );

Index ppath_what_background( const Ppath&   ppath );

void ppath_start_stepping(
              Ppath&          ppath,
        const Index&          atmosphere_dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        ConstTensor3View      z_field,
        ConstMatrixView       r_geoid,
        ConstMatrixView       z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        ConstVectorView       sensor_pos,
        ConstVectorView       sensor_los );

void ppath_step_geom_1d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        const Numeric&    r_geoid,
        const Numeric&    z_ground,
        const Numeric&    lmax );

void ppath_step_geom_2d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
        const Numeric&    lmax );

void ppath_step_geom_3d(
              Ppath&       ppath,
        ConstVectorView    p_grid,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        ConstTensor3View   z_field,
        ConstMatrixView    r_geoid,
        ConstMatrixView    z_ground,
	const Numeric&     lmax );

void ppath_step_refr_1d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        const Numeric&    r_geoid,
        const Numeric&    z_ground,
        const String&     rtrace_method,
        const Numeric&    lraytrace,
        const Numeric&    lmax,
        const String&     refrindex );

void ppath_step_refr_2d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
        const String&     rtrace_method,
        const Numeric&    lraytrace,
        const Numeric&    lmax,
        const String&     refrindex );

#endif  // ppath_h
