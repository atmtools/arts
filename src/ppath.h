/* Copyright (C) 2002-2008 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-02
  
  \brief  Propagation path structure and functions.
  
   This file contains the definition of the Ppath structure and the
   functions in ppath.cc that are of interest elsewhere.
*/


#ifndef ppath_h
#define ppath_h


#include "agenda_class.h"
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
  Numeric           constant;
  String            background;
  Matrix            pos;
  Matrix            los;
  Vector            r;
  Vector            l_step;
  Numeric           lspace;
  Vector            nreal;
  ArrayOfGridPos    gp_p;
  ArrayOfGridPos    gp_lat;
  ArrayOfGridPos    gp_lon;
};


/** An array of propagation paths. */
typedef Array<Ppath> ArrayOfPpath;



/*===========================================================================
  === Common precision variables
  ===========================================================================*/

// This variable defines the maximum allowed error tolerance for radius.
// The variable is, for example, used to check that a given a radius is
// consistent with the specified grid cell.
//
#ifdef USE_DOUBLE
const double   RTOL = 1e-2;
#else
const double   RTOL = 10;
#endif


// As RTOL but for latitudes and longitudes.
//
#ifdef USE_DOUBLE
const double   LATLONTOL = 1e-11;
#else
const double   LATLONTOL = 1e-6;
#endif


// This variable defines how much zenith and azimuth angles can
// deviate from 0, 90 and 180 degrees, but still be treated to be 0,
// 90 or 180.  For example, an azimuth angle of 180-ANGTOL/2 will
// be treated as a strictly southward observation.  However, the
// angles are not allowed to go outside their defined range.  This
// means, for example, that values above 180 are never allowed.
//
#ifdef USE_DOUBLE
const double   ANGTOL = 1e-7; 
#else
const double   ANGTOL = 1e-4; 
#endif


// Latitudes with an absolute value > POLELAT are considered to be on
// the south or north pole for 3D.
//
const double   POLELAT = 89.9999;


// Maximum tilt of pressure levels, in degrees
//
const double    PTILTMAX = 5;



/*===========================================================================
  === Functions from ppath.cc
  ===========================================================================*/

void map_daa(
             double&   za,
             double&   aa,
       const double&   za0,
       const double&   aa0,
       const double&   aa_grid );

double geometrical_ppc( const double& r, const double& za );

double geompath_za_at_r(
        const double&   ppc,
        const double&   a_za,
        const double&   r );

double geompath_lat_at_za(
        const double&   za0,
        const double&   lat0,
        const double&   za );

bool is_los_downwards( 
        const double&   za,
        const double&   tilt );

double plevel_slope_2d(
        ConstVectorView   lat_grid,           
        ConstVectorView   refellipsoid,
        ConstVectorView   z_surf,
        const GridPos&    gp,
        const double&     za );

double plevel_slope_3d(
        const double&   lat1,
        const double&   lat3,
        const double&   lon5,
        const double&   lon6,
        const double&   r15,
        const double&   r35,
        const double&   r36,
        const double&   r16,
        const double&   lat,
        const double&   lon,
        const double&   aa );

double plevel_slope_3d(
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,  
        ConstVectorView   refellipsoid,
        ConstMatrixView   z_surf,
        const GridPos&    gp_lat,
        const GridPos&    gp_lon,
        const double&     aa );

double plevel_angletilt(
        const double&   r,
        const double&   c );

void ppath_init_structure( 
              Ppath&      ppath,
        const Index&      atmosphere_dim,
        const Index&      np );

void ppath_set_background( 
              Ppath&      ppath,
        const Index&      case_nr );

Index ppath_what_background( const Ppath&   ppath );

void ppath_copy(
              Ppath&      ppath1,
        const Ppath&      ppath2 );

void ppath_start_stepping(
              Ppath&            ppath,
        const Index&            atmosphere_dim,
        ConstVectorView         p_grid,
        ConstVectorView         lat_grid,
        ConstVectorView         lon_grid,
        ConstTensor3View        z_field,
        ConstVectorView         refellipsoid,
        ConstMatrixView         z_surface,
        const Index &           cloudbox_on,
        const ArrayOfIndex &    cloudbox_limits,
        const bool &            outside_cloudbox,
        ConstVectorView         rte_pos,
        ConstVectorView         rte_los,
        const Verbosity&        verbosity);


void ppath_step_geom_1d(
              Ppath&      ppath,
        ConstVectorView   z_field,
        ConstVectorView   refellipsoid,
        const double&     z_surface,
        const double&     lmax );

void ppath_geom_updown_1d(
              Ppath&      ppath,
        ConstVectorView   z_field,
        ConstVectorView   refellipsoid,
        const double&     z_surface );

void ppath_step_geom_2d(
              Ppath&      ppath,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   refellipsoid,
        ConstVectorView   z_surface,
        const double&     lmax );

void ppath_step_geom_3d(
              Ppath&       ppath,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        ConstTensor3View   z_field,
        ConstVectorView    refellipsoid,
        ConstMatrixView    z_surface,
        const double&      lmax );

void ppath_step_refr_1d(
              Workspace&  ws,
              Ppath&      ppath,
              Numeric&    rte_pressure,
              Numeric&    rte_temperature,
              Vector&     rte_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        ConstMatrixView   vmr_field,
        ConstVectorView   refellipsoid,
        const double&     z_surface,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax );

void ppath_step_refr_2d(
              Workspace&  ws,
              Ppath&      ppath,
              Numeric&    rte_pressure,
              Numeric&    rte_temperature,
              Vector&     rte_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstTensor3View  vmr_field,
        ConstVectorView   refellipsoid,
        ConstVectorView   z_surface,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax );

void ppath_step_refr_3d(
              Workspace&  ws,
              Ppath&      ppath,
              Numeric&    rte_pressure,
              Numeric&    rte_temperature,
              Vector&     rte_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstTensor3View  z_field,
        ConstTensor3View  t_field,
        ConstTensor4View  vmr_field,
        ConstVectorView   refellipsoid,
        ConstMatrixView   z_surface,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax );

void ppath_calc(
              Workspace&      ws,
              Ppath&          ppath,
        const Agenda&         ppath_step_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Vector&         refellipsoid,
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         rte_pos,
        const Vector&         rte_los,
        const bool&           outside_cloudbox,
        const Verbosity&      verbosity);

#endif  // ppath_h
