/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   atm_funcs.h

   This file contains declerations of functions releated to atmospheric 
   physics or geometry.

   \author Patrick Eriksson
   \date 2000-09-18 
*/


#ifndef atmfuncs_h
#define atmfuncs_h



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "vecmat.h"



////////////////////////////////////////////////////////////////////////////
//   Physical functions
////////////////////////////////////////////////////////////////////////////

void planck (
              Matrix&     B, 
        const Vector&     f,
        const Vector&     t );

void planck (
              Vector&     B, 
        const Vector&     f,
        const Numeric&    t );

Numeric number_density (
       const Numeric&   p,
       const Numeric&   t );

Vector number_density (
       const Vector&    p,
       const Vector&    t );

Numeric g_of_z (
       const Numeric&   r_geoid,
       const Numeric&   g0,
       const Numeric&   z );


////////////////////////////////////////////////////////////////////////////
//   Core functions for RTE and BL
////////////////////////////////////////////////////////////////////////////

void rte_iterate (
             Vector&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const Matrix&   Tr,
       const Matrix&   S,
       const size_t    n_f );

void rte (
             Vector&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const Matrix&   Tr,
       const Matrix&   S,
       const Vector&   y_space,
       const Index&    ground,
       const Vector&   e_ground,
       const Vector&   y_ground );

void bl_iterate (
             Vector&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const Matrix&   Tr,
       const size_t    n_f );

void bl (
             Vector&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const Matrix&   Tr,
       const Index&    ground,
       const Vector&   e_ground );



////////////////////////////////////////////////////////////////////////////
//   Conversion and interpolation of pressure and altitude grids.
////////////////////////////////////////////////////////////////////////////

void z2p(
              Vector&     p,
        const Vector&     z0,
        const Vector&     p0,
        const Vector&     z );

void interpp(
              Vector&     x,
        const Vector&     p0,
        const Vector&     x0,
        const Vector&     p );

void interpp_cloud(
              Vector&     x,
        const Vector&     p0,
        const Vector&     x0,
        const Vector&     p );

void interpp(
              Matrix&  A,
        const Vector&  p0, 
        const Matrix&  A0, 
        const Vector&  p );

Numeric interpp(
        const Vector&     p0,
        const Vector&     x0,
        const Numeric&    p );

void interpz(
              Vector&     x, 
        const Vector&     p0,
        const Vector&     z0,
        const Vector&     x0,
        const Vector&     z );

Numeric interpz(
        const Vector&     p0,
        const Vector&     z0,
        const Vector&     x0,
        const Numeric&    z );



////////////////////////////////////////////////////////////////////////////
//   Tangent altitudes
////////////////////////////////////////////////////////////////////////////

Numeric ztan_geom(
        const Numeric&   za,
        const Numeric&   z_plat,
        const Numeric&   r_geoid );

Numeric n_for_z(
        const Numeric&      z,
        const Vector&       p_abs,
        const Vector&       z_abs,
        const Vector&       refr_index,
        const Numeric&      atm_limit );

Numeric refr_constant( 
        const Numeric&      r_geoid,
        const Numeric&      za,
        const Numeric&      z_plat,
        const Vector&       p_abs,
        const Vector&       z_abs,
        const Numeric&      atm_limit,
        const Vector&       refr_index );

Numeric ztan_refr(
        const Numeric&   c,
        const Numeric&   za,
        const Numeric&   z_plat,
        const Numeric&   z_ground,
        const Vector&    p_abs,
        const Vector&    z_abs,
        const Vector&    refr_index,
        const Numeric&   r_geoid );

#endif // atmfuncs_h
