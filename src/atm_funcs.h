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
              MATRIX&     B, 
        const VECTOR&     f,
        const VECTOR&     t );

void planck (
              VECTOR&     B, 
        const VECTOR&     f,
        const Numeric&    t );

Numeric number_density (
       const Numeric&   p,
       const Numeric&   t );

VECTOR number_density (
       const VECTOR&    p,
       const VECTOR&    t );

Numeric g_of_z (
       const Numeric&   r_geoid,
       const Numeric&   g0,
       const Numeric&   z );


////////////////////////////////////////////////////////////////////////////
//   Core functions for RTE and BL
////////////////////////////////////////////////////////////////////////////

void rte_iterate (
             VECTOR&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const MATRIX&   Tr,
       const MATRIX&   S,
       const size_t    n_f );

void rte (
             VECTOR&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const MATRIX&   Tr,
       const MATRIX&   S,
       const VECTOR&   y_space,
       const INDEX&    ground,
       const VECTOR&   e_ground,
       const VECTOR&   y_ground );

void bl_iterate (
             VECTOR&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const MATRIX&   Tr,
       const size_t    n_f );

void bl (
             VECTOR&   y,
       const size_t&   start_index,
       const size_t&   stop_index,
       const MATRIX&   Tr,
       const INDEX&    ground,
       const VECTOR&   e_ground );



////////////////////////////////////////////////////////////////////////////
//   Conversion and interpolation of pressure and altitude grids.
////////////////////////////////////////////////////////////////////////////

void z2p(
              VECTOR&     p,
        const VECTOR&     z0,
        const VECTOR&     p0,
        const VECTOR&     z );

void interpp(
              VECTOR&     x,
        const VECTOR&     p0,
        const VECTOR&     x0,
        const VECTOR&     p );

void interpp_cloud(
              VECTOR&     x,
        const VECTOR&     p0,
        const VECTOR&     x0,
        const VECTOR&     p );

void interpp(
              MATRIX&  A,
        const VECTOR&  p0, 
        const MATRIX&  A0, 
        const VECTOR&  p );

Numeric interpp(
        const VECTOR&     p0,
        const VECTOR&     x0,
        const Numeric&    p );

void interpz(
              VECTOR&     x, 
        const VECTOR&     p0,
        const VECTOR&     z0,
        const VECTOR&     x0,
        const VECTOR&     z );

Numeric interpz(
        const VECTOR&     p0,
        const VECTOR&     z0,
        const VECTOR&     x0,
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
        const VECTOR&       p_abs,
        const VECTOR&       z_abs,
        const VECTOR&       refr_index,
        const Numeric&      atm_limit );

Numeric refr_constant( 
        const Numeric&      r_geoid,
        const Numeric&      za,
        const Numeric&      z_plat,
        const VECTOR&       p_abs,
        const VECTOR&       z_abs,
        const Numeric&      atm_limit,
        const VECTOR&       refr_index );

Numeric ztan_refr(
        const Numeric&   c,
        const Numeric&   za,
        const Numeric&   z_plat,
        const Numeric&   z_ground,
        const VECTOR&    p_abs,
        const VECTOR&    z_abs,
        const VECTOR&    refr_index,
        const Numeric&   r_geoid );

#endif // atmfuncs_h
