/* Copyright (C) 2000, 2001 Patrick Eriksson <patrick@rss.chalmers.se>
                            Stefan Buehler   <sbuehler@uni-bremen.de>

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


#include "matpackI.h"



////////////////////////////////////////////////////////////////////////////
//   Physical functions
////////////////////////////////////////////////////////////////////////////

void planck (
              MatrixView     B, 
        ConstVectorView     f,
        ConstVectorView     t );

void planck (
              VectorView     B, 
        ConstVectorView     f,
        const Numeric    t );

Numeric number_density (
       const Numeric   p,
       const Numeric   t );

Vector number_density (
       ConstVectorView    p,
       ConstVectorView    t );

Numeric g_of_z (
       const Numeric   r_geoid,
       const Numeric   g0,
       const Numeric   z );


////////////////////////////////////////////////////////////////////////////
//   Core functions for RTE and BL
////////////////////////////////////////////////////////////////////////////

void rte_iterate (
             VectorView   y,
       const Index   start_index,
       const Index   stop_index,
       ConstMatrixView   Tr,
       ConstMatrixView   S,
       const Index    n_f );

void rte (
             VectorView   y,
       const Index   start_index,
       const Index   stop_index,
       ConstMatrixView   Tr,
       ConstMatrixView   S,
       ConstVectorView   y_space,
       const Index    ground,
       ConstVectorView   e_ground,
       ConstVectorView   y_ground );

void bl_iterate (
             VectorView   y,
       const Index   start_index,
       const Index   stop_index,
       ConstMatrixView   Tr,
       const Index    n_f );

void bl (
             Vector&   y,
       const Index   start_index,
       const Index   stop_index,
       ConstMatrixView   Tr,
       const Index    ground,
       ConstVectorView   e_ground );



////////////////////////////////////////////////////////////////////////////
//   Conversion and interpolation of pressure and altitude grids.
////////////////////////////////////////////////////////////////////////////

void z2p(
              VectorView     p,
        ConstVectorView     z0,
        ConstVectorView     p0,
        ConstVectorView     z );

void interpp(
              VectorView     x,
        ConstVectorView     p0,
        ConstVectorView     x0,
        ConstVectorView     p );

void interpp_cloud(
              VectorView     x,
        ConstVectorView     p0,
        ConstVectorView     x0,
        ConstVectorView     p );

void interpp(
              MatrixView  A,
        ConstVectorView  p0, 
        ConstMatrixView  A0, 
        ConstVectorView  p );

Numeric interpp(
        ConstVectorView     p0,
        ConstVectorView     x0,
        const Numeric    p );

void interpz(
              VectorView     x, 
        ConstVectorView     p0,
        ConstVectorView     z0,
        ConstVectorView     x0,
        ConstVectorView     z );

Numeric interpz(
        ConstVectorView     p0,
        ConstVectorView     z0,
        ConstVectorView     x0,
        const Numeric    z );



////////////////////////////////////////////////////////////////////////////
//   Tangent altitudes
////////////////////////////////////////////////////////////////////////////

Numeric ztan_geom(
        const Numeric   za,
        const Numeric   z_plat,
        const Numeric   r_geoid );

Numeric n_for_z(
        const Numeric      z,
        ConstVectorView       p_abs,
        ConstVectorView       z_abs,
        ConstVectorView       refr_index,
        const Numeric      atm_limit );

Numeric refr_constant( 
        const Numeric      r_geoid,
        const Numeric      za,
        const Numeric      z_plat,
        ConstVectorView       p_abs,
        ConstVectorView       z_abs,
        const Numeric      atm_limit,
        ConstVectorView       refr_index );

Numeric ztan_refr(
        const Numeric   c,
        const Numeric   za,
        const Numeric   z_plat,
        const Numeric   z_ground,
        ConstVectorView    p_abs,
        ConstVectorView    z_abs,
        ConstVectorView    refr_index,
        const Numeric   r_geoid );

#endif // atmfuncs_h
