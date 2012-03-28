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
   \file   geodetic.h
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date   2012-02-06 

   This file contains definitions of internal functions of geodetic character.
*/



#ifndef geodetic_h
#define geodetic_h

#include "interpolation.h"
#include "matpackI.h"

// 2D:

void cart2pol(
            double&   r,
            double&   lat,
      const double&   x,
      const double&   z,
      const double&   lat0,
      const double&   za0 );

void cart2poslos(
             double&   r,
             double&   lat,
             double&   za,
       const double&   x,
       const double&   z,
       const double&   dx,
       const double&   dz,
       const double&   ppc,
       const double&   lat0,
       const double&   za0 );

void distance2D(
            double&   l,
      const double&   r1,
      const double&   lat1,
      const double&   r2,
      const double&   lat2 );

/*
void geomtanpoint2d( 
             double&    r_tan,
             double&    lat_tan,
     ConstVectorView    refellipsoid,
       const double&    r,
       const double&    lat,
       const double&    za );
*/

void pol2cart(
            double&   x,
            double&   z,
      const double&   r,
      const double&   lat );

void poslos2cart(
             double&   x,
             double&   z,
             double&   dx,
             double&   dz,
       const double&   r,
       const double&   lat,
       const double&   za );



// 3D:

void cart2poslos(
             double&   r,
             double&   lat,
             double&   lon,
             double&   za,
             double&   aa,
       const double&   x,
       const double&   y,
       const double&   z,
       const double&   dx,
       const double&   dy,
       const double&   dz,
       const double&   ppc,
       const double&   lat0,
       const double&   lon0,
       const double&   za0,
       const double&   aa0 );

void cart2sph(
             double&   r,
             double&   lat,
             double&   lon,
       const double&   x,
       const double&   y,
       const double&   z,
       const double&   lat0,
       const double&   lon0,
       const double&   za0,
       const double&   aa0 );

void distance3D(
            double&   l,
      const double&   r1,
      const double&   lat1,
      const double&   lon1,
      const double&   r2,
      const double&   lat2,
      const double&   lon2 );

void geompath_tanpos_3d( 
             double&    r_tan,
             double&    lat_tan,
             double&    lon_tan,
             double&    l_tan,
       const double&    r,
       const double&    lat,
       const double&    lon,
       const double&    za,
       const double&    aa,
       const double&    ppc );

/*
void geomtanpoint( 
             double&    r_tan,
             double&    lat_tan,
             double&    lon_tan,
     ConstVectorView    refellipsoid,
       const double&    r,
       const double&    lat,
       const double&    lon,
       const double&    za,
       const double&    aa );
*/

void los2xyz( 
         Numeric&   za, 
         Numeric&   aa, 
   const Numeric&   r1,
   const Numeric&   lat1,    
   const Numeric&   lon1,
   const Numeric&   x1, 
   const Numeric&   y1, 
   const Numeric&   z1, 
   const Numeric&   x2, 
   const Numeric&   y2, 
   const Numeric&   z2 );

void poslos2cart(
              double&   x,
              double&   y,
              double&   z,
              double&   dx,
              double&   dy,
              double&   dz,
        const double&   r,
        const double&   lat,
        const double&   lon,
        const double&   za,
        const double&   aa );

double pos2refell_r(
       const Index&     atmosphere_dim,
       ConstVectorView  refellipsoid,
       ConstVectorView  lat_grid,
       ConstVectorView  rte_pos );

double refell2r(
       ConstVectorView  refellipsoid,
       const double&   lat );

double refell2d(
       ConstVectorView  refellipsoid,
       ConstVectorView  lat_grid,
       const GridPos    gp );

void sph2cart(
            double&   x,
            double&   y,
            double&   z,
      const double&   r,
      const double&   lat,
      const double&   lon );

#endif  // geodetic_h
