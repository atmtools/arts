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
#include "matpack_data.h"

// 2D:

void cart2pol(Numeric& r,
              Numeric& lat,
              const Numeric& x,
              const Numeric& z,
              const Numeric& lat0,
              const Numeric& za0);

void cart2poslos(Numeric& r,
                 Numeric& lat,
                 Numeric& za,
                 const Numeric& x,
                 const Numeric& z,
                 const Numeric& dx,
                 const Numeric& dz,
                 const Numeric& ppc,
                 const Numeric& lat0,
                 const Numeric& za0);

void distance2D(Numeric& l,
                const Numeric& r1,
                const Numeric& lat1,
                const Numeric& r2,
                const Numeric& lat2);

/*
void geomtanpoint2d( 
             Numeric&    r_tan,
             Numeric&    lat_tan,
     ConstVectorView    refellipsoid,
       const Numeric&    r,
       const Numeric&    lat,
       const Numeric&    za );
*/

void line_refellipsoid_intersect(Numeric& l,
                                 const Vector& refellipsoid,
                                 const Numeric& x,
                                 const Numeric& y,
                                 const Numeric& z,
                                 const Numeric& dx,
                                 const Numeric& dy,
                                 const Numeric& dz);

void line_circle_intersect(Numeric& x,
                           Numeric& z,
                           const Numeric& xl,
                           const Numeric& zl,
                           const Numeric& dx,
                           const Numeric& dz,
                           const Numeric& xc,
                           const Numeric& zc,
                           const Numeric& r);

void pol2cart(Numeric& x, Numeric& z, const Numeric& r, const Numeric& lat);

void poslos2cart(Numeric& x,
                 Numeric& z,
                 Numeric& dx,
                 Numeric& dz,
                 const Numeric& r,
                 const Numeric& lat,
                 const Numeric& za);

// 3D:

void cart2poslos(Numeric& r,
                 Numeric& lat,
                 Numeric& lon,
                 Numeric& za,
                 Numeric& aa,
                 const Numeric& x,
                 const Numeric& y,
                 const Numeric& z,
                 const Numeric& dx,
                 const Numeric& dy,
                 const Numeric& dz,
                 const Numeric& ppc,
                 const Numeric& x0,
                 const Numeric& y0,
                 const Numeric& z0,
                 const Numeric& lat0,
                 const Numeric& lon0,
                 const Numeric& za0,
                 const Numeric& aa0);

void cart2poslos_plain(Numeric& r,
                       Numeric& lat,
                       Numeric& lon,
                       Numeric& za,
                       Numeric& aa,
                       const Numeric& x,
                       const Numeric& y,
                       const Numeric& z,
                       const Numeric& dx,
                       const Numeric& dy,
                       const Numeric& dz);

void cart2sph(Numeric& r,
              Numeric& lat,
              Numeric& lon,
              const Numeric& x,
              const Numeric& y,
              const Numeric& z,
              const Numeric& lat0,
              const Numeric& lon0,
              const Numeric& za0,
              const Numeric& aa0);

void cart2sph_plain(Numeric& r,
                    Numeric& lat,
                    Numeric& lon,
                    const Numeric& x,
                    const Numeric& y,
                    const Numeric& z);

void distance3D(Numeric& l,
                const Numeric& r1,
                const Numeric& lat1,
                const Numeric& lon1,
                const Numeric& r2,
                const Numeric& lat2,
                const Numeric& lon2);

void geompath_tanpos_3d(Numeric& r_tan,
                        Numeric& lat_tan,
                        Numeric& lon_tan,
                        Numeric& l_tan,
                        const Numeric& r,
                        const Numeric& lat,
                        const Numeric& lon,
                        const Numeric& za,
                        const Numeric& aa,
                        const Numeric& ppc);

/*
void geomtanpoint( 
             Numeric&    r_tan,
             Numeric&    lat_tan,
             Numeric&    lon_tan,
     ConstVectorView    refellipsoid,
       const Numeric&    r,
       const Numeric&    lat,
       const Numeric&    lon,
       const Numeric&    za,
       const Numeric&    aa );
*/

void latlon_at_aa(Numeric& lat2,
                  Numeric& lon2,
                  const Numeric& lat1,
                  const Numeric& lon1,
                  const Numeric& aa,
                  const Numeric& ddeg);

void line_sphere_intersect(Numeric& x,
                           Numeric& y,
                           Numeric& z,
                           const Numeric& xl,
                           const Numeric& yl,
                           const Numeric& zl,
                           const Numeric& dx,
                           const Numeric& dy,
                           const Numeric& dz,
                           const Numeric& xc,
                           const Numeric& yc,
                           const Numeric& zc,
                           const Numeric& r);

void los2xyz(Numeric& za,
             Numeric& aa,
             const Numeric& r1,
             const Numeric& lat1,
             const Numeric& lon1,
             const Numeric& x1,
             const Numeric& y1,
             const Numeric& z1,
             const Numeric& x2,
             const Numeric& y2,
             const Numeric& z2);

void poslos2cart(Numeric& x,
                 Numeric& y,
                 Numeric& z,
                 Numeric& dx,
                 Numeric& dy,
                 Numeric& dz,
                 const Numeric& r,
                 const Numeric& lat,
                 const Numeric& lon,
                 const Numeric& za,
                 const Numeric& aa);

Numeric pos2refell_r(const Index& atmosphere_dim,
                     ConstVectorView refellipsoid,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstVectorView rte_pos);

Numeric refell2r(ConstVectorView refellipsoid, const Numeric& lat);

Numeric refell2d(ConstVectorView refellipsoid,
                 ConstVectorView lat_grid,
                 const GridPos gp);

Numeric sphdist(const Numeric& lat1,
                const Numeric& lon1,
                const Numeric& lat2,
                const Numeric& lon2);

void sph2cart(Numeric& x,
              Numeric& y,
              Numeric& z,
              const Numeric& r,
              const Numeric& lat,
              const Numeric& lon);

// longitude functions

void lon_shiftgrid(Vector& longrid_out,
                   ConstVectorView longrid_in,
                   const Numeric lon);

void cycle_lat_lon(Numeric& lat, Numeric& lon);

// Functions involving geodetic latitude

void geodetic2cart(Numeric& x,
                   Numeric& y,
                   Numeric& z,
                   const Numeric& h,
                   const Numeric& lat_gd,
                   const Numeric& lon_gd,
                   const Vector& refellipsoid );

void geodeticposlos2cart(Numeric& x,
                         Numeric& y,
                         Numeric& z,
                         Numeric& dx,
                         Numeric& dy,
                         Numeric& dz,
                         const Numeric& h,
                         const Numeric& lat_gd,
                         const Numeric& lon_gd,
                         const Numeric& za,
                         const Numeric& aa,
                         const Vector& refellipsoid );

void cart2geodeticposlos(Numeric& h,
                         Numeric& lat,
                         Numeric& lon,
                         Numeric& za,
                         Numeric& aa,
                         const Numeric& x,
                         const Numeric& y,
                         const Numeric& z,
                         const Numeric& dx,
                         const Numeric& dy,
                         const Numeric& dz,
                         const Vector& refellipsoid );
#endif  // geodetic_h
