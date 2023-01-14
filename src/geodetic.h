/* Copyright (C) 2021
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

/**
    @file    geodetic.h
    @author  Patrick Eriksson <Patrick.Eriksson@chalmers.se>
    @date    2021-07-29 

    @brief   This file contains the definition of functions associated
             with the reference ellipsoid, conversion between coordinate
             systems and similar stuff.
*/

#ifndef geodetic_h
#define geodetic_h

#include "matpack_data.h"


/** Threshold for non-spherical ellipsoid
  
    If the two radii of an ellipsoid differ with less than this value, it is
    treated as spherical for efficiency reasons.
*/
const Numeric ellipsoid_radii_threshold = 1e-3;


/** Size of north and south poles
  
    Latitudes with an absolute value > POLELAT are considered to be on
    the south or north pole. This is needed for definition of azimuth.
*/
const Numeric POLELATZZZ = 90 - 1e-8;   // Rename to POLELAT when other one removed


/** Conversion from ECEF to geocentric coordinates

    @param[out]  pos   Geocentric position (r,lat,lon)
    @param[in]   ecef  ECEF position (x,y,z)

    \author  Patrick Eriksson
    \date    2021-07-29
*/
void ecef2geocentric(VectorView pos,
                     ConstVectorView ecef);


/** Conversion from ECEF to geocentric coordinates including a LOS

    @param[out]  pos    Geocentric position (r,lat,lon)
    @param[out]  los    Geocentric line-of-sight (za,aa)
    @param[in]   ecef   ECEF position (x,y,z)
    @param[in]   decef  Normalised ECEF direction vector (dx,dy,dz)

    \author  Patrick Eriksson
    \date    2021-07-29
*/
void ecef2geocentric_los(VectorView pos,
                         VectorView los,
                         ConstVectorView ecef,
                         ConstVectorView decef);


/** Conversion from ECEF to geodetic coordinates

    @param[out]  pos           Geodetic position (h,lat,lon)
    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   refellipsoid  As the WSV with same name.

    @author  Patrick Eriksson
    @date    2021-07-29
*/
void ecef2geodetic(VectorView pos,
                   ConstVectorView ecef,
                   ConstVectorView refellipsoid);


/** Conversion from ECEF to geodetic coordinates including a LOS

    @param[out]  pos           Geodetic position (h,lat,lon)
    @param[out]  los           Geodetic line-of-sight (za,aa)
    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   refellipsoid  As the WSV with the same name.

    @author  Patrick Eriksson
    @date    2021-07-29
*/
void ecef2geodetic_los(VectorView pos,
                       VectorView los,
                       ConstVectorView ecef,
                       ConstVectorView decef,
                       ConstVectorView refellipsoid);


/** ECEF position at a given distance 

    *ecef* is *ecef0* is moved the distance *l* in the direction
    *specified by decef*. *ecef* and *ecef0* can be the argument.

    @param[in]  ecef   New position
    @param[in]  ecef0  Start ECEF position (x,y,z)
    @param[in]  decef  Normalised ECEF direction vector (dx,dy,dz)
    @param[in]  l      A distance

    @author  Patrick Eriksson
    @date    2021-07-29
*/
void ecef_at_distance(VectorView ecef,
                      ConstVectorView ecef0,
                      ConstVectorView decef,
                      const Numeric l);


/** The distance between two ECEF positions

   @param[in]  ecef1  A first ECEF position (x,y,z)
   @param[in]  ecef2  A second ECEF position (x,y,z)

   @return  The distance

   @author  Patrick Eriksson
   @date    2021-08-08
*/
Numeric ecef_distance(ConstVectorView ecef1,
                      ConstVectorView ecef2);


/** The vector between two ECEF positions

    Calculates the vector from ecef0 to ecef1. 

    @param[out]  ecef   The difference vector (x,y,z)
    @param[in]   ecef0  The reference ECEF position
    @param[in]   ecef1  Target ECEF position

    @author  Patrick Eriksson
    @date    2021-08-10
*/
void ecef_vector_distance(VectorView ecef,
                          ConstVectorView ecef0,
                          ConstVectorView ecef1);


/** Converts ENU unit vector vector to local zenith and azimuth

    This function and the sister function los2enu handles transformation of
    line-of-sights, from and to ENU (east-north-up). The ENU vector is
    normalised to have length 1.

    @param[out]  los  Geodetic line-of-sight (za,aa)
    @param[in]   enu  ENU vector (de,dn,du)

   @author  Patrick Eriksson
   @date    2020-09-17
*/
void enu2los(VectorView los,
             ConstVectorView enu);


/** Conversion from geocentric to ECEF coordinates

   @param[out]  ecef  ECEF position (x,y,z)
   @param[in]   pos   Geocentric position (r,lat,lon)

   @author  Patrick Eriksson
   @date    2021-07-29
*/
void geocentric2ecef(VectorView ecef,
                     ConstVectorView pos);


/** Conversion from geocentric to ECEF coordinates including a LOS

    @param[out]  ecef   ECEF position (x,y,z)
    @param[out]  decef  Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   pos    Geocentric position (r,lat,lon)
    @param[in]   los    Geocentric line-of-sight (za,aa)

    @author  Patrick Eriksson
    @date    2021-07-29
*/
void geocentric_los2ecef(VectorView ecef,
                         VectorView decef,
                         ConstVectorView pos,
                         ConstVectorView los);


/** Conversion from geodetic to ECEF coordinates.
 
    @param[out]  ecef          ECEF position (x,y,z)
    @param[in]   pos           Geodetic position (h,lat,lon)
    @param[in]   refellipsoid  As the WSV with the same name.
  
    @author Patrick Eriksson
    @date   2020-09-17
*/
void geodetic2ecef(VectorView ecef,
                   ConstVectorView pos,
                   ConstVectorView refellipsoid);


/** Conversion from geodetic to ECEF coordinates including a LOS
 
    @param[out]  ecef          ECEF position (x,y,z)
    @param[out]  decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   pos           Geodetic position (h,lat,lon)
    @param[in]   los           Geodetic line-of-sight (za,aa)
    @param[in]   refellipsoid  As the WSV with the same name.
  
    @author  Patrick Eriksson
    @date    2020-09-17
*/
void geodetic_los2ecef(VectorView ecef,
                       VectorView decef,
                       ConstVectorView pos,
                       ConstVectorView los,
                       ConstVectorView refellipsoid);


/** Calculates the geometrical tangent point, approximately

    The tangent point is defined as the lowest altitude above the reference
    ellipsiod (not the lowest radius) of the propagation path.

    The tangent point can be below the surface. If the zenith angle is abobe 90
    deg, the tangent point is behind the sensor.

    After testing it was found that the algorithm is not fully numerically
    stable and the method is "flagged" as approximative. It was found that the
    calculated tangent point could match zenith angles of 90.3 deg (it should be
    90 deg exactly).

    @param[in]  ecef_tan      ECEF position of the tangent point (x,y,z)
    @param[in]  ecef          ECEF position (x,y,z)
    @param[in]  decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]  refellipsoid  As the WSV with the same name.

    @author  Patrick Eriksson
    @date    2021-08-02
*/
void approx_geometrical_tangent_point(VectorView ecef_tan,
                                      ConstVectorView ecef,
                                      ConstVectorView decef,
                                      ConstVectorView refellipsoid);


/** Finds the distance to the intersection between an ECEF line and an ellipsoid

    A negative distance is returned if there is no intersection. 

    If multiple postive solutions, the smallest distance >= l_min is returned. 

    Note that there is no check if *pos* is below *altitude*, and the solution 
    can match a position on the other side of the planet.

    @param[in]  ecef          ECEF position (x,y,z)
    @param[in]  decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]  refellipsoid  As the WSV with the same name.
    @param[in]  altitude      The target altitude 
    @param[in]  l_min         Minimum distance of interest

    @return  The distance

    @author  Patrick Eriksson
    @date    2021-07-28
*/
Numeric intersection_altitude(ConstVectorView ecef,
                              ConstVectorView decef,
                              ConstVectorView refellipsoid,
                              const Numeric& altitude,
                              const Numeric& l_min = 0);


/** Finds the distance to the intersection between an ECEF line and a latitude

    A negative distance is returned if there is no intersection.

    @param[in]  ecef          ECEF position (x,y,z)
    @param[in]  decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]  pos           Geodetic position (h,lat,lon)
    @param[in]  los           Geodetic line-of-sight (za,aa)
    @param[in]  refellipsoid  As the WSV with the same name.
    @param[in]  lat           The target latitude

    @return  The distance

    @author  Patrick Eriksson
    @date    2021-08-10
*/
Numeric intersection_latitude(ConstVectorView ecef,
                              ConstVectorView decef,
                              ConstVectorView pos,
                              ConstVectorView los,
                              ConstVectorView refellipsoid,
                              const Numeric& lat);


/** Finds the distance to the intersection between an ECEF line and a longitude

    Looks only for local solutions. The azimuth angle in *los* must be in the
    direction of the target longitude. Otherwise it is considered that there is
    no intersection.

    A negative distance is returned if there is no intersection.

    @param[in]  ecef   ECEF position (x,y,z)
    @param[in]  decef  Normalised ECEF direction vector (dx,dy,dz)
    @param[in]  pos    Geodetic position (h,lat,lon)
    @param[in]  los    Geodetic line-of-sight (za,aa)
    @param[in]  lon    The target longitude

    @return  The distance

    @author  Patrick Eriksson
    @date    2021-07-30
*/
Numeric intersection_longitude(ConstVectorView ecef,
                               ConstVectorView decef,
                               ConstVectorView pos,
                               ConstVectorView los,
                               const Numeric& lon);


/** Determines if an ellipsoid can be treated as a sphere

   @param[in] ellipsoid Vector defining ellipsoid           

   @return   True or false

   @author Patrick Eriksson
   @date   2021-07-28
*/
bool is_ellipsoid_spherical(ConstVectorView ellipsoid);


/** Converts local zenith and azimuth angles to ENU unit vector.

    This function and the sister function enu2los handles transformation of
    line-of-sights, from and to ENU (east-north-up). The ENU vector is
    normalised to have length 1.

    @param[in]  enu ENU vector (de,dn,du)
    @param[in]  los  Geodetic line-of-sight (za,aa)

    @author  Patrick Eriksson
    @date    2020-09-17
*/
void los2enu(VectorView enu,
             ConstVectorView los);


/** Reverses a line-of-sight

    The new LOS is the one for going in the beckward direction.

    @param[out]  los_new  The line-of-sight for reversed direction.
    @param[in]   los      A line-of-sight

    @author  Patrick Eriksson 
    @date    2023-01-14
*/
void reverse_los(VectorView los_new,
                 ConstVectorView los);


/** Geodetic position and line-of-sightat a given distance

    @param[out]  pos           Geodetic position (h,lat,lon)
    @param[out]  los           Geodetic LOS (za,aa)
    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   refellipsoid  As the WSV with the same name.
    @param[in]   l             Distance from *ecef*

    @author  Patrick Eriksson
    @date    2021-08-12
*/
void poslos_at_distance(VectorView pos,
                        VectorView los,
                        ConstVectorView ecef,
                        ConstVectorView decef,
                        ConstVectorView refellipsoid,
                        const Numeric l);


/** Geodetic position at a given distance

    @param[out]  pos           Geodetic position (h,lat,lon)
    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   refellipsoid  As the WSV with the same name.
    @param[in]   l             Distance from *ecef*

    @author  Patrick Eriksson
    @date    2021-08-12
*/
void pos_at_distance(VectorView pos,
                     ConstVectorView ecef,
                     ConstVectorView decef,
                     ConstVectorView refellipsoid,
                     const Numeric l);


/** The prime vertical radius

    This radius is the distance to the polar axis from the reference
    ellipsoid surface in the local nadir direction.

    See further https://en.wikipedia.org/wiki/Earth_radius#Prime_vertical

    @param[in]  refellipsoid  As the WSV with the same name.
    @param[in]  lat           Latitude

    @return  Radius.

    @author  Patrick Eriksson
    @date    2023-01-06
*/
Numeric prime_vertical_radius(ConstVectorView refellipsoid,
                              const Numeric& lat);

#endif  // geodetic_h
