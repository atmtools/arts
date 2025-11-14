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

#include <matpack.h>

/** Conversion from ECEF to geocentric coordinates

    @param[in]   ecef  ECEF position (x,y,z)
    @return      pos   Geocentric position (r,lat,lon)

    \author  Patrick Eriksson
    \date    2021-07-29
*/
Vector3 ecef2geocentric(Vector3 ecef);

/** Conversion from ECEF to geocentric coordinates including a LOS

    @param[in]   ecef   ECEF position (x,y,z)
    @param[in]   decef  Normalised ECEF direction vector (dx,dy,dz)
    @return      pos    Geocentric position (r,lat,lon)
    @return      los    Geocentric line-of-sight (za,aa)

    \author  Patrick Eriksson
    \date    2021-07-29
*/
std::pair<Vector3, Vector2> ecef2geocentric_los(Vector3 ecef, Vector3 decef);

/** Conversion from ECEF to geodetic coordinates

    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   refellipsoid  As the WSV with same name.
    @return      pos           Geodetic position (h,lat,lon)

    @author  Patrick Eriksson
    @date    2021-07-29
*/
Vector3 ecef2geodetic(Vector3 ecef, Vector2 refellipsoid);

/** Conversion from ECEF to geodetic coordinates including a LOS

    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   refellipsoid  As the WSV with the same name.
    @return      pos           Geodetic position (h,lat,lon)
    @return      los           Geodetic line-of-sight (za,aa)

    @author  Patrick Eriksson
    @date    2021-07-29
*/
std::pair<Vector3, Vector2> ecef2geodetic_los(Vector3 ecef,
                                              Vector3 decef,
                                              Vector2 refellipsoid);

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
Vector3 ecef_at_distance(Vector3 ecef0, Vector3 decef, Numeric l);

/** The distance between two ECEF positions

   @param[in]  ecef1  A first ECEF position (x,y,z)
   @param[in]  ecef2  A second ECEF position (x,y,z)

   @return  The distance

   @author  Patrick Eriksson
   @date    2021-08-08
*/
Numeric ecef_distance(Vector3 ecef1, Vector3 ecef2);

/** The vector between two ECEF positions

    Calculates the vector from ecef0 to ecef1. 

    @param[in]   ecef0  The reference ECEF position
    @param[in]   ecef1  Target ECEF position
    @return      ecef   The difference vector (x,y,z)

    @author  Patrick Eriksson
    @date    2021-08-10
*/
Vector3 ecef_vector_distance(Vector3 ecef0, Vector3 ecef1);

/** Converts ENU unit vector vector to local zenith and azimuth

    This function and the sister function los2enu handles transformation of
    line-of-sights, from and to ENU (east-north-up). The ENU vector is
    normalised to have length 1.

    @param[in]   enu  ENU vector (de,dn,du)
    @return      los  Geodetic line-of-sight (za,aa)

   @author  Patrick Eriksson
   @date    2020-09-17
*/
Vector2 enu2los(Vector3 enu);

/** Conversion from geocentric to ECEF coordinates

   @param[in]   pos   Geocentric position (r,lat,lon)
   @return      ecef  ECEF position (x,y,z)

   @author  Patrick Eriksson
   @date    2021-07-29
*/
Vector3 geocentric2ecef(Vector3 pos);

/** Conversion from geocentric to ECEF coordinates including a LOS

    @param[in]   pos    Geocentric position (r,lat,lon)
    @param[in]   los    Geocentric line-of-sight (za,aa)
    @return      ecef   ECEF position (x,y,z)
    @return      decef  Normalised ECEF direction vector (dx,dy,dz)

    @author  Patrick Eriksson
    @date    2021-07-29
*/
std::pair<Vector3, Vector3> geocentric_los2ecef(Vector3 pos, Vector2 los);

/** Conversion from geodetic to ECEF coordinates.
 
    @param[in]   pos           Geodetic position (h,lat,lon)
    @param[in]   refellipsoid  As the WSV with the same name.
    @return      ecef          ECEF position (x,y,z)
  
    @author Patrick Eriksson
    @date   2020-09-17
*/
Vector3 geodetic2ecef(Vector3 pos, Vector2 refellipsoid);

/** Conversion from geodetic to ECEF coordinates including a LOS
 
    @param[in]   pos           Geodetic position (h,lat,lon)
    @param[in]   los           Geodetic line-of-sight (za,aa)
    @param[in]   refellipsoid  As the WSV with the same name.
    @return      ecef          ECEF position (x,y,z)
    @return      decef         Normalised ECEF direction vector (dx,dy,dz)

    @author  Patrick Eriksson
    @date    2020-09-17
*/
std::pair<Vector3, Vector3> geodetic_los2ecef(Vector3 pos,
                                              Vector2 los,
                                              Vector2 refellipsoid);

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
Vector3 approx_geometrical_tangent_point(Vector3 ecef,
                                         Vector3 decef,
                                         Vector2 refellipsoid);

/** Finds the distance to the intersection between an ECEF line and an ellipsoid

    A negative distance is returned if there is no intersection. 

    If multiple postive solutions, the smallest distance >= l_min is returned. 

    Note that there is no check if *pos* is below *alt*, and the solution 
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
Numeric intersection_altitude(Vector3 ecef,
                              Vector3 decef,
                              Vector2 refellipsoid,
                              Numeric altitude,
                              Numeric l_min = 0);

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
Numeric intersection_latitude(Vector3 ecef,
                              Vector3 decef,
                              Vector3 pos,
                              Vector2 los,
                              Vector2 refellipsoid,
                              Numeric lat);

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
Numeric intersection_longitude(
    Vector3 ecef, Vector3 decef, Vector3 pos, Vector2 los, Numeric lon);

/** Determines if an ellipsoid can be treated as a sphere

   @param[in] ellipsoid Vector defining ellipsoid           

   @return   True or false

   @author Patrick Eriksson
   @date   2021-07-28
*/
bool is_ellipsoid_spherical(Vector2 ellipsoid);

/** Converts local zenith and azimuth angles to ENU unit vector.

    This function and the sister function enu2los handles transformation of
    line-of-sights, from and to ENU (east-north-up). The ENU vector is
    normalised to have length 1.

    @param[in]  enu ENU vector (de,dn,du)
    @param[in]  los  Geodetic line-of-sight (za,aa)

    @author  Patrick Eriksson
    @date    2020-09-17
*/
Vector3 los2enu(Vector2 los);

/** Reverses a line-of-sight

    The new LOS is the one for going in the beckward direction.

    @param[in]   los      A line-of-sight
    @return      los_new  The line-of-sight for reversed direction.

    @author  Patrick Eriksson 
    @date    2023-01-14
*/
Vector2 reverse_los(Vector2 los);

/** Geodetic position and line-of-sightat a given distance

    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   refellipsoid  As the WSV with the same name.
    @param[in]   l             Distance from *ecef*
    @return      pos           Geodetic position (h,lat,lon)
    @return      los           Geodetic LOS (za,aa)

    @author  Patrick Eriksson
    @date    2021-08-12
*/
std::pair<Vector3, Vector2> poslos_at_distance(Vector3 ecef,
                                               Vector3 decef,
                                               Vector2 refellipsoid,
                                               Numeric l);

/** Geodetic position at a given distance

    @param[in]   ecef          ECEF position (x,y,z)
    @param[in]   decef         Normalised ECEF direction vector (dx,dy,dz)
    @param[in]   refellipsoid  As the WSV with the same name.
    @param[in]   l             Distance from *ecef*
    @return      pos           Geodetic position (h,lat,lon)

    @author  Patrick Eriksson
    @date    2021-08-12
*/
Vector3 pos_at_distance(Vector3 ecef,
                        Vector3 decef,
                        Vector2 refellipsoid,
                        Numeric l);

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
Numeric prime_vertical_radius(Vector2 refellipsoid, Numeric lat);

/** Converts geodetic coordinates to geocentric coordinates

    @param[in]  pos  Geodetic position (h,lat,lon)
    @param[in]  ell  Ellipsoid radii (a,b)

    @return  Geocentric position (r,lat,lon)

    @author  Patrick Eriksson
    @date    2021-07-29
*/
Vector3 geodetic2geocentric(Vector3 pos, Vector2 ell);

#endif  // geodetic_h
