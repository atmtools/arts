/* Copyright (C)  Patrick Eriksson <patrick@rss.chalmers.se>

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
  \file   ppath.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-02
  
  \brief  Functions releated to calculation of propagation paths.
  
  Functions to determine propagation paths for different atmospheric
  dimensionalities, with and without refraction.

  The term propagation path is here shortened to ppath.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "mystring.h"
#include "logic.h"
#include "poly_roots.h"
#include "ppath.h"
#include "refraction.h"
#include "special_interp.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



/*===========================================================================
  === Common precision variables
  ===========================================================================*/

// This variable defines the maximum allowed error tolerance for radius.
// The variable is, for example, used to check that a given a radius is
// consistent with the specified grid cell.
//
#ifdef USE_DOUBLE
const double   RTOL = 1e-6;
#else
const double   RTOL = 1;
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



/*===========================================================================
  === Functions related to geometrical propagation paths
  ===========================================================================*/

//! geometrical_ppc
/*! 
   Calculates the propagation path constant for pure geometrical calculations.

   Both positive and negative zenith angles are handled.

   \return         Path constant.
   \param   r      Radius of the sensor position.
   \param   za     Zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double geometrical_ppc( const double& r, const double& za )
{
  assert( r > 0 );
  assert( abs(za) <= 180 );

  return r * sin( DEG2RAD * abs(za) );
}



//! geompath_za_at_r
/*! 
   Calculates the zenith angle for a given radius along a geometrical 
   propagation path.

   For downlooking cases, the two points must be on the same side of 
   the tangent point.

   Both positive and negative zenith angles are handled.

   \return         Zenith angle at the point of interest.
   \param   ppc    Propagation path constant.
   \param   a_za   A zenith angle along the path on the same side of the 
                   tangent point as the point of interest.  
   \param   r      Radius of the point of interest.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double geompath_za_at_r(
       const double&   ppc,
       const double&   a_za,
       const double&   r )
{
  assert( ppc >= 0 );
  assert( abs(a_za) <= 180 );
  assert( r >= ppc );

  double za = RAD2DEG * asin( ppc / r );
  if( abs(a_za) > 90 )
    { za = 180 - za; }
  if( a_za < 0 )
    { za = -za; }
  return za;
}



//! geompath_r_at_za
/*! 
   Calculates the zenith angle for a given radius along a geometrical 
   propagation path.

   Both positive and negative zenith angles are handled.

   \return         Radius at the point of interest.
   \param   ppc    Propagation path constant.
   \param   za     Zenith angle at the point of interest.

   \author Patrick Eriksson
   \date   2002-06-05
*/
double geompath_r_at_za(
       const double&   ppc,
       const double&   za )
{
  assert( ppc >= 0 );
  assert( abs(za) <= 180 );

  return ppc / sin( DEG2RAD * abs(za) );
}



//! geompath_lat_at_za
/*!
   Calculates the latitude for a given zenith angle along a geometrical 
   propagation path.

   Positive and negative zenith angles are handled. A positive zenith angle
   means a movement towards higher latitudes.

   \return         The latitude of the second point.
   \param   za0    The zenith angle of the starting point.
   \param   lat0   The latitude of the starting point.
   \param   za     The zenith angle of the second point.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double geompath_lat_at_za(
       const double&   za0,
       const double&   lat0,
       const double&   za )
{
  assert( abs(za0) <= 180 );
  assert( abs(za) <= 180 );
  assert( ( za0 >= 0 && za >= 0 )  ||  ( za0 < 0 && za < 0 ) );

  return lat0 + za0 - za;
}



//! geompath_l_at_r
/*!
   Calculates the length from the tangent point for the given radius.

   The tangent point is either real or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         Length along the path from the tangent point.
   \param   ppc    Propagation path constant.
   \param   r      Radius of the point of concern.

   \author Patrick Eriksson
   \date   2002-05-20
*/
double geompath_l_at_r(
       const double&   ppc,
       const double&   r )
{
  assert( ppc >= 0 );
  assert( r >= ppc - RTOL );
  
  if( r > ppc )
    { return   sqrt( r*r - ppc*ppc ); }
  else
    { return   0; }
}



//! geompath_r_at_l
/*!
   Calculates the radius for a distance from the tangent point.

   The tangent point is either rwal or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         Radius. 
   \param   ppc    Propagation path constant.
   \param   l      Length from the tangent point.

   \author Patrick Eriksson
   \date   2002-05-20
*/
double geompath_r_at_l(
       const double&   ppc,
       const double&   l )
{
  assert( ppc >= 0 );
  assert( l >= 0 );
  
  return sqrt( l*l + ppc*ppc );
}



//! geompath_r_at_lat
/*!
   Calculates the radius for a given latitude.

   \return         Radius at the point of interest.
   \param   ppc    Propagation path constant.
   \param   lat0   Latitude at some other point of the path.
   \param   za0    Zenith angle for the point with latitude lat0.
   \param   lat    Latitude of the point of interest.

   \author Patrick Eriksson
   \date   2002-06-05
*/
double geompath_r_at_lat(
       const double&   ppc,
       const double&   lat0,
       const double&   za0,
       const double&   lat )
{
  assert( ppc >= 0 );
  assert( abs(za0) <= 180 );
  assert( ( za0 >= 0 && lat >= lat0 )  ||  ( za0 <= 0 && lat <= lat0 ) );

  // Zenith angle at the new latitude
  const double za = za0 + lat0 -lat;

  return geompath_r_at_za( ppc, za );
}



//! geompath_from_r1_to_r2
/*!
   Determines radii, latitudes and zenith angles between two points of a 
   propagation path.

   Both start and end point are included in the returned vectors.

   \param   r      Output: Radius of propagation path points.
   \param   lat    Output: Latitude of propagation path points.
   \param   za     Output: Zenith angle of propagation path points.
   \param   lstep  Output: Distance along the path between the points. 
   \param   ppc    Propagation path constant.
   \param   r1     Radius for first point.
   \param   lat1   Latitude for first point.
   \param   za1    Zenith angle for first point.
   \param   r2     Radius for second point.
   \param   lmax   Length criterion for distance between path points.
                   A negative value means no length criterion.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void geompath_from_r1_to_r2( 
             Vector&   r,
             Vector&   lat,
             Vector&   za,
             double&   lstep,
       const double&   ppc,
       const double&   r1,
       const double&   lat1,
       const double&   za1,
       const double&   r2,
       const double&   lmax )
{
  // Calculate length along the path for point 1 and 2.
  const double l1 =  geompath_l_at_r( ppc, r1 );
  const double l2 =  geompath_l_at_r( ppc, r2 );
  
  // Calculate needed number of steps, considering a possible length criterion
  Index n;
  if( lmax > 0 )
    {
      // The absolute value of the length distance is needed here
      n = Index( ceil( abs( l2 - l1 ) / lmax ) );
 
      // We can't accept n=0, which is the case if l1 = l2
      if( n < 1 )
        { n = 1; }
    }
  else
    { n = 1; }

  // Length of path steps (note that lstep here can be negative)
  lstep = ( l2 - l1 ) / n;

  // Allocate vectors and put in point 1
  r.resize(n+1);
  lat.resize(n+1);
  za.resize(n+1);
  r[0]   = r1;
  lat[0] = lat1;
  za[0]  = za1;

  // Loop steps (beside last) and calculate radius and zenith angle
  for( Index i=1; i<n; i++ )
    {
      r[i]   = geompath_r_at_l( ppc, l1 + lstep * i );
      za[i]  = geompath_za_at_r( ppc, za1, r[i] );
    }

  // For maximum accuracy, set last radius to be exactly r2.
  r[n]   = r2;
  za[n]  = geompath_za_at_r( ppc, za1, r2 );  // Don't use r[n] here

  // Ensure that zenith and nadir observations keep their zenith angle
  if( abs(za1) < ANGTOL  ||  abs(za1) > 180-ANGTOL )
    { za = za1; }

  // Calculate latitudes
  for( Index i=1; i<=n; i++ )
    { lat[i] = geompath_lat_at_za( za1, lat1, za[i] ); }

  // Take absolute value of lstep
  lstep = abs( lstep );
}



//! za_geom2other_point
/*!
   Calculates the zenith angle for the geometrical propagation path between
   two specified points.

   The returned zenith angle is valid for point 1. That is, the propagation
   path goes from point 1 to point 2.

   \return         Zenith angle.
   \param   r1     Radius for point 1.
   \param   lat1   Latiytude for point 1.
   \param   r2     Radius for point 2.
   \param   lat2   Latitude for point 2.

   \author Patrick Eriksson
   \date   2002-07-03
*/
double za_geom2other_point(
       const double&    r1,
       const double&   lat1,
       const double&    r2,
       const double&   lat2 )
{
  if( lat2 == lat1 )
    {
      if( r1 <= r2 )
        { return 0; }
      else
        { return 180; }
    }
  else
    {
      // Absolute value of the latitude difference
      const double dlat = abs( lat2 - lat1 );

      // The zenith angle is obtained by a combination of the lawes of sine
      // and cosine.
      double za = dlat + RAD2DEG * asin( r1 * sin( DEG2RAD * dlat ) / 
                 sqrt( r1*r1 + r2*r2 - 2 * r1 * r2 * cos( DEG2RAD * dlat ) ) );

      // Consider viewing direction
      if( lat2 < lat1 )
        { za = -za; }

      return za;
    }
}





/*===========================================================================
  === Basic (true) 3D function
  ===========================================================================*/

//! sph2cart
/*! 
   Conversion from spherical to cartesian coordinates.

   The cartesian coordinate system is defined such as the x-axis goes along
   lat=0 and lon=0, y-axis goes along lat=90, and the z-axis goes along lat=0
   and lon=90. 

   \param   x     Out: x position.
   \param   y     Out: y position.
   \param   z     Out: z position.
   \param   r     Radius.
   \param   lat   Latitude.
   \param   lon   Longitude.

   \author Patrick Eriksson
   \date   2002-12-17
*/
void sph2cart(
            double&   x,
            double&   y,
            double&   z,
      const double&   r,
      const double&   lat,
      const double&   lon )
{
  assert( r > 0 );
  assert( abs( lat ) <= 90 );
  assert( abs( lon ) <= 360 );

  const double   latrad = DEG2RAD * lat;
  const double   lonrad = DEG2RAD * lon;

  x = r * cos( latrad );   // Common term for x and z
  z = x * sin( lonrad );
  x = x * cos( lonrad );
  y = r * sin( latrad );
}



//! cart2sph
/*! 
   The inverse of *sph2cart*.

   \param   r     Out: Radius of observation position.
   \param   lat   Out: Latitude of observation position.
   \param   lon   Out: Longitude of observation position.
   \param   x     x-coordinate of observation position.
   \param   y     y-coordinate of observation position.
   \param   z     z-coordinate of observation position.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void cart2sph(
             double&    r,
             double&    lat,
             double&    lon,
       const double&    x,
       const double&    y,
       const double&    z )
{
  r   = sqrt( x*x + y*y + z*z );
  lat = RAD2DEG * asin( y / r );
  lon = RAD2DEG * atan2( z, x ); 
}



//! poslos2cart
/*! 
   Conversion from position and LOS to cartesian coordinates

   A position (in geographical coordinates) and LOS are converted to a
   cartesian position and a viewing vector. The viewing direction is the
   the vector [dx,dy,dz]. This vector is normalised (it has length 1).

   See the user guide for definition on the zenith and azimuth angles.

   \param   x     Out: x-coordinate of observation position.
   \param   y     Out: y-coordinate of observation position.
   \param   z     Out: z-coordinate of observation position.
   \param   dx    Out: x-part of LOS unit vector.
   \param   dy    Out: y-part of LOS unit vector.
   \param   dz    Out: z-part of LOS unit vector.
   \param   r     Radius of observation position.
   \param   lat   Latitude of observation position.
   \param   lon   Longitude of observation position.
   \param   za    LOS zenith angle at observation position.
   \param   aa    LOS azimuth angle at observation position.

   \author Patrick Eriksson
   \date   2002-12-30
*/
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
       const double&   aa )
{
  // lat = +-90 
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  if( abs( lat ) > POLELAT )
    {
      const double   s = sign( lat );

      x = 0;
      z = 0;
      y = s * r;

      dy = s * cos( DEG2RAD * za );
      dx = sin( DEG2RAD * za );
      dz = dx * sin( DEG2RAD * aa );
      dx = dx * cos( DEG2RAD * aa );
    }

  else
    {
      const double   latrad = DEG2RAD * lat;
      const double   lonrad = DEG2RAD * lon;
      const double   zarad  = DEG2RAD * za;
      const double   aarad  = DEG2RAD * aa;

      sph2cart( x, y, z, r, lat, lon );

      const double   coslat = cos( latrad );
      const double   sinlat = sin( latrad );
      const double   coslon = cos( lonrad );
      const double   sinlon = sin( lonrad );
      const double   cosza  = cos( zarad );
      const double   sinza  = sin( zarad );
      const double   cosaa  = cos( aarad );
      const double   sinaa  = sin( aarad );

      const double   dr   = cosza;
      const double   dlat = sinza * cosaa;         // r-terms cancel out below
      const double   dlon = sinza * sinaa / coslat; 

      dx = coslat*coslon * dr - sinlat*coslon * dlat - coslat*sinlon * dlon;
      dy =        sinlat * dr +        coslat * dlat;
      dz = coslat*sinlon * dr - sinlat*sinlon * dlat + coslat*coslon * dlon;
    }
}



//! cart2poslos
/*! 
   The inverse of *poslos2cart*.

   The azimuth angle is set to: <br> 
      0 when the zenith angle is 0 or 180.
      atan2(dz,dx) at the poles (lat = +- 90).

   The longitude is set to 0 at the poles (lat = +- 90).

   \param   r     Out: Radius of observation position.
   \param   lat   Out: Latitude of observation position.
   \param   lon   Out: Longitude of observation position.
   \param   za    Out: LOS zenith angle at observation position.
   \param   aa    Out: LOS azimuth angle at observation position.
   \param   x     x-coordinate of observation position.
   \param   y     y-coordinate of observation position.
   \param   z     z-coordinate of observation position.
   \param   dx    x-part of LOS unit vector.
   \param   dy    y-part of LOS unit vector.
   \param   dz    z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2002-12-30
*/
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
       const double&   dz )
{
  // Assert that LOS vector is normalised
  assert( abs( sqrt( dx*dx + dy*dy + dz*dz ) - 1 ) < 1e-6 );

  // Spherical coordinates
  cart2sph( r, lat, lon, x, y, z );

  // Spherical derivatives
  const double   coslat = cos( DEG2RAD * lat );
  const double   sinlat = sin( DEG2RAD * lat );
  const double   coslon = cos( DEG2RAD * lon );
  const double   sinlon = sin( DEG2RAD * lon );
  const double   dr   = coslat*coslon*dx + sinlat*dy + coslat*sinlon*dz;
  const double   dlat = -sinlat*coslon/r*dx + coslat/r*dy -sinlat*sinlon/r*dz;
  const double   dlon = -sinlon/coslat/r*dx + coslon/coslat/r*dz;

  // LOS angles
  za = RAD2DEG * acos( dr );
  //
  if( za < ANGTOL  ||  za > 180-ANGTOL  )
    { aa = 0; }

  else if( abs( lat ) <= POLELAT )
    {
      aa = RAD2DEG * acos( r * dlat / sin( DEG2RAD * za ) );

      if( isnan( aa ) )
        {
          if( dlat >= 0 )
            { aa = 0; }
          else
            { aa = 180; }
        }
      else if( dlon < 0 )
        { aa = -aa; }
    }

  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  else
    { aa = RAD2DEG * atan2( dz, dx ); }
}



//! resolve_lon
/*! 
   Resolves which longitude angle that shall be used.

   Longitudes are allowed to vary between -360 and 360 degress, while the
   inverse trigonomtric functions returns values between -180 and 180.
   This function determines if the longitude shall be shifted -360 or
   +360 to fit the longitudes set by the user.
   
   The argument *lon* as input is a value calculated by some inverse
   trigonometric function. The arguments *lon5* and *lon6* are the
   lower and upper limit for the probable range for *lon*. The longitude
   *lon* will be shifted with -360 or +360 degrees if such a shift better
   fit *lon5* and *lon6*. No error is given if it is not possible to
   obtain a value between *lon5* and *lon6*. 

   The function exists both in a float and double version to avoid
   unnecessary copying of data for Numeric=float.

   \param   lon    In/Out: Longitude, possible shifted when returned.
   \param   lon5   Lower limit of probable range for lon.
   \param   lon6   Upper limit of probable range for lon

   \author Patrick Eriksson
   \date   2003-01-05
*/
void resolve_lon(
              double&   lon,
        const double&   lon5,
        const double&   lon6 )
{
  assert( lon6 >= lon5 );
  assert( ( lon6 - lon5 ) < 360 );

  if( lon < lon5  || lon > lon6 )
    {
      const double   meanlon = ( lon5 + lon6 ) / 2;
      const double   diff0   = abs( meanlon - lon );

      if( abs( lon + 360 - meanlon ) < diff0 )
        { lon += 360; }
      else if( abs( lon - 360 - meanlon ) < diff0 )
        { lon -= 360; }
    }
}
#ifndef USE_DOUBLE
void resolve_lon(
              float&    lon,
        const double&   lon5,
        const double&   lon6 )
{
  assert( lon6 >= lon5 );
  assert( ( lon6 - lon5 ) < 360 );

  if( lon < lon5  || lon > lon6 )
    {
      const double   meanlon = ( lon5 + lon6 ) / 2;
      const double   diff0   = abs( meanlon - lon );

      if( abs( lon + 360 - meanlon ) < diff0 )
        { lon += 360; }
      else if( abs( lon - 360 - meanlon ) < diff0 )
        { lon -= 360; }
    }
}
#endif



//! geompath_tanpos_3d
/*! 
   Position of the tangent point for 3D cases.

   The zenith angle must be >= 90.

   \param   r_tan       Out: Radius of tangent point.
   \param   lat_tan     Out: Latitude of tangent point.
   \param   lon_tan     Out: Longitude of tangent point.
   \param   l_tan       Out: Distance along path to the tangent point.
   \param   r           Radius of observation position.
   \param   lat         Latitude of observation position.
   \param   lon         Longitude of observation position.
   \param   za          LOS zenith angle at observation position.
   \param   aa          LOS azimuth angle at observation position.
   \param   ppc         FIXME: Add documentation.

   \author Patrick Eriksson
   \date   2002-12-31
*/
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
       const double&    ppc )
{
  assert( za >= 90 );
  assert( r >= ppc );

  double   x, y, z, dx, dy, dz; 

  poslos2cart( x, y, z, dx, dy, dz, r, lat, lon, za, aa );

  l_tan = sqrt( r*r - ppc*ppc );

  cart2sph( r_tan, lat_tan, lon_tan, x+dx*l_tan, y+dy*l_tan, z+dz*l_tan );
}





/*===========================================================================
  === Functions related to slope and tilt of the ground and pressure surfaces
  ===========================================================================*/

//! psurface_slope_2d
/*!
   Calculates the radial slope of the ground or a pressure surface for 2D.

   The radial slope is here the derivative of the radius with respect to the
   latitude. The unit is accordingly m/degree.

   Note that the radius is defined to change linearly between grid points,
   and the slope is constant between to points of the latitude grid.

   Note also that the slope is always calculated with respect to increasing
   latitudes, independently of the zenith angle. The zenith angle is
   only used to determine which grid range that is of interest when the
   position is exactly on top of a grid point. 

   \return              The radial slope [m/degree]
   \param   lat_grid    The latitude grid.
   \param   r_geoid     Radius of the geoid for the latitude dimension.
   \param   z_surf      Geometrical altitude of the ground, or the pressure
                        surface of interest, for the latitide dimension
   \param   gp          Latitude grid position for the position of interest
   \param   za          LOS zenith angle.

   \author Patrick Eriksson
   \date   2002-06-03
*/
double psurface_slope_2d(
        ConstVectorView   lat_grid,           
        ConstVectorView   r_geoid,
        ConstVectorView   z_surf,
        const GridPos&    gp,
        const double&     za )
{
  Index i1 = gridpos2gridrange( gp, abs( za ) >= 0 );
  const double r1 = r_geoid[i1] + z_surf[i1];
  const double r2 = r_geoid[i1+1] + z_surf[i1+1];
  return ( r2 - r1 ) / ( lat_grid[i1+1] - lat_grid[i1] );
}



//! psurface_slope_2d
/*!
   Calculates the radial slope of the ground or a pressure surface for 2D.

   This function returns the same quantity as the function above, but takes
   the radius and latitude at two points of the pressure surface, instead
   of vector input. That is, for this function the interesting latitude range
   is known when calling the function.

   \return         The radial slope [m/degree]
   \param   lat1   A latitude.
   \param   lat2   Another latitude.
   \param   r1     Radius at *lat1*.
   \param   r2     Radius at *lat2*.

   \author Patrick Eriksson
   \date   2002-12-21
*/
double psurface_slope_2d(
        const double&   lat1,
        const double&   lat2,
        const double&   r1,
        const double&   r2 )
{
  return   ( r2 - r1 ) / ( lat2 -lat1 );
}



//! rsurf_at_latlon
/*!
   Determines the radius of a pressure surface or the ground given the
   radius at the corners of a 3D grid cell.

   \return         Radius at the given latitude and longitude.
   \param   lat1   Lower latitude of grid cell.
   \param   lat3   Upper latitude of grid cell.
   \param   lon5   Lower longitude of grid cell.
   \param   lon6   Upper longitude of grid cell.
   \param   r15    Radius at crossing of *lat1* and *lon5*.
   \param   r35    Radius at crossing of *lat3* and *lon5*.
   \param   r36    Radius at crossing of *lat3* and *lon6*.
   \param   r16    Radius at crossing of *lat1* and *lon6*.
   \param   lat    Latitude for which radius shall be determined.
   \param   lon    Longitude for which radius shall be determined.

   \author Patrick Eriksson
   \date   2002-12-30
*/
double rsurf_at_latlon(
       const double&   lat1,
       const double&   lat3,
       const double&   lon5,
       const double&   lon6,
       const double&   r15,
       const double&   r35,
       const double&   r36,
       const double&   r16,
       const double&   lat,
       const double&   lon )
{
  // We can't have any assert of *lat* and *lon* here as we can go outside
  // the ranges when called from *psurface_slope_3d*.

  if( lat == lat1 )
    { return   r15 + ( lon - lon5 ) * ( r16 - r15 ) / ( lon6 -lon5 ); }
  else if( lat == lat3 )
    { return   r35 + ( lon - lon5 ) * ( r36 - r35 ) / ( lon6 -lon5 ); }
  else if( lon == lon5 )
    { return   r15 + ( lat - lat1 ) * ( r35 - r15 ) / ( lat3 -lat1 ); }
  else if( lon == lon6 )
    { return   r16 + ( lat - lat1 ) * ( r36 - r16 ) / ( lat3 -lat1 ); }
  else
    {
      const double   fdlat = ( lat - lat1 ) / ( lat3 - lat1 );
      const double   fdlon = ( lon - lon5 ) / ( lon6 - lon5 );
      return   (1-fdlat)*(1-fdlon)*r15 + fdlat*(1-fdlon)*r35 + 
                                         (1-fdlat)*fdlon*r16 + fdlat*fdlon*r36;
    }
}



//! psurface_slope_3d
/*!
   Calculates the local radial slope of the ground or a pressure surface 
   for 3D.

   The function works basically as the non-vector version of
   *psurface_slope_2d*, but the position and viewing direction must
   here be specicified as the slope varies inside the cell grid, in
   constrast to a 2D latitude grid range.

   See further the other version of the function below.

   \return         The radial slope [m/degree]
   \param   lat1   Lower latitude of grid cell.
   \param   lat3   Upper latitude of grid cell.
   \param   lon5   Lower longitude of grid cell.
   \param   lon6   Upper longitude of grid cell.
   \param   r15    Radius at crossing of *lat1* and *lon5*.
   \param   r35    Radius at crossing of *lat3* and *lon5*.
   \param   r36    Radius at crossing of *lat3* and *lon6*.
   \param   r16    Radius at crossing of *lat1* and *lon6*.
   \param   lat    Latitude for which slope shall be determined.
   \param   lon    Longitude for which slope shall be determined.
   \param   aa     Azimuth angle for which slope shall be determined.

   \author Patrick Eriksson
   \date   2002-12-30
*/
double psurface_slope_3d(
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
        const double&   aa )
{
  // Size of test angular distance. Unit is degrees.
  const double   dang = 1e-4;

  // Radius at point of interest
  const double   r0 = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                                r15, r35, r36, r16, lat, lon );

  // Convert position and an imaginary LOS to cartesian coordinates
  // double is hard-coded to match declaration of poslos2cart and for
  // maximum precision  
  double   x, y, z, dx, dy, dz;
  poslos2cart( x, y, z, dx, dy, dz, r0, lat, lon, 90, aa );

  // Calculate the distance corresponding to *dang*
  const double   l = r0 * DEG2RAD * dang;

  // Calculate radius of the pressure surface at the lat/lon, the distance
  // *l* away.
  double   r2, lat2, lon2, za2, aa2;
  cart2poslos( r2, lat2, lon2, za2, aa2, x+dx*l, y+dy*l, z+dz*l, dx, dy, dz );
  r2 = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15, r35, r36, r16, lat2,lon2);
  
  // Return slope
  return   ( r2 - r0 ) / dang;
}



//! psurface_slope_3d
/*!
   Calculates the radial slope of the ground or a pressure surface for 3D.

   The radial slope is here the derivative of the radius with respect
   to an angular change (in degrees) along the great circle along the
   given azimuth angle. That is, how much the radius would change for a
   movement of r*pi/180 in the given azimuth angle (if the
   slope where constant along the distance). The unit is m/degree.

   For a point exactly on a grid value it is not clear if it is the
   range below or above that is of interest. The azimuth angle is used
   to resolve such cases.

   This function is in practice another way to call the non-vector version
   of the function above.

   \return              The radial slope [m/degree]
   \param   lat_grid    The latitude grid.
   \param   lon_grid    The longitude grid.
   \param   r_geoid     As the WSV with the same name.
   \param   z_surf      Geometrical altitude of the ground, or the pressure
                        surface of interest.
   \param   gp_lat      Latitude grid position for the position of interest.
   \param   gp_lon      Longitude grid position for the position of interest.
   \param   aa          Azimuth angle.

   \author Patrick Eriksson
   \date   2002-06-03
*/
double psurface_slope_3d(
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,  
        ConstMatrixView   r_geoid,
        ConstMatrixView   z_surf,
        const GridPos&    gp_lat,
        const GridPos&    gp_lon,
        const double&    aa )
{
  Index ilat = gridpos2gridrange( gp_lat, abs( aa ) >= 0 );
  Index ilon = gridpos2gridrange( gp_lon, aa >= 0 );

  // Restore latitude and longitude values
  Vector  itw(2);
  double  lat, lon;
  interpweights( itw, gp_lat );
  lat = interp( itw, lat_grid, gp_lat );
  interpweights( itw, gp_lon );
  lon = interp( itw, lon_grid, gp_lon );

  // Extract values that defines the grid cell
  const double   lat1 = lat_grid[ilat];
  const double   lat3 = lat_grid[ilat+1];
  const double   lon5 = lon_grid[ilon];
  const double   lon6 = lon_grid[ilon+1];
  const double   r15  = r_geoid(ilat,ilon) + z_surf(ilat,ilon);
  const double   r35  = r_geoid(ilat+1,ilon) + z_surf(ilat+1,ilon);
  const double   r36  = r_geoid(ilat+1,ilon+1) + z_surf(ilat+1,ilon+1);
  const double   r16  = r_geoid(ilat,ilon+1) + z_surf(ilat,ilon+1);

  return   psurface_slope_3d( lat1, lat3, lon5, lon6, r15, r35, r36, r16, 
                                                                lat, lon, aa );
}



//! psurface_angletilt
/*!
   Calculates the angular tilt of the ground or a pressure surface.

   Note that the tilt value is a local value. The tilt for a constant
   slope value, is different for different radii.

   \return        The angular tilt.
   \param    r    The radius for the surface at the point of interest.
   \param    c    The radial slope, as returned by e.g. psurface_slope_2d.

   \author Patrick Eriksson
   \date   2002-06-03
*/
double psurface_angletilt(
        const double&   r,
        const double&   c )
{
  // The tilt (in radians) is c/r if c is converted to m/radian. So we get
  // conversion RAD2DEG twice
  return   RAD2DEG * RAD2DEG * c / r;
}



//! is_los_downwards
/*!
   Determines if a line-of-sight is downwards compared to the angular tilt
   of the ground or a pressure surface.

   For example, this function can be used to determine if the line-of-sight
   goes into the ground for a starting point exactly on the ground radius.
  
   As the radius of the ground and pressure surfaces varies as a function of
   latitude, it is not clear if a zenith angle of 90 is above or below e.g.
   the ground.
 
   \return         Boolean that is true if LOS is downwards.
   \param   za     Zenith angle of line-of-sight.
   \param   tilt   Angular tilt of the ground or the pressure surface (as
                   returned by psurface_angletilt)

   \author Patrick Eriksson
   \date   2002-06-03
*/
bool is_los_downwards( 
        const double&   za,
        const double&   tilt )
{
  assert( abs(za) <= 180 );

  // Yes, it shall be -tilt in both cases, if you wonder.
  if( za > (90-tilt)  ||  za < (-90-tilt) )
    { return true; }
  else
    { return false; }
}



//! psurface_crossing_2d
/*!
   Calculates the angular distance to a crossing of a 2D pressure surface
   or the ground.

   The function solves the problem mentioned above for a pressure
   surface, or the ground, where the radius changes linearly as a
   function of latitude. No analytical solution to the original
   problem has been found. The problem involves sine and cosine of the
   latitude difference and these functions are replaced with their
   Taylor expansions where the two first terms are kept. This should
   be OK as in practical situations, the latitude difference inside a
   grid cell should not exceed 2 degrees, and the accuracy should be
   sufficient for values up to 3 degrees.

   The problem and its solution is further described in AUG. See the
   chapter on propagation paths.

   Both positive and negative zenith angles are handled. 

   The function only looks for crossings in the forward direction of
   the given zenith angle. This means that if r>r0 and the absolute
   value of the zenith angle is < 90, no crossing will be found (if
   not the slope of the pressure surface happen to be very strong).

   For downlooking cases, only the part down to the tangent point is
   considered.
   
   If the given path point is on the pressure surface (r=r0), the
   solution 0 is rejected.
 
   The latitude difference is set to 999 if no crossing exists.

   The variable names below are the same as in AUG.

   \return         The angular distance to the crossing.
   \param   rp     Radius of a point of the path inside the grid cell
   \param   za     Zenith angle of the path at r.
   \param   r0     Radius of the pressure surface or the ground at the
                   latitude of r.
   \param   c      Linear slope term, as returned by psurface_slope_2d.

   \author Patrick Eriksson
   \date   2002-06-07
*/
double psurface_crossing_2d(
        const double&   rp,
        const double&   za,
        const double&   r0,
              double    c )
{
  assert( abs(za) <= 180 );

  const double no_crossing = 999;

  // Handle the cases of za=0 and za=180. 
  if( abs(za) < ANGTOL )
    {
      if( rp < r0 )
        { return 0; }
      else
        { return no_crossing; }
    }
  if( abs(za) > 180-ANGTOL )
    {
      if( rp > r0 )
        { return 0; }
      else
        { return no_crossing; }
    }

  // If rp>=r0, check if the given LOS goes in the direction towards
  // the pressure surface.  If not, return 999.
  //
  if( rp >= r0   &&  !is_los_downwards( za,  psurface_angletilt( r0, c ) ) )
    { return no_crossing; }

  // The case with c=0 can be handled analytically
  if( c == 0 )
    {
      double   ppc = geometrical_ppc( rp, za );
      if( ( rp < r0  &&  abs(za) <= 90 )  ||  
                                 ( rp > r0  &&  abs(za) > 90  &&  r0 >= ppc ) )
        { return geompath_lat_at_za( za, 0, geompath_za_at_r( ppc, za, r0 ) );}
      else
        { return no_crossing; }
    }


  // Approximative solution:
  else
    {
      // If r0=rp, numerical inaccuracy can give a false solution, very close
      // to 0, that we must throw away.
      double   dmin = 0;
      if( r0 == rp )
        { dmin = 1e-12; }

      // The nadir angle in radians, and cosine and sine of that angle
      const double   beta = DEG2RAD * ( 180 - abs(za) );
      const double   cv = cos( beta );
      const double   sv = sin( beta );

      // Convert slope to m/radian and consider viewing direction
      c *= RAD2DEG;
      if( za < 0 )
        { c = -c; }

      // The vector of polynomial coefficients
      Vector p(5);
      //
      p[0] = ( r0 - rp ) * sv;
      p[1] = r0 * cv + c * sv;
      p[2] = -r0 * sv / 2 + c * cv;
      p[3] = -r0 * cv / 6 - c * sv / 2;
      p[4] = -c * cv / 6;

      // Calculate roots of the polynomial
      Matrix roots(4,2);
      poly_root_solve( roots, p );

      // Find the smallest root with imaginary part = 0, and real part > 0.
      //
      double dlat = no_crossing / RAD2DEG;
      //
      for( Index i=0; i<4; i++ )
        {
          if( roots(i,1) == 0  &&   roots(i,0) > dmin  &&  roots(i,0) < dlat )
            { dlat = roots(i,0); }
        }  

      // Change sign if zenith angle is negative
      if( dlat < no_crossing  &&  za < 0 )
        { dlat = -dlat; }

      return   RAD2DEG * dlat;
    }
}



//! psurface_crossing_3d
/*!
   Calculates the radius of a crossing of a 3D pressure surface or the ground.

   The function solves the problem mentioned above for a pressure
   surface, or the ground, for 3D cases. The problem is solved by making
   calculations for five radii, between the min and max values among
   *r15*, *r35*, *r36* and *r16*. For each test radius, the latitude and 
   longitude for the crossing of the path and the assumed radius are 
   calculated. The test radius is then compared to the radius of the 
   pressure surface, or the ground, for the found latitude and longitude.
   These two radii shall ideally be identical. A radius is selected by
   an interpolation between the test radii.

   The problem and its solution is further described in AUG. See the
   chapter on propagation paths.

   The function only looks for crossings in the forward direction of
   the given zenith angle (neglecting all solutions giving *l* = 0).
 
   The radius *r* is set to -1 if no crossing is found.

   \param   r         Out: Radius of found crossing.
   \param   lat       Out: Latitude of found crossing.
   \param   lon       Out: Longitude of found crossing.
   \param   l         Out: Length along the path to the crossing.
   \param   r_surf    Radius of the surface
   \param   r_start   Radius of observation point.
   \param   lat_start Latitude of start point.
   \param   lon_start Longitude of start point.
   \param   za_start  Zenith angle at observation point.
   \param   x         x-coordinate of observation position.
   \param   y         y-coordinate of observation position.
   \param   z         z-coordinate of observation position.
   \param   dx        x-part of LOS unit vector.
   \param   dy        y-part of LOS unit vector.
   \param   dz        z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void psurface_crossing_3d(
             double&    r,
             double&    lat,
             double&    lon,
             double&    l,
       const double     r_surf,
       const double     r_start,         // No & to avoid problems if these
       const double     lat_start,       // variables are the same as the
       const double     lon_start,       // output ones.
       const double     za_start,
       const double&    x,
       const double&    y,
       const double&    z,
       const double&    dx,
       const double&    dy,
       const double&    dz )
{
  assert( za_start >=   0 );
  assert( za_start <= 180 );

  // Set values for no crossing found
  r   = -1;
  lat = 999;
  lon = 999;
  l   = -1;

  // Handle the cases of za=0 and za=180. 
  if( ( za_start < ANGTOL  &&  r_start < r_surf )  || 
                              ( za_start > 180-ANGTOL  &&  r_start > r_surf ) )
    {
      r   = r_surf;
      lat = lat_start;
      lon = lon_start;
      l   = abs( r_surf - r_start );
    }

  else
    {
      const double   p  = x*dx + y*dy + z*dz;
      const double   pp = p * p;
      const double   q  = x*x + y*y + z*z - r_surf*r_surf;

      const double   l1 = -p + sqrt( pp - q );
      const double   l2 = -p - sqrt( pp - q );

      if( l1 < 0  &&  l2 > 0 )
        { l = l2; }
      else if( l1 > 0  &&  l2 < 0 )
        { l = l1; }
      else if( l1 < l2 )
        { l = l1; }
      else
        { l = l2; }

      if( l > 0 )
        {
          r   = r_surf;
          lat = RAD2DEG * asin( ( y+dy*l ) / r );
          lon = RAD2DEG * atan2( z+dz*l, x+dx*l );
        }

    }
}




/*===========================================================================
  === Path through grid cells and 1D pressure range
  ===========================================================================*/

//! do_gridrange_1d
/*!
   Calculates the geometrical path through a 1D grid range.

   This function works as *do_gridcell_2d*, but is valid for 1D cases.

   The coding of variables and end face is as for *do_gridcell_2d*, with
   the excpetion that end faces 2 and 4 do not exist here.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   tanpoint    Out: Set to 1 if end point is a tangent point.
   \param   r_start0    Radius of start point.
   \param   lat_start   Latitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   ra          Radius of lower pressure surface.
   \param   rb          Radius of upper pressure surface (rb > ra);
   \param   rground     Radius for the ground.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void do_gridrange_1d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   za_v,
              double&   lstep,
              Index&    endface,
              Index&    tanpoint,
        const double&   r_start0,
        const double&   lat_start,
        const double&   za_start,
        const double&   ppc,
        const double&   lmax,
        const double&   ra,
        const double&   rb,
        const double&   rground )
{
  double   r_start = r_start0;

  assert( rb > ra );
  assert( r_start >= ra - RTOL );
  assert( r_start <= rb + RTOL );

  // Shift radius if outside
  if( r_start < ra )
    { r_start = ra; }
  else if( r_start > rb )
    { r_start = rb; }

  // Get end radius of the path step (r_end). If looking downwards, it must 
  // be checked if:
  //    a tangent point is passed
  //    there is an intersection with the ground
  //
  double r_end;
  //
  endface  = 0;
  tanpoint = 0;
  //
  if( za_start <= 90 )
    { 
      r_end   = rb; 
      endface = 3;
    }
  else
    {
      // The tangent radius equals here ppc.

      if( ( ra > rground )  &&  ( ra > ppc ) )
        {
          r_end   = ra;
          endface = 1;
        }
      else if( ppc >= rground )
        {
          r_end    = ppc;
          tanpoint = 1;
          if( ppc != ra )
            { endface = 0; }
        }
      else
        {
          r_end   = rground;
          endface = 7;
        }
    }

  // Calculate basic variables from r_start to r_end.
  //
  geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, lat_start, 
                                                       za_start, r_end, lmax );

  // Force end zenith angle to be exact when we know the correct value
  if( tanpoint )
    { za_v[za_v.nelem()-1] = 90; }
}



//! do_gridcell_2d
/*!
   Calculates the geometrical path through a 2D grid cell.

   The function determines the geometrical path from the given start
   point to the boundary of the grid cell. The face where the path
   exits the grid cell is denoted as the end face. The following
   number coding is used for the variable *endface*: <br>
   1: The face at the lower latitude point. <br>
   2: The face at the lower (geometrically) pressure surface. <br>
   3: The face at the upper latitude point. <br>
   4: The face at the upper (geometrically) pressure surface. <br>
   7: The end point is an intersection with the ground. 

   The corner points are names r[lat][a,b]. For example: r3b.
   The latitudes are numbered to match the end faces. This means that
   the lower latitude has number 1, and the upper number 3. The pressure
   surfaces are named as a and b: <br>
   a: Lower pressure surface (highest pressure). <br>
   b: Upper pressure surface (lowest pressure).

   Path points are included if *lmax*>0 and the distance to the end
   point is > than *lmax*.

   The return vectors (*r_v* etc.) can have any length when handed to
   the function.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   tanpoint    Out: Set to 1 if end point is a tangent point.
   \param   r_start0    Radius of start point.
   \param   lat_start   Latitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   lat1        Latitude of left end face (face 1) of the grid cell.
   \param   lat3        Latitude of right end face (face 3) of the grid cell.
   \param   r1a         Radius of lower-left corner of the grid cell.
   \param   r3a         Radius of lower-right corner of the grid cell.
   \param   r3b         Radius of upper-right corner of the grid cell (r3b>r3a)
   \param   r1b         Radius of upper-left corner of the grid cell (r1b>r1a).
   \param   rground1    Radius for the ground at *lat1*.
   \param   rground3    Radius for the ground at *lat3*.

   \author Patrick Eriksson
   \date   2002-11-28
*/
void do_gridcell_2d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   za_v,
              double&   lstep,
              Index&    endface,
              Index&    tanpoint,
        const double&   r_start0,
        const double&   lat_start,
        const double&   za_start,
        const double&   ppc,
        const double&   lmax,
        const double&   lat1,
        const double&   lat3,
        const double&   r1a,
        const double&   r3a,
        const double&   r3b,
        const double&   r1b,
        const double&   rground1,
        const double&   rground3 )
{
  double   r_start = r_start0;

  // Assert latitudes
  assert( lat_start >= lat1 );
  assert( lat_start <= lat3 );

  // Slopes of pressure surfaces
  const double  c2 = psurface_slope_2d( lat1, lat3, r1a, r3a );
  const double  c4 = psurface_slope_2d( lat1, lat3, r1b, r3b );
  const double  cground = psurface_slope_2d( lat1, lat3, rground1, rground3 );

  // Latitude distance between start point and left grid cell boundary
  const double dlat_left  = lat_start - lat1;

  // Radius of lower and upper pressure surface at *lat_start*.
  const double   rlow = r1a + c2*dlat_left;
  const double   rupp = r1b + c4*dlat_left;

  // Latitude distance to latitude end face in the viewing direction
  double dlat_endface;
  if( za_start >= 0 )
    { dlat_endface = lat3 - lat_start; }
  else
    { dlat_endface = -dlat_left; }

  // Assert radius (some extra tolerance is needed for radius)
  assert( r_start >= rlow - RTOL );
  assert( r_start <= rupp + RTOL );

  // Shift radius if outside
  if( r_start < rlow )
    { r_start = rlow; }
  else if( r_start > rupp )
    { r_start = rupp; }

  // The end point is most easily determined by the latitude difference to 
  // the start point. For some cases there exist several crossings and we want 
  // the one closest in latitude to the start point. The latitude distance 
  // for the crossing shall not exceed dlat2end.
  //
  double   dlat2end = 999;
  double   abs_za_start = abs( za_start );
           endface = 0;


  // --- Lower face 
  //
  // We don't need to consider the face if we are standing
  // on the pressure surface.
  //
  if( r_start > rlow  &&  abs(za_start) > ANGTOL )
    {
      dlat2end = psurface_crossing_2d( r_start, za_start, rlow, c2 );
      endface = 2;  // This variable will be re-set if there was no crossing
    }


  // --- The ground.
  //
  // Check shall be done only if the ground is, at least partly, inside 
  // the grid cell.
  //
  if( rground1 >= r1a  ||  rground3 >= r3a )
    {
      double r_ground = rground1 + cground * dlat_left;
     
      assert( r_start >= r_ground );

      double dlat2ground = psurface_crossing_2d( r_start, za_start, r_ground,
                                                                     cground );
      if( abs(dlat2ground) <= abs(dlat2end) )
        {
          dlat2end = dlat2ground;
          endface  = 7;
        }
    }


  // If dlat2end <= dlat_endface we are ready. Otherwise we have to check
  // remaining cell faces. The same applies after testing upper face.


  // --- Upper face  (pressure surface ip+1).
  //
  if( abs(dlat2end) > abs(dlat_endface)  &&  abs_za_start < 180-ANGTOL )
    {
      // For cases when the tangent point is in-between *r_start* and
      // the pressure surface, 999 is returned. This case will anyhow
      // be handled correctly.

      double dlattry = psurface_crossing_2d( r_start, za_start, rupp, c4 );

      if( abs(dlattry) < abs(dlat2end) )
        {
          dlat2end = dlattry;
          endface  = 4;
        }
    }


  // Left or right end face
  if( abs(dlat2end) > abs(dlat_endface) )
    { 
      dlat2end = dlat_endface; 
      if( za_start >= 0 )
        { endface  = 3; }
      else
        { endface  = 1; }
    }

  assert( endface );

  // Check there is a tangent point inside the grid cell
  if( abs_za_start > 90  &&  ( abs_za_start - abs(dlat2end) ) <= 90 ) 
    { 
      tanpoint = 1;

      // Check if the tangent point is closer than the end point
      if( ( abs_za_start - 90 ) < abs( dlat2end ) )
        { 
          endface  = 0; 
          dlat2end = sign(za_start) * ( abs_za_start - 90 ); 
        }
    }
  else
    { tanpoint = 0; }


  // Calculate radius for end point.
  // To obtain best possible accuracy it is calculated to match found end face,
  // and not based on dlat2end.
  //
  double   r_end = -1;
  //
  if( tanpoint )
    { r_end = geompath_r_at_za( ppc, sign(za_start) * 90 ); }
  else if( endface == 1 )
    { r_end = geompath_r_at_lat( ppc, lat_start, za_start, lat1 ); }
  else if( endface == 2 )
    { r_end = r1a + c2 * ( dlat_left + dlat2end ); }
  else if( endface == 3 )
    { r_end = geompath_r_at_lat( ppc, lat_start, za_start, lat3 ); }
  else if( endface == 4 )
    { r_end = r1b + c4 * ( dlat_left + dlat2end ); }
  else if( endface == 7 )
    { r_end = rground1 + cground * ( dlat_left + dlat2end ); }

  // Fill the return vectors
  //
  geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, lat_start, 
                                                       za_start, r_end, lmax );

  // Force end latitude and zenith angle to be exact as possible
  if( tanpoint )
    {
      if( za_start >= 0 )
        { za_v[za_v.nelem()-1] = 90; }
      else
        { za_v[za_v.nelem()-1] = -90; }
    }
  //
  if( endface == 1 )
    { lat_v[lat_v.nelem()-1] = lat1; }
  else if( endface == 3 )
    { lat_v[lat_v.nelem()-1] = lat3; }
  else
    { lat_v[lat_v.nelem()-1] = lat_start + dlat2end; }
}



//! do_gridcell_3d
/*!
   Calculates the geometrical path through a 2D grid cell.

   The function determines the geometrical path from the given start
   point to the boundary of the grid cell. The face where the path
   exits the grid cell is denoted as the end face. The same number
   coding as in *do_gridcell_2d* is used, where the additional longitude
   end faces are numbered as: <br>
   5: The face at the lower longitude point. <br>
   6: The face at the upper longitude point.

   The corner points are numbered as *do_gridcell_2d*, but 5 or 6 is added
   after the latitude number to indicate the longitude. This means that
   r16a, is the corner at lat1, lon6 and pressure surface a.

   See further *do_gridcell_2d*.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   lon_v       Out: Vector with longitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   aa_v        FIXME: Add documentation.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   tanpoint    Out: Set to 1 if end point is a tangent point.
   \param   r_start0    Radius of start point.
   \param   lat_start0  Latitude of start point.
   \param   lon_start0  Longitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   aa_start    LOS azimuth angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   lat1        Latitude of left end face (face 1) of the grid cell.
   \param   lat3        Latitude of right end face (face 3) of the grid cell.
   \param   lon5        Lower longitude limit of the grid cell.
   \param   lon6        Upper longitude limit of the grid cell.
   \param   r15a        Radius of corner: lower p-surface,*lat1* and *lon5*.
   \param   r35a        Radius of corner: lower p-surface,*lat3* and *lon5*.
   \param   r36a        Radius of corner: lower p-surface,*lat3* and *lon6*.
   \param   r16a        Radius of corner: lower p-surface,*lat1* and *lon6*.
   \param   r15b        Radius of corner: upper p-surface,*lat1* and *lon5*.
   \param   r35b        Radius of corner: upper p-surface,*lat3* and *lon5*.
   \param   r36b        Radius of corner: upper p-surface,*lat3* and *lon6*.
   \param   r16b        Radius of corner: upper p-surface,*lat1* and *lon6*.
   \param   rground15   Radius for the ground at *lat1* and *lon5*.
   \param   rground35   Radius for the ground at *lat3* and *lon5*.
   \param   rground36   Radius for the ground at *lat3* and *lon6*.
   \param   rground16   Radius for the ground at *lat1* and *lon6*.

   \author Patrick Eriksson
   \date   2002-11-28
*/
void do_gridcell_3d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   lon_v,
              Vector&   za_v,
              Vector&   aa_v,
              double&   lstep,
              Index&    endface,
              Index&    tanpoint,
        const double&   r_start0, 
        const double&   lat_start0,
        const double&   lon_start0,
        const double&   za_start,
        const double&   aa_start,
        const double&   ppc,
        const double&   lmax,
        const double&   lat1,
        const double&   lat3,
        const double&   lon5,
        const double&   lon6,
        const double&   r15a,
        const double&   r35a,
        const double&   r36a,
        const double&   r16a,
        const double&   r15b,
        const double&   r35b,
        const double&   r36b,
        const double&   r16b,
        const double&   rground15,
        const double&   rground35,
        const double&   rground36,
        const double&   rground16 )
{
  double   r_start   = r_start0;
  double   lat_start = lat_start0;
  double   lon_start = lon_start0;

  // Assert latitude and longitude
  assert( lat_start >= lat1 - LATLONTOL );
  assert( lat_start <= lat3 + LATLONTOL );
  assert( !( abs( lat_start) < 90  &&  lon_start < lon5 - LATLONTOL ) );
  assert( !( abs( lat_start) < 90  &&  lon_start > lon6 + LATLONTOL ) );

  // Shift latitude and longitude if outside
  if( lat_start < lat1 )
    { lat_start = lat1; }
  else if( lat_start > lat3 )
    { lat_start = lat3; }
  if( lon_start < lon5 )
    { lon_start = lon5; }
  else if( lon_start > lon6 )
    { lon_start = lon6; }

  // Radius of lower and upper pressure surface at the start position
  double   rlow = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                r15a, r35a, r36a, r16a, lat_start, lon_start );
  double   rupp = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                r15b, r35b, r36b, r16b, lat_start, lon_start );

  // Assert radius (some extra tolerance is needed for radius)
  assert( r_start >= rlow - RTOL );
  assert( r_start <= rupp + RTOL );


  // Shift radius if outside
  if( r_start < rlow )
    { r_start = rlow; }
  else if( r_start > rupp )
    { r_start = rupp; }


  // Assert that not standing at the edge looking out.
  // As these asserts are rather costly, they should maybe be removed
  // when everything is carefully tested.
  assert( !( lat_start==lat1 && lat_start!=-90 && abs( aa_start ) > 90 ) );
  assert( !( lat_start==lat3 && lat_start!= 90 && abs( aa_start ) < 90 ) );
  assert( !( lon_start==lon5 && abs(lat_start)!=90 &&  aa_start < 0 ) );
  assert( !( lon_start==lon6 && abs(lat_start)!=90 &&  aa_start > 0 ) );
  //
  if( r_start == rlow )
    {
      const double c = psurface_slope_3d( lat1, lat3, lon5, lon6, 
                      r15a, r35a, r36a, r16a, lat_start, lon_start, aa_start );
      const double tilt = psurface_angletilt( r_start, c );

      assert( !is_los_downwards( za_start, tilt ) );
    }
  else if( r_start == rupp )
    {
      const double c = psurface_slope_3d( lat1, lat3, lon5, lon6, 
                      r15b, r35b, r36b, r16b, lat_start, lon_start, aa_start );
      const double tilt = psurface_angletilt( r_start, c );

      assert( is_los_downwards( za_start, tilt ) );
    }


  // Position and LOS in cartesian coordinates
  double   x, y, z, dx, dy, dz;
  poslos2cart( x, y, z, dx, dy, dz, r_start, lat_start, lon_start, 
                                                          za_start, aa_start );


  // Determine the position of the end point
  //
  endface  = 0;
  tanpoint = 0;
  //
  double   r_end, lat_end, lon_end, l_end;


  // Zenith and nadir looking are handled as special cases

  // Zenith looking
  if( za_start < ANGTOL )
    {
      r_end   = rupp;
      lat_end = lat_start;
      lon_end = lon_start;
      l_end   = rupp - r_start;
      endface  = 4;
    }

  // Nadir looking
  else if( za_start > 180-ANGTOL )
    {
      const double   rground = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                   rground15, rground35, rground36, rground16, 
                                                lat_start, lon_start );
      
      if( rlow > rground )
        {
          r_end  = rlow;
          endface = 2;
        }
      else
        {
          r_end  = rground;
          endface = 7;
        }
      lat_end = lat_start;
      lon_end = lon_start;
      l_end   = r_start - r_end;
    }

  // Check faces in number order
  else
    {
      // Calculate correction terms for the position to compensate for
      // numerical inaccuracy
      //
      double   r_corr, lat_corr, lon_corr;
      //
      cart2sph( r_corr, lat_corr, lon_corr, x, y, z );
      //
      r_corr   -= r_start;
      lat_corr -= lat_start;
      lon_corr -= lon_start;

      // The end point is found by testing different step lengths until the
      // step length has been determined by a precision of *l_acc*.
      //
      // The first step is to found a length that goes outside the grid cell, 
      // to find an upper limit. The lower limit is set to 0. The upper
      // limit is found be testing lengths of 1, 10 100 ... m.
      // The search algorith is bisection, the next length to test is the
      // mean value of the minimum and maximum length limits.

               l_end  = 1;

      double   l_acc  = 1e-3;
      double   l_in   = 0, l_out = l_end;
      bool     ready  = false, startup = true;
      double   abs_aa = abs( aa_start );

      double   l_tan = 99e6;
      if( za_start > 90 )
        { l_tan = sqrt( r_start*r_start - ppc*ppc ); }

      bool   do_ground = false;
      if( rground15 >= r15a  ||  rground35 >= r35a  ||  
                                     rground36 >= r36a  ||  rground16 >= r16a )
        { do_ground = true; }


      while( !ready )
        {
          cart2sph( r_end, lat_end, lon_end, 
                                          x+dx*l_end, y+dy*l_end, z+dz*l_end );
          r_end   -= r_corr;
          lat_end -= lat_corr;
          lon_end -= lon_corr;

          // Special fixes for north-south observations
          if( abs( lat_start ) < 90  &&  ( abs(aa_start) < ANGTOL  ||  
                                           abs(aa_start) > 180-ANGTOL ) )
            { lon_end = lon_start; }
          else
            { resolve_lon( lon_end, lon5, lon6 ); }
          
          // Special fixes for west-east observations
          if( abs( abs_aa - 90 ) < ANGTOL )
            {
              if( lat_start == 0 )
                { lat_end = 0; }
              else if( lat_start > 0  &&  lat_end > lat_start )
                { lat_end = lat_start; }
              else if( lat_start < 0  &&  lat_end < lat_start )
                { lat_end = lat_start; }
            }

          bool   inside = true;

          if( lat_end < lat1 )
            { inside = false;   endface = 1; }
          else if( lat_end > lat3 )
            { inside = false;   endface = 3; }
          else if( lon_end < lon5 )
            { inside = false;   endface = 5; }
          else if( lon_end > lon6 )
            { inside = false;   endface = 6; }

          if( inside ) 
            {
              rlow = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                    r15a, r35a, r36a, r16a, lat_end, lon_end );
              rupp = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                    r15b, r35b, r36b, r16b, lat_end, lon_end );
              
              if( r_end < rlow )
                { inside = false;   endface = 2; }
              else if( r_end > rupp )
                { inside = false;   endface = 4; }

              if( do_ground )
                {
                  const double   r_ground = 
                                  rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                  rground15, rground35, rground36, rground16, 
                                                            lat_end, lon_end );

                  if( r_ground >= rlow  &&  r_end < r_ground )
                    { inside = false;   endface = 7; }
                }
            }              

          if( startup )
            {
              if( inside )
                { 
                  if( l_end < l_tan  &&  l_end*10 > l_tan )
                    { l_end = l_tan; }
                  else
                    { l_end *= 10; }
                }
              else
                { 
                  l_out = l_end;   
                  l_end = ( l_out + l_in ) / 2;
                  startup = false; 
                }
            }
          else
            {
              if( inside )
                { l_in = l_end; }
              else
                { l_out = l_end; }
                            
              if( ( l_out - l_in ) < l_acc )
                { ready = true; }
              else
                { l_end = ( l_out + l_in ) / 2; }
            }
        }

      // Now when we are ready, we remove the correction terms. Otherwise
      // we can end up in an infinite loop if the step length is smaller
      // than the correction.
      r_end   += r_corr;
      lat_end += lat_corr;
      lon_end += lon_corr;

      // Set the relevant coordinate to be consistent with found endface.
      //
      if( endface == 1 )
        { lat_end = lat1; }
      else if( endface == 2 )
        { r_end = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15a, r35a, r36a, 
                                                    r16a, lat_end, lon_end ); }
      else if( endface == 3 )
        { lat_end = lat3; }
      else if( endface == 4 )
        { r_end = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15b, r35b, r36b, 
                                                    r16b, lat_end, lon_end ); }
      else if( endface == 5 )
        { lon_end = lon5; }
      else if( endface == 6 )
        { lon_end = lon6; }
      else if( endface == 7 )
        { r_end = rsurf_at_latlon( lat1, lat3, lon5, lon6, rground15, 
                         rground35, rground36, rground16, lat_end, lon_end ); }

      //--- Tangent point?
      //
      if( za_start > 90 )
        {
          if( l_tan < l_end )
            {
              geompath_tanpos_3d( r_end, lat_end, lon_end, l_end,
                                  r_start, lat_start, lon_start, 
                                  za_start, aa_start, ppc );
              endface  = 0;
              tanpoint = 1;
            }
          else if( l_tan == l_end )
            { tanpoint = 1; }
        }
    }


  //--- Create return vectors
  //
  Index n = 1;
  //
  if( lmax > 0 )
    {
      n = Index( ceil( abs( l_end / lmax ) ) );
      if( n < 1 )
        { n = 1; }
    }
  //
  r_v.resize( n+1 );
  lat_v.resize( n+1 );
  lon_v.resize( n+1 );
  za_v.resize( n+1 );
  aa_v.resize( n+1 );
  //
  r_v[0]   = r_start;
  lat_v[0] = lat_start;
  lon_v[0] = lon_start;
  za_v[0]  = za_start;
  aa_v[0]  = aa_start;
  //
  double  ldouble = l_end / n;
  lstep = ldouble;
  // 
  for( Index j=1; j<=n; j++ )
    {
      const double   l  = ldouble * j;

      // We seperate here float/double to avoid unnecessary copying for
      // Numeric=double. A solution like this is needed to avoid compiler
      // warnings.
#ifdef USE_DOUBLE
      cart2poslos( r_v[j], lat_v[j], lon_v[j], za_v[j], aa_v[j],
                                          x+dx*l, y+dy*l, z+dz*l, dx, dy, dz );
#else
      double   rtmp, lattmp, lontmp, zatmp, aatmp;
      cart2poslos( rtmp, lattmp, lontmp, zatmp, aatmp,
                                          x+dx*l, y+dy*l, z+dz*l, dx, dy, dz );
      r_v[j]   = rtmp;
      lat_v[j] = lattmp;
      lon_v[j] = lontmp;
      za_v[j]  = zatmp;
      aa_v[j]  = aatmp;
#endif
   }


  //--- Set last point especially, which should improve the accuracy
  r_v[n]   = r_end;
  lat_v[n] = lat_end;
  lon_v[n] = lon_end;

  //--- Set last zenith angle to be as accurate as possible
  if( za_start < ANGTOL  ||  za_start > 180-ANGTOL )
    { za_v[n] = za_start; }
  else if( tanpoint )
    { za_v[n] = 90; }
  else
    { za_v[n] = geompath_za_at_r( ppc, za_start, r_v[n] ); }

  //--- Set last azimuth angle to be as accurate as possible for
  //    zenith and nadir observations
  if( abs( lat_start ) < 90  &&  
          ( abs(aa_start) < ANGTOL  ||  abs( aa_start) > 180-ANGTOL ) )
    {  
      aa_v[n]  = aa_start; 
      lon_v[n] = lon_start;
    }

  // Shall lon values be shifted (value 0 and n+1 are already OK)?
  for( Index j=1; j<n; j++ )
    { resolve_lon( lon_v[j], lon5, lon6 ); }
}





/*===========================================================================
  === Functions operating on the Ppath structure
  ===========================================================================*/

//! ppath_init_structure
/*!
   Initiates a Ppath structure to hold the given number of points.

   All fields releated with the ground, symmetry and tangent point are set
   to 0 or empty. The background field is set to background case 0. The
   constant field is set to -1. The refraction field is set to 0.

   The length of the l_step field is set to np-1.

   \param   ppath            Output: A Ppath structure.
   \param   atmosphere_dim   The atmospheric dimensionality.
   \param   np               Number of points of the path.

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_init_structure( 
              Ppath&      ppath,
        const Index&      atmosphere_dim,
        const Index&      np )
{
  assert( atmosphere_dim >= 1 );
  assert( atmosphere_dim <= 3 );

  ppath.dim        = atmosphere_dim;
  ppath.np         = np;
  ppath.refraction = 0;
  ppath.method     = "-";
  ppath.constant   = -1;   
  if( atmosphere_dim < 3 )
    {
      ppath.pos.resize( np, 2 );
      ppath.los.resize( np, 1 );
    }
  else
    {
      ppath.pos.resize( np, atmosphere_dim );
      ppath.los.resize( np, 2 );
      ppath.gp_lon.resize( np );
    }
  ppath.gp_p.resize( np );
  if( atmosphere_dim >= 2 )
    { ppath.gp_lat.resize( np ); }
  ppath.z.resize( np );
  if( np > 0 )
    { ppath.l_step.resize( np-1 ); }
  else
    { ppath.l_step.resize( 0 ); }
  ppath_set_background( ppath, 0 );
  ppath.tan_pos.resize(0);
  ppath.geom_tan_pos.resize(0);
}



//! ppath_set_background 
/*!
   Sets the background field of a Ppath structure.

   The different background cases have a number coding to simplify a possible
   change of the strings and checking of the what case that is valid.

   The case numbers are:                    <br>
      0. Not yet set.                       <br>
      1. Space.                             <br>
      2. The ground.                        <br>
      3. The surface of the cloud box.      <br>
      4. The interior of the cloud box.     

   \param   ppath            Output: A Ppath structure.
   \param   case_nr          Case number (see above)

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_set_background( 
              Ppath&      ppath,
        const Index&      case_nr )
{
  switch ( case_nr )
    {
    case 0:
      ppath.background = "";
      break;
    case 1:
      ppath.background = "space";
      break;
    case 2:
      ppath.background = "ground";
      break;
    case 3:
      ppath.background = "cloud box surface";
      break;
    case 4:
      ppath.background = "cloud box interior";
      break;
    default:
      ostringstream os;
      os << "Case number " << case_nr << " is not defined.";
      throw runtime_error(os.str());
    }
}



//! ppath_what_background
/*!
   Returns the case number for the radiative background.

   See further the function *ppath_set_background*.

   \return                   The case number.
   \param   ppath            A Ppath structure.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Index ppath_what_background( const Ppath&   ppath )
{
  if( ppath.background == "" )
    { return 0; }
  else if( ppath.background == "space" )
    { return 1; }
  else if( ppath.background == "ground" )
    { return 2; }
  else if( ppath.background == "cloud box surface" )
    { return 3; }
  else if( ppath.background == "cloud box interior" )
    { return 4; }
  else
    {
      ostringstream os;
      os << "The string " << ppath.background 
         << " is not a valid background case.";
      throw runtime_error(os.str());
    }
}



//! ppath_copy
/*!
   Copy the content in ppath2 to ppath1.

   The ppath1 structure must be allocated before calling the function. The
   structure can be allocated to hold more points than found in ppath2.
   The data of ppath2 is placed in the first positions of ppath1.

   \param   ppath1    Output: PPath structure.
   \param   ppath2    The Ppath structure to be copied.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_copy(
           Ppath&      ppath1,
     const Ppath&      ppath2 )
{
  assert( ppath1.np >= ppath2.np ); 

  // The field np shall not be copied !!!

  ppath1.dim        = ppath2.dim;
  ppath1.refraction = ppath2.refraction;
  ppath1.method     = ppath2.method;
  ppath1.constant   = ppath2.constant;
  ppath1.background = ppath2.background;

  for( Index i=0; i<ppath2.np; i++ )
    {
      ppath1.pos(i,0)      = ppath2.pos(i,0);
      ppath1.pos(i,1)      = ppath2.pos(i,1);
      ppath1.los(i,0)      = ppath2.los(i,0);
      ppath1.z[i]          = ppath2.z[i];
      gridpos_copy( ppath1.gp_p[i], ppath2.gp_p[i] );
      
      if( ppath1.dim >= 2 )
        { gridpos_copy( ppath1.gp_lat[i], ppath2.gp_lat[i] ); }
      
      if( ppath1.dim == 3 )
        {
          ppath1.pos(i,2)        = ppath2.pos(i,2);
          ppath1.los(i,1)        = ppath2.los(i,1);
          gridpos_copy( ppath1.gp_lon[i], ppath2.gp_lon[i] ); 
        }
      
      if( i > 0 )
        { ppath1.l_step[i-1] = ppath2.l_step[i-1]; }
     }

  if( ppath2.tan_pos.nelem() )
    {
      ppath1.tan_pos.resize( ppath2.tan_pos.nelem() );
      ppath1.tan_pos = ppath2.tan_pos; 
    }
  if( ppath2.geom_tan_pos.nelem() )
    {
      ppath1.geom_tan_pos.resize( ppath2.geom_tan_pos.nelem() );
      ppath1.geom_tan_pos = ppath2.geom_tan_pos; 
    }
}



//! ppath_append
/*!
   Combines two Ppath structures   

   The function appends a Ppath structure to another structure. 
 
   All the data of ppath1 is kept.

   The first point in ppath2 is assumed to be the same as the last in ppath1.
   Only data in ppath from the fields pos, los, z, l_step, gp_XXX and 
   background are considered.

   \param   ppath1    Output: Ppath structure to be expanded.
   \param   ppath2    The Ppath structure to include in ppath.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_append(
           Ppath&   ppath1,
     const Ppath&   ppath2 )
{
  const Index n1 = ppath1.np;
  const Index n2 = ppath2.np;

  Ppath   ppath;
  ppath_init_structure( ppath, ppath1.dim, n1 );
  ppath_copy( ppath, ppath1 );

  ppath_init_structure( ppath1, ppath1.dim, n1 + n2 - 1 );
  ppath_copy( ppath1, ppath );

  // Append data from ppath2
  Index i1;
  for( Index i=1; i<n2; i++ )
    {
      i1 = n1 + i - 1;

      ppath1.pos(i1,0)      = ppath2.pos(i,0);
      ppath1.pos(i1,1)      = ppath2.pos(i,1);
      ppath1.los(i1,0)      = ppath2.los(i,0);
      ppath1.z[i1]          = ppath2.z[i];
      gridpos_copy( ppath1.gp_p[i1], ppath2.gp_p[i] );

      if( ppath1.dim >= 2 )
        { gridpos_copy( ppath1.gp_lat[i1], ppath2.gp_lat[i] ); }
      
      if( ppath1.dim == 3 )
        {
          ppath1.pos(i1,2)        = ppath2.pos(i,2);
          ppath1.los(i1,1)        = ppath2.los(i,1);
          gridpos_copy( ppath1.gp_lon[i1], ppath2.gp_lon[i] ); 
        }
      
      ppath1.l_step[i1-1] = ppath2.l_step[i-1];
    }

  if( ppath_what_background( ppath2 ) )
    { ppath1.background = ppath2.background; }
}



//! ppath_fill_1d
/*!
   Fills a 1D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step and gp_p.
   The pressure grid positions (gp_p) are filtered through gridpos_check_fd.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath      Output: Ppath structure.
   \param   r          Vector with radius for the path points.
   \param   lat        Vector with latitude for the path points.
   \param   za         Vector with zenith angle for the path points.
   \param   lstep      Length along the path between the points.
   \param   r_geoid    Geoid radii.
   \param   z_field    Geometrical altitudes.
   \param   ip         Pressure grid range.

   \author Patrick Eriksson
   \date   2002-07-18
*/
void ppath_fill_1d(
           Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   za,
     ConstVectorView   lstep,
     const double&     r_geoid,
     ConstVectorView   z_field,
     const Index&      ip )
{
  // Help variables that are common for all points.
  const double   r1 = r_geoid + z_field[ip];
  const double   dr = z_field[ip+1] - z_field[ip];

  for( Index i=0; i<r.nelem(); i++ )
    {
      ppath.pos(i,0) = r[i];
      ppath.pos(i,1) = lat[i];
      ppath.los(i,0) = za[i];
      
      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( r[i] - r1 ) / dr;
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      ppath.z[i] = r[i] - r_geoid;

      if( i > 0 )
        { ppath.l_step[i-1] = lstep[i-1]; }
    }
}



//! ppath_fill_2d
/*!
   Fills a 2D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step, gp_p and gp_lat.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath      Output: Ppath structure.
   \param   r          Vector with radius for the path points.
   \param   lat        Vector with latitude for the path points.
   \param   za         Vector with zenith angle for the path points.
   \param   lstep      Length along the path between the points.
   \param   r_geoid    Geoid radii.
   \param   z_field    Geometrical altitudes
   \param   lat_grid   Latitude grid.
   \param   ip         Pressure grid range.
   \param   ilat       Latitude grid range.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_fill_2d(
           Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   za,
     const double&     lstep,
     ConstVectorView   r_geoid,
     ConstMatrixView   z_field,
     ConstVectorView   lat_grid,
     const Index&      ip,
     const Index&      ilat )
{
  // Help variables that are common for all points.
  const double   dlat  = lat_grid[ilat+1] - lat_grid[ilat];
  const double   r1low = r_geoid[ilat] + z_field(ip,ilat);
  const double   drlow = r_geoid[ilat+1] + z_field(ip,ilat+1) -r1low;
  const double   r1upp = r_geoid[ilat] + z_field(ip+1,ilat);
  const double   drupp = r_geoid[ilat+1] + z_field(ip+1,ilat+1) - r1upp;
  const double   z1low =  z_field(ip,ilat);
  const double   dzlow =  z_field(ip,ilat+1) -z1low;
  const double   z1upp =  z_field(ip+1,ilat);
  const double   dzupp =  z_field(ip+1,ilat+1) - z1upp;

  for( Index i=0; i<r.nelem(); i++ )
    {
      ppath.pos(i,0) = r[i];
      ppath.pos(i,1) = lat[i];
      ppath.los(i,0) = za[i];
      
      // Weigt in the latitude direction
      double w = ( lat[i] - lat_grid[ilat] ) / dlat;

      // Radius of lower and upper face at present latitude
      const double rlow = r1low + w * drlow;
      const double rupp = r1upp + w * drupp;

      // Geometrical altitude of lower and upper face at present latitude
      const double zlow = z1low + w * dzlow;
      const double zupp = z1upp + w * dzupp;

      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( r[i] - rlow ) / ( rupp - rlow );
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      ppath.z[i] = zlow + ppath.gp_p[i].fd[0] * ( zupp -zlow );

      ppath.gp_lat[i].idx   = ilat;
      ppath.gp_lat[i].fd[0] = ( lat[i] - lat_grid[ilat] ) / dlat;
      ppath.gp_lat[i].fd[1] = 1 - ppath.gp_lat[i].fd[0];
      gridpos_check_fd( ppath.gp_lat[i] );

      if( i > 0 )
        { ppath.l_step[i-1] = lstep; }
    }
}



//! ppath_fill_3d
/*!
   Fills a 3D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step, gp_p, gp_lat and gp_lon.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath      Output: Ppath structure.
   \param   r          Vector with radius for the path points.
   \param   lat        Vector with latitude for the path points.
   \param   lon        Vector with longitude for the path points.
   \param   za         Vector with zenith angle for the path points.
   \param   aa         Vector with azimuth angle for the path points.
   \param   lstep      Length along the path between the points.
   \param   r_geoid    Geoid radii.
   \param   z_field    Geometrical altitudes
   \param   lat_grid   Latitude grid.
   \param   lon_grid   Longitude grid.
   \param   ip         Pressure grid range.
   \param   ilat       Latitude grid range.
   \param   ilon       Longitude grid range.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_fill_3d(
           Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   lon,
     ConstVectorView   za,
     ConstVectorView   aa,
     const double&     lstep,
     ConstMatrixView   r_geoid,
     ConstTensor3View  z_field,
     ConstVectorView   lat_grid,
     ConstVectorView   lon_grid,
     const Index&      ip,
     const Index&      ilat,
     const Index&      ilon )
{
  // Help variables that are common for all points.
  const double   lat1  = lat_grid[ilat];
  const double   lat3  = lat_grid[ilat+1];
  const double   lon5  = lon_grid[ilon];
  const double   lon6  = lon_grid[ilon+1];
  const double   r15a  = r_geoid(ilat,ilon) + z_field(ip,ilat,ilon);
  const double   r35a  = r_geoid(ilat+1,ilon) + z_field(ip,ilat+1,ilon); 
  const double   r36a  = r_geoid(ilat+1,ilon+1) + z_field(ip,ilat+1,ilon+1); 
  const double   r16a  = r_geoid(ilat,ilon+1) + z_field(ip,ilat,ilon+1);
  const double   r15b  = r_geoid(ilat,ilon) + z_field(ip+1,ilat,ilon);
  const double   r35b  = r_geoid(ilat+1,ilon) + z_field(ip+1,ilat+1,ilon); 
  const double   r36b  = r_geoid(ilat+1,ilon+1) + z_field(ip+1,ilat+1,ilon+1);
  const double   r16b  = r_geoid(ilat,ilon+1) + z_field(ip+1,ilat,ilon+1);
  const double   rg15  = r_geoid(ilat,ilon);
  const double   rg35  = r_geoid(ilat+1,ilon); 
  const double   rg36  = r_geoid(ilat+1,ilon+1);
  const double   rg16  = r_geoid(ilat,ilon+1);     
  const double   dlat  = lat3 - lat1;
  const double   dlon  = lon6 - lon5;

  for( Index i=0; i<r.nelem(); i++ )
    {
      // Radius of pressure surfaces at present lat and lon
      const double   rlow = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                      r15a, r35a, r36a, r16a, lat[i], lon[i] );
      const double   rupp = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                      r15b, r35b, r36b, r16b, lat[i], lon[i] );

      // Position and LOS. 
      ppath.pos(i,0) = r[i];
      ppath.pos(i,1) = lat[i];
      ppath.pos(i,2) = lon[i];
      ppath.los(i,0) = za[i];
      ppath.los(i,1) = aa[i];
      
      // Pressure grid index
      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( ppath.pos(i,0) - rlow ) / ( rupp - rlow );
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      // Geometrical altitude
      const double   rgeoid = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                      rg15, rg35, rg36, rg16, lat[i], lon[i] );
      const double   zlow = rlow - rgeoid;
      const double   zupp = rupp - rgeoid;
      //
      ppath.z[i] = zlow + ppath.gp_p[i].fd[0] * ( zupp -zlow );

      // Latitude grid index
      ppath.gp_lat[i].idx   = ilat;
      ppath.gp_lat[i].fd[0] = ( lat[i] - lat1 ) / dlat;
      ppath.gp_lat[i].fd[1] = 1 - ppath.gp_lat[i].fd[0];
      gridpos_check_fd( ppath.gp_lat[i] );

      // Longitude grid index
      //
      // The longitude  is undefined at the poles. The grid index is set to
      // the start point.
      //
      if( abs( lat[i] ) < 90 )
        {
          ppath.gp_lon[i].idx   = ilon;
          ppath.gp_lon[i].fd[0] = ( lon[i] - lon5 ) / dlon;
          ppath.gp_lon[i].fd[1] = 1 - ppath.gp_lon[i].fd[0];
          gridpos_check_fd( ppath.gp_lon[i] );
        }
      else
        {
          ppath.gp_lon[i].idx   = 0;
          ppath.gp_lon[i].fd[0] = 0;
          ppath.gp_lon[i].fd[1] = 1;
        }

      if( i > 0 )
        { ppath.l_step[i-1] = lstep; }
    }
}





/*===========================================================================
  === Functions related to propagation paths with refraction
  ===========================================================================*/

//! refraction_ppc
/*! 
   Calculates the propagation path constant for cases with refraction.

   Both positive and negative zenith angles are handled.

   \return               Path constant.
   \param   r            Radius.
   \param   za           LOS Zenith angle.
   \param   refr_index   Refractive index.

   \author Patrick Eriksson
   \date   2002-05-17
*/
double refraction_ppc( 
        const double&   r, 
        const double&   za, 
        const double&   refr_index )
{
  assert( r > 0 );
  assert( abs(za) <= 180 );

  return r * refr_index * sin( DEG2RAD * abs(za) );
}





/*===========================================================================
  === Help functions for the *ppath_step* functions found below
  === These functions are mainly pieces of code that are common for at least
  === two functions (or two places in some function) and for this reason 
  === the headers are not complete. 
  ===========================================================================*/


//! ppath_start_1d
/*! 
   Internal help function for 1D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-13
*/
void ppath_start_1d(
              double&     r_start,
              double&     lat_start,
              double&     za_start,
              Index&      ip,
        const Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        const double&     r_geoid,
        const double&     z_ground )
{
  // Number of points in the incoming ppath
  const Index   imax = ppath.np - 1;

  // Number of pressure levels
  const Index   npl = p_grid.nelem();

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.pos(imax,0);
  lat_start = ppath.pos(imax,1);
  za_start  = ppath.los(imax,0);

  // Asserts
  assert( npl >= 2 );
  assert( is_decreasing( p_grid ) );
  assert( is_size( z_field, npl ) );
  assert( is_increasing( z_field ) );
  assert( r_geoid > 0 );
  //
  assert( ppath.dim == 1 );
  assert( ppath.np >= 1 );
  assert( ppath.gp_p[imax].idx >= 0 );
  assert( ppath.gp_p[imax].idx <= ( npl - 2 ) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] <= 1 );
  //
  assert( r_start >= r_geoid + z_ground - 1e-6 );
  assert( za_start >= 0  &&  za_start <= 180 );


  // Determine index of the pressure surface being the lower limit for the
  // grid range of interest.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], za_start<=90 );

  out3 << "  pressure grid range  : " << ip << "\n";
}



//! ppath_end_1d
/*! 
   Internal help function for 1D path calculations.

   The function performs the end part of the calculations, that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-27
*/
void ppath_end_1d(
              Ppath&      ppath,
        ConstVectorView   r_v,
        ConstVectorView   lat_v,
        ConstVectorView   za_v,
        const double&     lstep,
        ConstVectorView   z_field,
        const double&     r_geoid,
        const Index&      ip,
        const Index&      endface,
        const Index&      tanpoint,
        const String&     method,
        const Index&      refraction,
        const double&     ppc )
{
  out3 << "  end face number code : " << endface << "\n";

  // Number of path points
  const Index   np = r_v.nelem();

  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure(  ppath, 1, np );
  //
  ppath.method     = method;
  ppath.refraction = refraction;
  ppath.constant   = ppc;
  //
  ppath_fill_1d( ppath, r_v, lat_v, za_v, Vector(np-1,lstep), r_geoid, 
                                                                 z_field, ip );


  // Different options depending on position of end point of step:

  //--- End point is the ground
  if( endface == 7 )
    { ppath_set_background( ppath, 2 ); }

  //--- End point is a tangent point
  else if( tanpoint )
    {
      ppath.tan_pos.resize(2);
      ppath.tan_pos[0] = r_v[np-1];
      ppath.tan_pos[1] = lat_v[np-1];
    }

  //--- End point is on top of a pressure surface
  else
    {
      gridpos_force_end_fd( ppath.gp_p[np-1] );
    }
}



//! ppath_start_2d
/*! 
   Internal help function for 2D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-18
*/
void ppath_start_2d(
              double&     r_start,
              double&     lat_start,
              double&     za_start,
              Index&      ip,
              Index&      ilat,
              double&     lat1,
              double&     lat3,
              double&     r1a,
              double&     r3a,
              double&     r3b,
              double&     r1b,
              double&     rground1,
              double&     rground3,
        const Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground )
{
  // Number of points in the incoming ppath
  const Index imax = ppath.np - 1;

  // Number of pressure levels and latitudes
  const Index npl = p_grid.nelem();
  const Index nlat = lat_grid.nelem();

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.pos(imax,0);
  lat_start = ppath.pos(imax,1);
  za_start  = ppath.los(imax,0);

  // First asserts (more below)
  assert( npl >= 2 );
  assert( is_decreasing( p_grid ) );
  assert( nlat >= 2 );
  assert( is_increasing( lat_grid ) );
  assert( is_size( z_field, npl, nlat ) );
  assert( is_size( r_geoid, nlat ) );
  assert( is_size( z_ground, nlat ) );
  //
  assert( ppath.dim == 2 );
  assert( ppath.np >= 1 );
  assert( ppath.gp_p[imax].idx >= 0 );
  assert( ppath.gp_p[imax].idx <= ( npl - 2 ) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] <= 1 );
  assert( ppath.gp_lat[imax].idx >= 0 );
  assert( ppath.gp_lat[imax].idx <= ( nlat - 2 ) );
  assert( ppath.gp_lat[imax].fd[0] >= 0 );
  assert( ppath.gp_lat[imax].fd[0] <= 1 );
  //
  assert( za_start >= -180  );
  assert( za_start <= 180 );
  //
  // more asserts below ...


  // Determine interesting latitude grid range and latitude end points of 
  // the range.
  //
  ilat = gridpos2gridrange( ppath.gp_lat[imax], za_start >= 0 );
  //
  lat1 = lat_grid[ilat];
  lat3 = lat_grid[ilat+1];

  // Determine interesting pressure grid range. Do this first assuming that
  // the pressure surfaces are not tilted (that is, abs(za_start<=90) always
  // mean upward observation). 
  // Set radius for the corners of the grid cell and the radial slope of 
  // pressure surface limits of the grid cell to match the found ip.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], abs(za_start) <= 90);
  //
  r1a = r_geoid[ilat] + z_field(ip,ilat);        // lower-left
  r3a = r_geoid[ilat+1] + z_field(ip,ilat+1);    // lower-right
  r3b = r_geoid[ilat+1] + z_field(ip+1,ilat+1);  // upper-right
  r1b = r_geoid[ilat] + z_field(ip+1,ilat);      // upper-left
  
  // Slopes of pressure surfaces
  double   c2 = psurface_slope_2d( lat1, lat3, r1a, r3a );
  double   c4 = psurface_slope_2d( lat1, lat3, r1b, r3b );

  // Check if the LOS zenith angle happen to be between 90 and the zenith angle
  // of the pressure surface (that is, 90 + tilt of pressure surface), and in
  // that case if ip must be changed. This check is only needed when the
  // start point is on a pressure surface.
  //
  if( is_gridpos_at_index_i( ppath.gp_p[imax], ip )  )
    {
      double tilt = psurface_angletilt( r_start, c2 );

      if( is_los_downwards( za_start, tilt ) )
        {
          ip--;
          r1b = r1a;   r3b = r3a;   c4 = c2;
          r1a = r_geoid[ilat] + z_field(ip,ilat);
          r3a = r_geoid[ilat+1] + z_field(ip,ilat+1);
          c2 = psurface_slope_2d( lat1, lat3, r1a, r3a );
        }
    }
  else if( is_gridpos_at_index_i( ppath.gp_p[imax], ip+1 )  )
    {
      double tilt = psurface_angletilt( r_start, c4 );

      if( !is_los_downwards( za_start, tilt ) )
        {
          ip++;
          r1a = r1b;   r3a = r3b;   c2 = c4;
          r3b = r_geoid[ilat+1] + z_field(ip+1,ilat+1);
          r1b = r_geoid[ilat] + z_field(ip+1,ilat);    
          c4 = psurface_slope_2d( lat1, lat3, r1b, r3b );
        }
    }

  out3 << "  pressure grid range  : " << ip << "\n";
  out3 << "  latitude grid range  : " << ilat << "\n";

  // Ground radius at latitude end points
  rground1 = r_geoid[ilat] + z_ground[ilat];
  rground3 = r_geoid[ilat+1] + z_ground[ilat+1];
}



//! ppath_end_2d
/*! 
   Internal help function for 2D path calculations.

   The function performs the end part of the calculations, that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-29
*/
void ppath_end_2d(
              Ppath&      ppath,
        ConstVectorView   r_v,
        ConstVectorView   lat_v,
        ConstVectorView   za_v,
        const double&     lstep,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   r_geoid,
        const Index&      ip,
        const Index&      ilat,
        const Index&      endface,
        const Index&      tanpoint,
        const String&     method,
        const Index&      refraction,
        const double&     ppc )
{
  // Number of path points
  const Index   np = r_v.nelem();

  out3 << "  end face number code : " << endface << "\n";

  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure(  ppath, 2, np );
  //
  ppath.method     = method;
  ppath.refraction = refraction;
  ppath.constant   = ppc;
  //
  ppath_fill_2d( ppath, r_v, lat_v, za_v, lstep, r_geoid, z_field, lat_grid, 
                                                                    ip, ilat );

  // Do end-face specific tasks
  if( endface == 7 )
    { ppath_set_background( ppath, 2 ); }
  else if( tanpoint )
    {
      ppath.tan_pos.resize(2);
      ppath.tan_pos[0] = r_v[np-1];
      ppath.tan_pos[1] = lat_v[np-1];
    }

  // Set fractional distance for end point
  //
  if( endface == 1  ||  endface == 3 )
    { gridpos_force_end_fd( ppath.gp_lat[np-1] ); }
  else if( endface == 2  ||  endface == 4 )
    { gridpos_force_end_fd( ppath.gp_p[np-1] ); }
}



//! ppath_start_3d
/*! 
   Internal help function for 3D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_start_3d(
              double&     r_start,
              double&     lat_start,
              double&     lon_start,
              double&     za_start,
              double&     aa_start,
              Index&      ip,
              Index&      ilat,
              Index&      ilon,
              double&     lat1,
              double&     lat3,
              double&     lon5,
              double&     lon6,
              double&     r15a,
              double&     r35a,
              double&     r36a,
              double&     r16a,
              double&     r15b,
              double&     r35b,
              double&     r36b,
              double&     r16b,
              double&     rground15,
              double&     rground35,
              double&     rground36,
              double&     rground16,
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstTensor3View  z_field,
        ConstMatrixView   r_geoid,
        ConstMatrixView   z_ground )
{
  // Number of points in the incoming ppath
  const Index imax = ppath.np - 1;

  // Number of pressure levels and latitudes
  const Index npl = p_grid.nelem();
  const Index nlat = lat_grid.nelem();
  const Index nlon = lon_grid.nelem();

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.pos(imax,0);
  lat_start = ppath.pos(imax,1);
  lon_start = ppath.pos(imax,2);
  za_start  = ppath.los(imax,0);
  aa_start  = ppath.los(imax,1);

  // First asserts (more below)
  assert( npl >= 2 );
  assert( is_decreasing( p_grid ) );
  assert( nlat >= 2 );
  assert( nlon >= 2 );
  assert( is_increasing( lat_grid ) );
  assert( is_increasing( lon_grid ) );
  assert( is_size( z_field, npl, nlat, nlon ) );
  assert( is_size( r_geoid, nlat, nlon ) );
  assert( is_size( z_ground, nlat, nlon ) );
  //
  assert( ppath.dim == 3 );
  assert( ppath.np >= 1 );
  assert( ppath.gp_p[imax].idx >= 0 );
  assert( ppath.gp_p[imax].idx <= ( npl - 2 ) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] <= 1 );
  assert( ppath.gp_lat[imax].idx >= 0 );
  assert( ppath.gp_lat[imax].idx <= ( nlat - 2 ) );
  assert( ppath.gp_lat[imax].fd[0] >= 0 );
  assert( ppath.gp_lat[imax].fd[0] <= 1 );
  assert( ppath.gp_lon[imax].idx >= 0 );
  assert( ppath.gp_lon[imax].idx <= ( nlon - 2 ) );
  assert( ppath.gp_lon[imax].fd[0] >= 0 );
  assert( ppath.gp_lon[imax].fd[0] <= 1 );
  //
  assert( za_start >= 0 );
  assert( za_start <= 180 );
  assert( aa_start >= -180 );
  assert( aa_start <= 180 );


  // Lower index of lat and lon ranges of interest
  //
  // The longitude is undefined at the poles and as the azimuth angle
  // is defined in other way at the poles.
  //
  if( lat_start == 90 )
    { 
      ilat = nlat - 2;
      GridPos   gp_tmp;
      gridpos( gp_tmp, lon_grid, aa_start );
      if( aa_start < 180 )
        { ilon = gridpos2gridrange( gp_tmp, 1 ); }
      else
        { ilon = gridpos2gridrange( gp_tmp, 0 ); }
    }
  else if( lat_start == -90 )
    { 
      ilat = 0; 
      GridPos   gp_tmp;
      gridpos( gp_tmp, lon_grid, aa_start );
      if( aa_start < 180 )
        { ilon = gridpos2gridrange( gp_tmp, 1 ); }
      else
        { ilon = gridpos2gridrange( gp_tmp, 0 ); }
    }
  else
    { 
      if( lat_start > 0 )
        ilat = gridpos2gridrange( ppath.gp_lat[imax], abs( aa_start ) < 90 );
      else
        ilat = gridpos2gridrange( ppath.gp_lat[imax], abs( aa_start ) <= 90 );
      if( lon_start < lon_grid[nlon-1] )
        { ilon = gridpos2gridrange( ppath.gp_lon[imax], aa_start >= 0 ); }
      else
        { ilon = nlon - 2; }
    }
  //
  lat1 = lat_grid[ilat];
  lat3 = lat_grid[ilat+1];
  lon5 = lon_grid[ilon];
  lon6 = lon_grid[ilon+1];

  // Determine interesting pressure grid range. Do this first assuming that
  // the pressure surfaces are not tilted (that is, abs(za_start<=90) always
  // mean upward observation). 
  // Set radius for the corners of the grid cell and the radial slope of 
  // pressure surface limits of the grid cell to match the found ip.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], za_start <= 90 );
  //
  r15a = r_geoid(ilat,ilon) + z_field(ip,ilat,ilon);
  r35a = r_geoid(ilat+1,ilon) + z_field(ip,ilat+1,ilon); 
  r36a = r_geoid(ilat+1,ilon+1) + z_field(ip,ilat+1,ilon+1); 
  r16a = r_geoid(ilat,ilon+1) + z_field(ip,ilat,ilon+1);
  r15b = r_geoid(ilat,ilon) + z_field(ip+1,ilat,ilon);
  r35b = r_geoid(ilat+1,ilon) + z_field(ip+1,ilat+1,ilon); 
  r36b = r_geoid(ilat+1,ilon+1) + z_field(ip+1,ilat+1,ilon+1); 
  r16b = r_geoid(ilat,ilon+1) + z_field(ip+1,ilat,ilon+1);

  // This part is a fix to catch start postions on top of a pressure surface
  // that does not have an end fractional distance for the first step.
  {
    // Radius of lower and upper pressure surface at the start position
    const double   rlow = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
				     r15a, r35a, r36a, r16a, lat_start, lon_start );
    const double   rupp = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
				     r15b, r35b, r36b, r16b, lat_start, lon_start );
    if (abs(r_start-rlow) < RTOL || abs(r_start-rupp) < RTOL)
      { gridpos_force_end_fd( ppath.gp_p[imax] ); }
  }

  // Check if the LOS zenith angle happen to be between 90 and the zenith angle
  // of the pressure surface (that is, 90 + tilt of pressure surface), and in
  // that case if ip must be changed. This check is only needed when the
  // start point is on a pressure surface.
  //
  if( is_gridpos_at_index_i( ppath.gp_p[imax], ip )  )
    {
      // Slope and angular tilt of lower pressure surface
      double c2 = psurface_slope_3d( lat1, lat3, lon5, lon6, 
                      r15a, r35a, r36a, r16a, lat_start, lon_start, aa_start );
      double tilt = psurface_angletilt( r_start, c2 );

      if( is_los_downwards( za_start, tilt ) )
        {
          ip--;
          r15b = r15a;   r35b = r35a;   r36b = r36a;   r16b = r16a;
          r15a = r_geoid(ilat,ilon) + z_field(ip,ilat,ilon);
          r35a = r_geoid(ilat+1,ilon) + z_field(ip,ilat+1,ilon); 
          r36a = r_geoid(ilat+1,ilon+1) + z_field(ip,ilat+1,ilon+1); 
          r16a = r_geoid(ilat,ilon+1) + z_field(ip,ilat,ilon+1);
        }
    }
  else if( is_gridpos_at_index_i( ppath.gp_p[imax], ip+1 )  )
    {
      // Slope and angular tilt of lower pressure surface
      double c4 = psurface_slope_3d( lat1, lat3 ,lon5, lon6, 
                      r15b, r35b, r36b, r16b, lat_start, lon_start, aa_start );
      double tilt = psurface_angletilt( r_start, c4 );

      if( !is_los_downwards( za_start, tilt ) )
        {
          ip++;
          r15a = r15b;   r35a = r35b;   r36a = r36b;   r16a = r16b;
          r15b = r_geoid(ilat,ilon) + z_field(ip+1,ilat,ilon);
          r35b = r_geoid(ilat+1,ilon) + z_field(ip+1,ilat+1,ilon); 
          r36b = r_geoid(ilat+1,ilon+1) + z_field(ip+1,ilat+1,ilon+1); 
          r16b = r_geoid(ilat,ilon+1) + z_field(ip+1,ilat,ilon+1);
        }
    }

  out3 << "  pressure grid range  : " << ip << "\n";
  out3 << "  latitude grid range  : " << ilat << "\n";
  out3 << "  longitude grid range : " << ilon << "\n";

  // Ground radius at latitude/longitude corner points
  rground15 = r_geoid(ilat,ilon) + z_ground(ilat,ilon);
  rground35 = r_geoid(ilat+1,ilon) + z_ground(ilat+1,ilon);
  rground36 = r_geoid(ilat+1,ilon+1) + z_ground(ilat+1,ilon+1);
  rground16 = r_geoid(ilat,ilon+1) + z_ground(ilat,ilon+1);
}



//! ppath_end_3d
/*! 
   Internal help function for 3D path calculations.

   The function performs the end part of the calculations, that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_end_3d(
              Ppath&      ppath,
        ConstVectorView   r_v,
        ConstVectorView   lat_v,
        ConstVectorView   lon_v,
        ConstVectorView   za_v,
        ConstVectorView   aa_v,
        const double&     lstep,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstTensor3View  z_field,
        ConstMatrixView   r_geoid,
        const Index&      ip,
        const Index&      ilat,
        const Index&      ilon,
        const Index&      endface,
        const Index&      tanpoint,
        const String&     method,
        const Index&      refraction,
        const double&     ppc )
{
  // Number of path points
  const Index   np = r_v.nelem();

  out3 << "  end face number code : " << endface << "\n";

  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure(  ppath, 3, np );
  //
  ppath.method     = method;
  ppath.refraction = refraction;
  ppath.constant   = ppc;
  //
  ppath_fill_3d( ppath, r_v, lat_v, lon_v, za_v, aa_v, lstep, 
                        r_geoid, z_field, lat_grid, lon_grid, ip, ilat, ilon );

  // Do end-face specific tasks
  if( endface == 7 )
    { ppath_set_background( ppath, 2 ); }
  if( tanpoint )
    {
      ppath.tan_pos.resize(3);
      ppath.tan_pos[0] = r_v[np-1];
      ppath.tan_pos[1] = lat_v[np-1];
      ppath.tan_pos[2] = lon_v[np-1];
    }

  // Set fractional distance for end point
  //
  if( endface == 1  ||  endface == 3 )
    { gridpos_force_end_fd( ppath.gp_lat[np-1] ); }
  else if( endface == 2  ||  endface == 4 )
    { gridpos_force_end_fd( ppath.gp_p[np-1] ); }
  else if( endface == 5  ||  endface == 6 )
    { gridpos_force_end_fd( ppath.gp_lon[np-1] ); }
}



//! interpolate_raytracing_points
/*! 
   Interpolates a set of ray tracing points to a set of points linearly
   spaced along the path.

   All quantities are interpolated linearly.

   Empty vectors can be sent as input for *lon_rt* and *aa_rt* for 1D and 2D.
   The output vectors *lon_v* and *aa_v* are then not filled.

   \author Patrick Eriksson
   \date   2002-11-27
*/
void interpolate_raytracing_points(
             Vector&     r_v,
             Vector&     lat_v,
             Vector&     lon_v,
             Vector&     za_v,
             Vector&     aa_v,
             double&     lstep,
       ConstVectorView   r_rt,
       ConstVectorView   lat_rt,
       ConstVectorView   lon_rt,
       ConstVectorView   za_rt,
       ConstVectorView   aa_rt,
       ConstVectorView   l_rt,
       const double&     lmax )
{
  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. If *lmax* <= 0, then only the end points
  // shall be kept.
  //
  const Index     nrt = r_rt.nelem();
        Index     n = 2;
        double   ltotsum = l_rt.sum();
  //
  if( lmax > 0 )
    { n = Index( ceil( ltotsum / lmax ) ) + 1; }
  //
  r_v.resize(n);
  lat_v.resize(n);
  za_v.resize(n);
  //
  if( n == 2 )
    {
      r_v[0] = r_rt[0];     r_v[1] = r_rt[nrt-1];
      za_v[0] = za_rt[0];   za_v[1] = za_rt[nrt-1];
      lat_v[0] = lat_rt[0]; lat_v[1] = lat_rt[nrt-1];
      lstep = ltotsum;
      if( lon_rt.nelem() > 0 )
        {
            lon_v.resize(n); lon_v[0] = lon_rt[0]; lon_v[1] = lon_rt[nrt-1];
            aa_v.resize(n);  aa_v[0] = aa_rt[0];   aa_v[1] = aa_rt[nrt-1];
        }
    }
  else
    {
      Vector           lsum(nrt), llin(n);
      ArrayOfGridPos   gp(n);
      Matrix           itw(n,2);
      
      lsum[0] = 0;
      for( Index i=1; i<nrt; i++ )
        { lsum[i] = lsum[i-1] + l_rt[i-1]; }

      nlinspace( llin, 0, ltotsum, n );

      gridpos( gp, lsum, llin );
      gridpos_force_end_fd( gp[0] );
      gridpos_force_end_fd( gp[n-1] );

      interpweights( itw, gp );

      interp( r_v, itw, r_rt, gp );
      interp( za_v, itw, za_rt, gp );
      interp( lat_v, itw, lat_rt, gp );
      lstep = llin[1] - llin[0];

      if( lon_rt.nelem() > 0 )
        {
            lon_v.resize(n); interp( lon_v, itw, lon_rt, gp );
            aa_v.resize(n);  interp( aa_v, itw, aa_rt, gp );
        }
    }
}



//! from_raytracingarrays_to_ppath_vectors_1d_and_2d
/*! 
   A small help function to convert arrays with ray tracing points t
   interpolated values along the path.

   This function is common for 1D and 2D.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void from_raytracingarrays_to_ppath_vectors_1d_and_2d(
             Vector&          r_v,
             Vector&          lat_v,
             Vector&          za_v,
             double&          lstep,
       const Array<double>&   r_array,
       const Array<double>&   lat_array,
       const Array<double>&   za_array,
       const Array<double>&   l_array,
       const Index&           reversed,
       const double&          lmax )
{
  // Copy arrays to vectors for later interpolation
  //
  // Number of ray tracing points
  const Index   nrt=r_array.nelem();
  //
  Vector    r_rt(nrt), lat_rt(nrt), za_rt(nrt), l_rt(nrt-1);    
  //
  if( !reversed )
    {
      for( Index i=0; i<nrt; i++ )
        {
          r_rt[i]   = r_array[i];
          lat_rt[i] = lat_array[i]; 
          za_rt[i]  = za_array[i];
          if( i < (nrt-1) )
            { l_rt[i]   = l_array[i]; }
        }
    }
  else
    {
      for( Index i=0; i<nrt; i++ )
        {
          const Index i1 = nrt - 1 - i;
          r_rt[i]   = r_array[i1];
          lat_rt[i] = lat_array[0]+ lat_array[nrt-1] - lat_array[i1]; 
          za_rt[i]  = 180 - za_array[i1];
          if( i < (nrt-1) )
            { l_rt[i]   = l_array[i1-1]; }
        }
    }

  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. If *lmax* <= 0, then only the end points
  //
  Vector    dummy;   // Dummy vector for output variables not used
  //
  interpolate_raytracing_points( r_v, lat_v, dummy, za_v, dummy, lstep,
      r_rt, lat_rt, Vector(0), za_rt, Vector(0), l_rt, lmax );
}



//! from_raytracingarrays_to_ppath_vectors_3d
/*! 
   A small help function to convert arrays with ray tracing points to
   interpolated values along the path.

   \author Patrick Eriksson
   \date   2003-01-18
*/
void from_raytracingarrays_to_ppath_vectors_3d(
             Vector&          r_v,
             Vector&          lat_v,
             Vector&          lon_v,
             Vector&          za_v,
             Vector&          aa_v,
             double&          lstep,
       const Array<double>&   r_array,
       const Array<double>&   lat_array,
       const Array<double>&   lon_array,
       const Array<double>&   za_array,
       const Array<double>&   aa_array,
       const Array<double>&   l_array,
       const double&          lmax )
{
  // Copy arrays to vectors for later interpolation
  //
  // Number of ray tracing points
  const Index   nrt=r_array.nelem();
  //
  Vector    r_rt(nrt), lat_rt(nrt), lon_rt(nrt);
  Vector    za_rt(nrt), aa_rt(nrt),l_rt(nrt-1);    
  //
  for( Index i=0; i<nrt; i++ )
    {
      r_rt[i]   = r_array[i];
      lat_rt[i] = lat_array[i]; 
      lon_rt[i] = lon_array[i]; 
      za_rt[i]  = za_array[i];
      aa_rt[i]  = aa_array[i];
      if( i < (nrt-1) )
        { l_rt[i]   = l_array[i]; }
    }

  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. If *lmax* <= 0, then only the end points
  //
  interpolate_raytracing_points( r_v, lat_v, lon_v, za_v, aa_v, lstep,
                              r_rt, lat_rt, lon_rt, za_rt, aa_rt, l_rt, lmax );
}



/*===========================================================================
  === Core functions for geometrical *ppath_step* functions
  ===========================================================================*/

//! ppath_step_geom_1d
/*! 
   Calculates 1D geometrical propagation path steps.

   This is the core function to determine 1D propagation path steps by pure
   geometrical calculations. Path points are included for crossings with the 
   grids, tangent points and points of ground intersections. In addition,
   points are included in the propgation path to ensure that the distance
   along the path between the points does not exceed the selected maximum 
   length (lmax). If lmax is <= 0, this means that no length criterion shall
   be applied.

   Note that the input variables are here compressed to only hold data for
   a 1D atmosphere. For example, z_field is z_field(:,0,0).

   For more information read the chapter on propagation paths in AUG.

   \param   ppath             Output: A Ppath structure.
   \param   p_grid            Pressure grid.
   \param   z_field           Geometrical altitudes corresponding to p_grid.
   \param   r_geoid           Geoid radius.
   \param   z_ground          Ground altitude.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-05-20
*/
void ppath_step_geom_1d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        const double&     r_geoid,
        const double&     z_ground,
        const double&     lmax )
{
  // Starting radius, zenith angle and latitude
  double r_start, lat_start, za_start;

  // Index of the pressure surface being the lower limit for the
  // grid range of interest.
  Index ip;

  // Determine the variables defined above, and make asserts of input
  ppath_start_1d( r_start, lat_start, za_start, ip,
                                   ppath, p_grid, z_field, r_geoid, z_ground );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  double ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // The path is determined by another function. Determine some variables
  // needed bý that function and call the function.
  //
  // Vars to hold found path points, path step length and coding for end face
  Vector    r_v, lat_v, za_v;
  double   lstep;
  Index     endface, tanpoint;
  //
  do_gridrange_1d( r_v, lat_v, za_v, lstep, endface, tanpoint,
                   r_start, lat_start, za_start, ppc, lmax, 
                r_geoid+z_field[ip], r_geoid+z_field[ip+1], r_geoid+z_ground );


  // Fill *ppath*
  //
  String   method;
  if( lmax < 0 )
    { method     = "1D geometrical"; }
  else
    { method     = "1D geometrical with length criterion"; }
  //
  ppath_end_1d( ppath, r_v, lat_v, za_v, lstep, z_field, r_geoid, ip, endface, 
                                                    tanpoint, method, 0, ppc );


  // Make part from a tangent point and up to the starting pressure level.
  if( endface == 0  &&  tanpoint )
    {
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      ppath_step_geom_1d( ppath2, p_grid, z_field, r_geoid, z_ground, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }
}



//! ppath_step_geom_2d
/*! 
   Calculates 2D geometrical propagation path steps.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Output: A Ppath structure.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   z_field           Geometrical altitudes
   \param   r_geoid           Geoid radii.
   \param   z_ground          Ground altitudes.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_step_geom_2d(
              Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
        const double&     lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, za_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat;

  // Radii and latitudes set by *ppath_start_2d*.
  double   lat1, lat3, r1a, r3a, r3b, r1b, rground1, rground3;

  // Determine the variables defined above and make all possible asserts
  ppath_start_2d( r_start, lat_start, za_start, ip, ilat, 
                  lat1, lat3, r1a, r3a, r3b, r1b, rground1, rground3,
                         ppath, p_grid, lat_grid, z_field, r_geoid, z_ground );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  double ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // Vars to hold found path points, path step length and coding for end face
  Vector    r_v, lat_v, za_v;
  double   lstep;
  Index     endface, tanpoint;

  do_gridcell_2d( r_v, lat_v, za_v, lstep, endface, tanpoint,
                  r_start, lat_start, za_start, ppc, lmax, lat1, lat3, 
                                      r1a, r3a, r3b, r1b, rground1, rground3 );


  // Fill *ppath*
  //
  String   method;
  if( lmax < 0 )
    { method     = "2D geometrical"; }
  else
    { method     = "2D geometrical with length criterion"; }
  //
  ppath_end_2d( ppath, r_v, lat_v, za_v, lstep, lat_grid, z_field, r_geoid, 
                                 ip, ilat, endface, tanpoint, method, 0, ppc );


  // Make part after a tangent point.
  //
  if( endface == 0  &&  tanpoint )
    {
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      // Call this function recursively
      ppath_step_geom_2d( ppath2, p_grid, lat_grid, z_field,
                                                     r_geoid, z_ground, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }
}



//! ppath_step_geom_3d
/*! 
   Calculates 3D geometrical propagation path steps.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Output: A Ppath structure.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   lon_grid          Longitude grid.
   \param   z_field           Geometrical altitudes
   \param   r_geoid           Geoid radii.
   \param   z_ground          Ground altitudes.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void ppath_step_geom_3d(
              Ppath&       ppath,
        ConstVectorView    p_grid,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        ConstTensor3View   z_field,
        ConstMatrixView    r_geoid,
        ConstMatrixView    z_ground,
        const double&      lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, lon_start, za_start, aa_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat, ilon;

  // Radius for corner points, latitude and longitude of the grid cell
  //
  double   lat1, lat3, lon5, lon6;
  double   r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b;
  double   rground15, rground35, rground36, rground16;

  // Determine the variables defined above and make all possible asserts
  ppath_start_3d( r_start, lat_start, lon_start, za_start, aa_start, 
                  ip, ilat, ilon, lat1, lat3, lon5, lon6,
                  r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b, 
                  rground15, rground35, rground36, rground16,
               ppath, p_grid, lat_grid, lon_grid, z_field, r_geoid, z_ground );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  double ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // Vars to hold found path points, path step length and coding for end face
  Vector    r_v, lat_v, lon_v, za_v, aa_v;
  double   lstep;
  Index     endface, tanpoint;

  do_gridcell_3d( r_v, lat_v, lon_v, za_v, aa_v, lstep, endface, tanpoint,
                  r_start, lat_start, lon_start, za_start, aa_start, ppc, lmax,
                  lat1, lat3, lon5, lon6, 
                  r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b,
                                  rground15, rground35, rground36, rground16 );

  // Fill *ppath*
  //
  String   method;
  if( lmax < 0 )
    { method     = "3D geometrical"; }
  else
    { method     = "3D geometrical with length criterion"; }
  //
  ppath_end_3d( ppath, r_v, lat_v, lon_v, za_v, aa_v, lstep, lat_grid, 
                lon_grid, z_field, r_geoid, ip, ilat, ilon, endface, tanpoint,
                                                              method, 0, ppc );


  // Make part after a tangent point.
  //
  if( endface == 0  &&  tanpoint )
    {
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      // Call this function recursively
      ppath_step_geom_3d( ppath2, p_grid, lat_grid, lon_grid, z_field,
                                                     r_geoid, z_ground, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }
}





/*===========================================================================
  === Ray tracing functions
  ===========================================================================*/

//! raytrace_1d_linear_euler
/*! 
   Performs ray tracing for 1D with linear Euler steps.

   A geometrical step with length of *lraytrace* is taken from each
   point. The zenith angle for the end point of that step is
   calculated exactly by the expression c = r*n*sin(theta), and a new
   step is taken. The length of the last ray tracing step to reach the
   end radius is adopted to the distance to the end radius.

   The refractive index is assumed to vary linearly between the pressure
   surfaces.

   As the ray tracing is performed from the last end point, the found path
   will not be symmetric around the tangent point.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of AUG.

   The array variables *r_array*, *lat_array* and *za_array* shall include
   the start position when calling the function. The length of *l_array*
   will be one smaller than the length of the other arrays.

   \param   r_array      Out: Radius of ray tracing points.
   \param   lat_array    Out: Latitude of ray tracing points.
   \param   za_array     Out: LOS zenith angle at ray tracing points.
   \param   l_array      Out: Distance along the path between ray tracing 
                         points.
   \param   endface      FIXME: Add documentation.
   \param   tanpoint     FIXME: Add documentation.
   \param   r            Start radius for ray tracing.
   \param   lat          Start latitude for ray tracing.
   \param   za           Start zenith angle for ray tracing.
   \param   rte_pressure     FIXME: Add documentation.
   \param   rte_temperature  FIXME: Add documentation.
   \param   rte_vmr_list     FIXME: Add documentation.
   \param   refr_index   Refractive index.
   \param   refr_index_agenda   FIXME: Add documentation.
   \param   ppc          Propagation path constant.
   \param   lraytrace    Maximum allowed length for ray tracing steps.
   \param   r1           Radius of lower pressure surface.
   \param   r3           Radius of upper pressure surface (r3 > r1).
   \param   r_ground     Radius of the ground.
   \param   r_geoid      Geoid radii.
   \param   p_grid       Pressure grid.
   \param   z_field      Geometrical altitudes corresponding to p_grid.
   \param   t_field      Temperatures corresponding to p_grid.
   \param   vmr_field    VMR values corresponding to p_grid.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void raytrace_1d_linear_euler(
              Array<double>&   r_array,
              Array<double>&   lat_array,
              Array<double>&   za_array,
              Array<double>&   l_array,
              Index&           endface,
              Index&           tanpoint,
              double           r,
              double           lat,
              double           za,
              Numeric&         rte_pressure,
              Numeric&         rte_temperature,
              Vector&          rte_vmr_list,
              Numeric&         refr_index,
        const Agenda&          refr_index_agenda,
        const double&          ppc,
        const double&          lraytrace,
        const double&          r1,
        const double&          r3,
        const double&          r_ground,
        const double&          r_geoid,
        ConstVectorView        p_grid,
        ConstVectorView        z_field,
        ConstVectorView        t_field,
        ConstMatrixView        vmr_field )
{
  // Loop boolean
  bool ready = false;

  // Variables for output from do_gridrange_1d
  Vector    r_v, lat_v, za_v;
  double   lstep, dlat = 9999;
  Index     agenda_verb = 0;

  while( !ready )
    {
      // Constant for the geometrical step to make
      const double   ppc_step = geometrical_ppc( r, za );

      // Where will a geometric path exit the grid cell?
      do_gridrange_1d( r_v, lat_v, za_v, lstep, endface, tanpoint, r, lat, za,
                                              ppc_step, -1, r1, r3, r_ground );

      assert( r_v.nelem() == 2 );

      // If *lstep* is < *lraytrace*, extract the found end point and
      // we are ready. 
      // Otherwise, we make a geometrical step with length *lraytrace*.
      // To avoid that the refrcation calculations end up outside the grid
      // cell when *lstep* is just below *lraytrace*, we allow *lstep*
      // to be somewhat larger than *lraytrace*.

      if( lstep <= 1.01 * lraytrace )
        {
          r     = r_v[1];
          dlat  = lat_v[1] - lat;
          lat   = lat_v[1];
          ready = true;
        }
      else
        {
          if( abs(za) <= 90 )
            { lstep = lraytrace; }
          else
            { lstep = -lraytrace; }

          const double   r_new = geompath_r_at_l( ppc_step, 
                                         geompath_l_at_r(ppc_step,r) + lstep );

          dlat = RAD2DEG * acos( ( r_new*r_new + r*r - 
                                             lstep*lstep ) / ( 2 * r_new*r ) );
          r     = r_new;
          lat   = lat + dlat;
          lstep = abs( lstep );
        }

      // Calculate LOS zenith angle at found point.
      //
      // Refractive index at *r*
      get_refr_index_1d( rte_pressure, rte_temperature, rte_vmr_list, 
                         refr_index_agenda, agenda_verb, p_grid, r_geoid, 
                         z_field, t_field, vmr_field, r );

      agenda_verb = 1;

      const double   ppc_local = ppc / refr_index; 

      if( ready  &&  tanpoint )
        { za = 90; }
      else if( r >= ppc_local )
        { za = geompath_za_at_r( ppc_local, za, r ); }
      else
        {  
          r = ppc_local;
          if( za > 90 )
            { 
              ready    = true;  
              tanpoint = 1;  
              if( r > r1 )
                { endface  = 0; }
            }
          za = 90;
        }
  
      // Store found point
      r_array.push_back( r );
      lat_array.push_back( lat );
      za_array.push_back( za );
      l_array.push_back( lstep );
    }  
}



//! raytrace_2d_linear_euler
/*! 
   Performs ray tracing for 2D with linear Euler steps.

   A geometrical step with length of *lraytrace* is taken from each
   point. The zenith angle for the end point of that step is
   calculated considering the gradient of the refractive index. The
   length of the last ray tracing step to reach the end radius is
   adopted to the distance to the end radius.

   The refractive index is assumed to vary linearly along the pressure
   surfaces and the latitude grid points.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of AUG.

   The array variables *r_array*, *lat_array* and *za_array* shall include
   the start position when calling the function. The length of *l_array*
   will be one smaller than the length of the other arrays.

   \param   r_array         Out: Radius of ray tracing points.
   \param   lat_array       Out: Latitude of ray tracing points.
   \param   za_array        Out: LOS zenith angle at ray tracing points.
   \param   l_array         Out: Distance along the path between ray tracing 
                            points.
   \param   endface         Out: Number coding of exit face.
   \param   tanpoint        Out: True if last point is a tangent point.
   \param   r               Out: Start radius for ray tracing.
   \param   lat             Out: Start latitude for ray tracing.
   \param   za              Out: Start zenith angle for ray tracing.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda    The WSV with the same name.
   \param   lraytrace       Maximum allowed length for ray tracing steps.
   \param   lat1            Latitude of left end face of the grid cell.
   \param   lat3            Latitude of right end face  of the grid cell.
   \param   r1a             Radius of lower-left corner of the grid cell.
   \param   r3a             Radius of lower-right corner of the grid cell.
   \param   r3b             Radius of upper-right corner of the grid cell.
   \param   r1b             Radius of upper-left corner of the grid cell.
   \param   rground1        Radius for the ground at *lat1*.
   \param   rground3        Radius for the ground at *lat3*.
   \param   p_grid          The WSV with the same name.
   \param   lat_grid        The WSV with the same name.
   \param   r_geoid         The WSV with the same name.
   \param   z_field         The WSV with the same name.
   \param   t_field         The WSV with the same name.
   \param   vmr_field       The WSV with the same name.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void raytrace_2d_linear_euler(
              Array<double>&   r_array,
              Array<double>&   lat_array,
              Array<double>&   za_array,
              Array<double>&   l_array,
              Index&           endface,
              Index&           tanpoint,
              double           r,
              double           lat,
              double           za,
              Numeric&         rte_pressure,
              Numeric&         rte_temperature,
              Vector&          rte_vmr_list,
              Numeric&         refr_index,
        const Agenda&          refr_index_agenda,
        const double&          lraytrace,
        const double&          lat1,
        const double&          lat3,
        const double&          r1a,
        const double&          r3a,
        const double&          r3b,
        const double&          r1b,
        const double&          rground1,
        const double&          rground3,
        ConstVectorView        p_grid,
        ConstVectorView        lat_grid,
        ConstVectorView        r_geoid,
        ConstMatrixView        z_field,
        ConstMatrixView        t_field,
        ConstTensor3View       vmr_field )
{
  // Loop boolean
  bool ready = false;

  // Variables for output from do_gridcell_2d
  Vector   r_v, lat_v, za_v;
  double   lstep, dlat = 9999, r_new, lat_new;
  Index    agenda_verb = 0;

  while( !ready )
    {
      // Constant for the geometrical step to make
      const double   ppc_step = geometrical_ppc( r, za );

      // Where will a geometric path exit the grid cell?
      do_gridcell_2d( r_v, lat_v, za_v, lstep, endface, tanpoint,
                     r, lat, za, ppc_step, -1, lat1, lat3, r1a, r3a, r3b, r1b, 
                                                          rground1, rground3 );
      assert( r_v.nelem() == 2 );

      // If *lstep* is < *lraytrace*, extract the found end point and
      // we are ready.
      // Otherwise, we make a geometrical step with length *lraytrace*.
      // To avoid that the refrcation calculations end up outside the grid
      // cell when *lstep* is just below *lraytrace*, we allow *lstep*
      // to be somewhat larger than *lraytrace*.

      if( lstep <= 1.01 * lraytrace )
        {
          r_new   = r_v[1];
          dlat    = lat_v[1] - lat;
          lat_new = lat_v[1];
          ready   = true;
        }
      else
        {
          if( abs(za) <= 90 )
            { lstep = lraytrace; }
          else
            { lstep = -lraytrace; }

          r_new = geompath_r_at_l( ppc_step, 
                                         geompath_l_at_r(ppc_step,r) + lstep );
          dlat = RAD2DEG * acos( ( r_new*r_new + r*r - 
                                             lstep*lstep ) / ( 2 * r_new*r ) );
          if( za < 0 )
            { dlat = -dlat; }

          lat_new = lat + dlat;
          lstep   = abs( lstep );

          // For paths along the latitude end faces we can end up outside the
          // grid cell. We simply look for points outisde the grid cell.
          if( lat_new < lat1 )
            { lat_new = lat1; }
          else if( lat_new > lat3 )
            { lat_new = lat3; }
        }

      // Calculate LOS zenith angle at found point.
      if( ready  &&  tanpoint )
        { 
          // It is not totally correct to set *za* to 90 here. We neglect
          // then the curvature of this ray tracing step, but we don't care
          // about this for simplicity. This point has been flagged as the
          // the tangent point and we treat it then it that way.
          za = sign(za) * 90; 
        }
      else
        {
                Numeric   dndr, dndlat;
          const double   za_rad = DEG2RAD * za;

          refr_gradients_2d( refr_index, dndr, dndlat, rte_pressure, 
                             rte_temperature, rte_vmr_list, refr_index_agenda,
                             agenda_verb,
                             p_grid, lat_grid, r_geoid, z_field, t_field, 
                             vmr_field, r, lat );

          agenda_verb = 1;

          za += -dlat + RAD2DEG * lstep / refr_index * ( -sin(za_rad) * dndr +
                                                        cos(za_rad) * dndlat );

          // Make sure that obtained *za* is inside valid range
          if( za < -180 )
            { za += 360; }
          else if( za > 180 )
            { za -= 360; }
        }

      r   = r_new;
      lat = lat_new;

      // If the path is zenith/nadir along a latitude end face, we must check
      // that the path does not exit with new *za*.
      if( lat == lat1  &&  za < 0 )
        { endface = 1;   ready = 1; }
      else if( lat == lat3  &&  za > 0 )
        { endface = 3;   ready = 1; }
      
      // Store found point
      r_array.push_back( r );
      lat_array.push_back( lat );
      za_array.push_back( za );
      l_array.push_back( lstep );
    }  
}



//! raytrace_3d_linear_euler
/*! 
   Performs ray tracing for 3D with linear Euler steps.

   A geometrical step with length of *lraytrace* is taken from each
   point. The zenith angle for the end point of that step is
   calculated considering the gradient of the refractive index. The
   length of the last ray tracing step to reach the end radius is
   adopted to the distance to the end radius.

   The refractive index is assumed to vary linearly along the pressure
   surfaces and the latitude grid points.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of AUG.

   The array variables *r_array*, *lat_array*, *lon_array*, *za_array*
   and *aa_array* shall include the start position when calling the
   function. The length of *l_array* will be one smaller than the
   length of the other arrays.

   \param   r_array        Out: Radius of ray tracing points.
   \param   lat_array      Out: Latitude of ray tracing points.
   \param   lon_array      Out: Longitude of ray tracing points.
   \param   za_array       Out: LOS zenith angle at ray tracing points.
   \param   aa_array       Out: LOS azimuth angle at ray tracing points.
   \param   l_array        Out: Distance along the path between ray tracing 
                           points.
   \param   endface        Out: Number coding of exit face.
   \param   tanpoint       Out: True if last point is a tangent point.
   \param   r              Out: Start radius for ray tracing.
   \param   lat            Out: Start latitude for ray tracing.
   \param   lon            Out: Start longitude for ray tracing.
   \param   za             Out: Start zenith angle for ray tracing.
   \param   aa             Out: Start azimuth angle for ray tracing.
   \param   rte_pressure     Out: The WSV with the same name.
   \param   rte_temperature  Out: The WSV with the same name.
   \param   rte_vmr_list     Out: The WSV with the same name.
   \param   refr_index       Out: The WSV with the same name.
   \param   refr_index_agenda    The WSV with the same name.
   \param   lraytrace      Maximum allowed length for ray tracing steps.
   \param   lat1           Latitude of left end face of the grid cell.
   \param   lat3           Latitude of right end face of the grid cell.
   \param   lon5           Lower longitude of the grid cell.
   \param   lon6           Upper longitude of the grid cell.
   \param   r15a           Radius of corner: lower p-surface,*lat1* and *lon5*.
   \param   r35a           Radius of corner: lower p-surface,*lat3* and *lon5*.
   \param   r36a           Radius of corner: lower p-surface,*lat3* and *lon6*.
   \param   r16a           Radius of corner: lower p-surface,*lat1* and *lon6*.
   \param   r15b           Radius of corner: upper p-surface,*lat1* and *lon5*.
   \param   r35b           Radius of corner: upper p-surface,*lat3* and *lon5*.
   \param   r36b           Radius of corner: upper p-surface,*lat3* and *lon6*.
   \param   r16b           Radius of corner: upper p-surface,*lat1* and *lon6*.
   \param   rground15      Radius for the ground at *lat1* and *lon5*.
   \param   rground35      Radius for the ground at *lat3* and *lon5*.
   \param   rground36      Radius for the ground at *lat3* and *lon6*.
   \param   rground16      Radius for the ground at *lat1* and *lon6*.   
   \param   p_grid         The WSV with the same name.
   \param   lat_grid       The WSV with the same name.
   \param   lon_grid       The WSV with the same name.
   \param   r_geoid        The WSV with the same name.
   \param   z_field        The WSV with the same name.
   \param   t_field        The WSV with the same name.
   \param   vmr_field      The WSV with the same name.

   \author Patrick Eriksson
   \date   2003-01-18
*/
void raytrace_3d_linear_euler(
              Array<double>&   r_array,
              Array<double>&   lat_array,
              Array<double>&   lon_array,
              Array<double>&   za_array,
              Array<double>&   aa_array,
              Array<double>&   l_array,
              Index&           endface,
              Index&           tanpoint,
              double           r,
              double           lat,
              double           lon,
              double           za,
              double           aa,
              Numeric&         rte_pressure,
              Numeric&         rte_temperature,
              Vector&          rte_vmr_list,
              Numeric&         refr_index,
        const Agenda&          refr_index_agenda,
        const double&          lraytrace,
        const double&          lat1,
        const double&          lat3,
        const double&          lon5,
        const double&          lon6,
        const double&          r15a,
        const double&          r35a,
        const double&          r36a,
        const double&          r16a,
        const double&          r15b,
        const double&          r35b,
        const double&          r36b,
        const double&          r16b,
        const double&          rground15,
        const double&          rground35,
        const double&          rground36,
        const double&          rground16,
        ConstVectorView        p_grid,
        ConstVectorView        lat_grid,
        ConstVectorView        lon_grid,
        ConstMatrixView        r_geoid,
        ConstTensor3View       z_field,
        ConstTensor3View       t_field,
        ConstTensor4View       vmr_field )
{
  // Loop boolean
  bool ready = false;

  // Variables for output from do_gridcell_2d
  Vector    r_v, lat_v, lon_v, za_v, aa_v;
  double    r_new, lat_new, lon_new, za_new, aa_new;
  double    lstep;
  Index     agenda_verb = 0;

  while( !ready )
    {
      // Constant for the geometrical step to make
      const double   ppc_step = geometrical_ppc( r, za );

      // Where will a geometric path exit the grid cell?
      do_gridcell_3d( r_v, lat_v, lon_v, za_v, aa_v, lstep, endface, tanpoint,
                    r, lat, lon, za, aa, ppc_step, -1, lat1, lat3, lon5, lon6, 
                    r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b,
                                  rground15, rground35, rground36, rground16 );

      assert( r_v.nelem() == 2 );


      // If *lstep* is < *lraytrace*, extract the found end point and
      // we are ready.
      // Otherwise, we make a geometrical step with length *lraytrace*.
      // To avoid that the refrcation calculations end up outside the grid
      // cell when *lstep* is just below *lraytrace*, we allow *lstep*
      // to be somewhat larger than *lraytrace*.

      if( lstep <= 1.01 * lraytrace )
        {
          r_new   = r_v[1];
          lat_new = lat_v[1];
          lon_new = lon_v[1];
          za_new  = za_v[1];
          aa_new  = aa_v[1];
          ready   = true;
        }
      else
        {
          // Sensor pos and LOS in cartesian coordinates
          double   x, y, z, dx, dy, dz;
          poslos2cart( x, y, z, dx, dy, dz, r, lat, lon, za, aa ); 

          lstep = lraytrace;

          cart2poslos( r_new, lat_new, lon_new, za_new, aa_new, 
                              x+dx*lstep, y+dy*lstep, z+dz*lstep, dx, dy, dz );

          // For paths along some end face we can end up outside the
          // grid cell. We simply look for points outisde the grid cell.
          if( lat_new < lat1 )
            { lat_new = lat1; }
          else if( lat_new > lat3 )
            { lat_new = lat3; }
          if( lon_new < lon5 )
            { lon_new = lon5; }
          else if( lon_new > lon6 )
            { lon_new = lon6; }

          // Checks to improve the accuracy for speciel cases.
          // The checks are only needed for values cacluated locally as
          // the same checks are made inside *do_gridcell_3d*.

          //--- Set zenith angle to be as accurate as possible
          if( za < ANGTOL  ||  za > 180-ANGTOL )
            { za_new = za; }
          else
            { za_new = geompath_za_at_r( ppc_step, za, r_new ); }

          //--- Set azimuth angle and lon. to be as accurate as possible for
          //    north-south observations
          if( abs( lat ) < 90  &&  
                     ( abs(aa) < ANGTOL  ||  abs( aa ) > 180-ANGTOL ) )
            { 
              aa_new  = aa; 
              lon_new = lon;
            }

          // Shall lon values be shifted?
          resolve_lon( lon_new, lon5, lon6 );
        }

      // Calculate LOS zenith angle at found point.
      if( ready  &&  tanpoint )
        {
          // It is not totally correct to set *za* to 90 here. We neglect
          // then the curvature of this ray tracing step, but we don't care
          // about this for simplicity. This point has been flagged as the
          // the tangent point and we treat it then it that way.
          za_new = 90; 
        }
      else
        {
          Numeric   dndr, dndlat, dndlon;

          refr_gradients_3d( refr_index, dndr, dndlat, dndlon, rte_pressure, 
                             rte_temperature, rte_vmr_list, refr_index_agenda, 
                             agenda_verb,
                             p_grid, lat_grid, lon_grid, r_geoid, z_field, 
                             t_field, vmr_field, r, lat, lon );

          agenda_verb = 1;

          const double   aterm = RAD2DEG * lstep / refr_index;
          const double   za_rad = DEG2RAD * za;
          const double   aa_rad = DEG2RAD * aa;
          const double   sinza = sin( za_rad );
          const double   sinaa = sin( aa_rad );
          const double   cosaa = cos( aa_rad );


          if( za < ANGTOL  ||  za > 180-ANGTOL )
            { 
              za_new += aterm * ( cos(za_rad) * 
                                       ( cosaa * dndlat + sinaa * dndlon ) );
              aa_new = RAD2DEG * atan2( dndlon, dndlat); 
            }
          else
            { 
              za_new += aterm * ( -sinza * dndr + cos(za_rad) * 
                                       ( cosaa * dndlat + sinaa * dndlon ) );
              aa_new += aterm * sinza * ( cosaa * dndlon - sinaa * dndlat ); 
            }
          
          // Make sure that obtained angles are inside valid ranges
          if( za_new > 180 )
            { za_new = 180 - za_new; }
          else if( za_new < 0 )
            { za_new = -za_new; }
          //
          if( aa_new > 180 )
            { aa_new -= 360; }
          else if( aa_new < -180 )
            { aa_new += 360; }
        }

      r   = r_new;
      lat = lat_new;
      lon = lon_new;
      za  = za_new;
      aa  = aa_new;

      // For some cases where the path goes along an end face, 
      // it could be the case that the refraction bends the path out 
      // of the grid cell.
      if( lon == lon5  &&  aa < 0 )
        { endface = 5;   ready = 1; }
      else if( lon == lon6  &&  aa > 0 )
        { endface = 6;   ready = 1; }
      else if( za > 0  &&  za < 180 )
        {
          if( lat == lat1  &&  lat != -90  &&  abs( aa ) > 90 ) 
            { endface = 1;   ready = 1; }
          else if( lat == lat3  &&  lat != 90  &&  abs( aa ) < 90 ) 
            { endface = 3;   ready = 1; }
        }

      // Store found point
      r_array.push_back( r );
      lat_array.push_back( lat );
      lon_array.push_back( lon );
      za_array.push_back( za );
      aa_array.push_back( aa );
      l_array.push_back( lstep );
    }  
}





/*===========================================================================
  === Core functions for refraction *ppath_step* functions
  ===========================================================================*/

//! ppath_step_refr_1d
/*! 
   Calculates 1D propagation path steps including effects of refraction.

   This function works as the function *ppath_step_geom_1d* but considers
   also refraction. The upper length of the ray tracing steps is set by
   the argument *lraytrace*. This argument controls only the internal
   calculations. The maximum distance between the path points is still
   determined by *lmax*.

   \param   ppath             Out: A Ppath structure.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda The WSV with the same name.
   \param   p_grid            Pressure grid.
   \param   z_field           Geometrical altitudes corresponding to p_grid.
   \param   t_field           Temperatures corresponding to p_grid.
   \param   vmr_field         VMR values corresponding to p_grid.
   \param   r_geoid           Geoid radius.
   \param   z_ground          Ground altitude.
   \param   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-11-26
*/
void ppath_step_refr_1d(
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
        const double&     r_geoid,
        const double&     z_ground,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax )
{
  // Starting radius, zenith angle and latitude
  double   r_start, lat_start, za_start;

  // Index of the pressure surface being the lower limit for the
  // grid range of interest.
  Index   ip;

  // Determine the variables defined above, and make asserts of input
  ppath_start_1d( r_start, lat_start, za_start, ip, 
                                   ppath, p_grid, z_field, r_geoid, z_ground );

  // Assert not done for geometrical calculations
  assert( t_field.nelem() == p_grid.nelem() );
  assert( vmr_field.ncols() == p_grid.nelem() );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  // If the sensor is placed outside the atmosphere, the constant is
  // already set.
  double ppc;
  if( ppath.constant < 0 )
    { 
      get_refr_index_1d( rte_pressure, rte_temperature, rte_vmr_list, 
                         refr_index_agenda, 1, p_grid, r_geoid, z_field, 
                         t_field, vmr_field, r_start );
      ppc = refraction_ppc( r_start, za_start, refr_index ); 
    }
  else
    { ppc = ppath.constant; }


  // Perform the ray tracing
  //
  // Arrays to store found ray tracing points
  // (Vectors don't work here as we don't know how many points there will be)
  Array<double>   r_array, lat_array, za_array, l_array;
  //
  // Store the start point
  r_array.push_back( r_start );
  lat_array.push_back( lat_start );
  za_array.push_back( za_start );
  //
  // String to store description of ray tracing method
  String   method;
  //
  // Number coding for end face
  Index   endface, tanpoint;
  //
  if( rtrace_method  == "linear_euler" )
    {
      if( lmax < 0 )
        { method = "1D linear Euler"; }
      else
        { method = "1D linear Euler, with length criterion"; }

      raytrace_1d_linear_euler( r_array, lat_array, za_array, l_array, endface,
        tanpoint, r_start, lat_start, za_start, rte_pressure, rte_temperature, 
            rte_vmr_list, refr_index, refr_index_agenda, ppc, lraytrace, 
            r_geoid+z_field[ip], r_geoid+z_field[ip+1], r_geoid + z_ground, 
                                r_geoid, p_grid, z_field, t_field, vmr_field );
    }
  else
    {
      bool   known_ray_trace_method = false;
      assert( known_ray_trace_method );
    }


  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. 
  //
  Vector    r_v, lat_v, za_v;
  double   lstep;
  //
  from_raytracingarrays_to_ppath_vectors_1d_and_2d( r_v, lat_v, za_v, lstep, 
                              r_array, lat_array, za_array, l_array, 0, lmax );


  // Fill *ppath*
  ppath_end_1d( ppath, r_v, lat_v, za_v, lstep, z_field, r_geoid, ip, endface, 
                                                    tanpoint, method, 1, ppc );


  //--- End point is a tangent point
  if( endface == 0  &&  tanpoint )
    {
      // Make part from tangent point and up to the starting pressure level.
      //
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      ppath_step_refr_1d( ppath2, rte_pressure, rte_temperature, rte_vmr_list,
                      refr_index, refr_index_agenda, p_grid, z_field, t_field,
                vmr_field, r_geoid, z_ground, rtrace_method, lraytrace, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }
}



//! ppath_step_refr_2d
/*! 
   Calculates 2D propagation path steps, with refraction, using a simple
   and fast ray tracing scheme.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Out: A Ppath structure.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda The WSV with the same name.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   z_field           Geometrical altitudes.
   \param   t_field           Atmospheric temperatures.
   \param   vmr_field         VMR values.
   \param   r_geoid           Geoid radii.
   \param   z_ground          Ground altitudes.
   \param   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void ppath_step_refr_2d(
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
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, za_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat;

  // Radii and latitudes set by *ppath_start_2d*.
  double   lat1, lat3, r1a, r3a, r3b, r1b, rground1, rground3;

  // Determine the variables defined above and make all possible asserts
  ppath_start_2d( r_start, lat_start, za_start, ip, ilat, 
                  lat1, lat3, r1a, r3a, r3b, r1b, rground1, rground3,
                         ppath, p_grid, lat_grid, z_field, r_geoid, z_ground );

  // Assert not done for geometrical calculations
  assert( t_field.nrows() == p_grid.nelem() );
  assert( t_field.ncols() == lat_grid.nelem() );
  assert( vmr_field.nrows() == p_grid.nelem() );
  assert( vmr_field.ncols() == lat_grid.nelem() );

  // No constant for the path is valid here.


  // Perform the ray tracing
  //
  // Arrays to store found ray tracing points
  // (Vectors don't work here as we don't know how many points there will be)
  Array<double>   r_array, lat_array, za_array, l_array;
  //
  // Store the start point
  r_array.push_back( r_start );
  lat_array.push_back( lat_start );
  za_array.push_back( za_start );
  //
  // String to store description of ray tracing method
  String   method;
  //
  // Number coding for end face
  Index   endface, tanpoint;
  //
  if( rtrace_method  == "linear_euler" )
    {
      if( lmax < 0 )
        { method = "2D linear Euler"; }
      else
        { method = "2D linear Euler, with length criterion"; }

      raytrace_2d_linear_euler( r_array, lat_array, za_array, l_array, endface,
        tanpoint, r_start, lat_start, za_start, rte_pressure, rte_temperature, 
                    rte_vmr_list, refr_index, refr_index_agenda, lraytrace, 
                          lat1, lat3, r1a, r3a, r3b, r1b, rground1, rground3,
                      p_grid, lat_grid, r_geoid, z_field, t_field, vmr_field );
    }
  else
    {
      bool   known_ray_trace_method = false;
      assert( known_ray_trace_method );
    }


  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. 
  //
  Vector    r_v, lat_v, za_v;
  double   lstep;
  //
  from_raytracingarrays_to_ppath_vectors_1d_and_2d( r_v, lat_v, za_v, lstep, 
                              r_array, lat_array, za_array, l_array, 0, lmax );


  // Fill *ppath*
  ppath_end_2d( ppath, r_v, lat_v, za_v, lstep, lat_grid, z_field, r_geoid, 
                                  ip, ilat, endface, tanpoint, method, 1, -1 );


  // Make part after a tangent point.
  //
  if( endface == 0  &&  tanpoint )
    {
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      // Call this function recursively
      ppath_step_refr_2d( ppath2, rte_pressure, rte_temperature, rte_vmr_list,
                    refr_index, refr_index_agenda, p_grid, lat_grid, z_field, 
                    t_field, vmr_field, 
                           r_geoid, z_ground, rtrace_method, lraytrace, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }
}



//! ppath_step_refr_3d
/*! 
   Calculates 3D propagation path steps, with refraction, using a simple
   and fast ray tracing scheme.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Out: A Ppath structure.
   \param   rte_pressure      Out: The WSV with the same name.
   \param   rte_temperature   Out: The WSV with the same name.
   \param   rte_vmr_list      Out: The WSV with the same name.
   \param   refr_index        Out: The WSV with the same name.
   \param   refr_index_agenda The WSV with the same name.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   lon_grid          Longitude grid.
   \param   z_field           Geometrical altitudes.
   \param   t_field           Atmospheric temperatures.
   \param   vmr_field         VMR values.
   \param   r_geoid           Geoid radii.
   \param   z_ground          Ground altitudes.
   \param   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2003-01-08
*/
void ppath_step_refr_3d(
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
        ConstMatrixView   r_geoid,
        ConstMatrixView   z_ground,
        const String&     rtrace_method,
        const double&     lraytrace,
        const double&     lmax )
{
  // Radius, zenith angle and latitude of start point.
  double   r_start, lat_start, lon_start, za_start, aa_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat, ilon;

  // Radius for corner points, latitude and longitude of the grid cell
  //
  double   lat1, lat3, lon5, lon6;
  double   r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b;
  double   rground15, rground35, rground36, rground16;

  // Determine the variables defined above and make all possible asserts
  ppath_start_3d( r_start, lat_start, lon_start, za_start, aa_start, 
                  ip, ilat, ilon, lat1, lat3, lon5, lon6,
                  r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b, 
                  rground15, rground35, rground36, rground16,
               ppath, p_grid, lat_grid, lon_grid, z_field, r_geoid, z_ground );

  // Assert not done for geometrical calculations
  assert( t_field.npages() == p_grid.nelem() );
  assert( t_field.nrows() == lat_grid.nelem() );
  assert( t_field.ncols() == lon_grid.nelem() );
  assert( vmr_field.npages() == p_grid.nelem() );
  assert( vmr_field.nrows() == lat_grid.nelem() );
  assert( vmr_field.ncols() == lon_grid.nelem() );

  // No constant for the path is valid here.


  // Perform the ray tracing
  //
  // Arrays to store found ray tracing points
  // (Vectors don't work here as we don't know how many points there will be)
  Array<double>   r_array, lat_array, lon_array, za_array, aa_array, l_array;
  //
  // Store the start point
  r_array.push_back( r_start );
  lat_array.push_back( lat_start );
  lon_array.push_back( lon_start );
  za_array.push_back( za_start );
  aa_array.push_back( aa_start );
  //
  // String to store description of ray tracing method
  String   method;
  //
  // Number coding for end face
  Index   endface, tanpoint;
  //
  if( rtrace_method  == "linear_euler" )
    {
      if( lmax < 0 )
        { method = "3D linear Euler"; }
      else
        { method = "3D linear Euler, with length criterion"; }

      raytrace_3d_linear_euler( r_array, lat_array, lon_array, za_array, 
                                aa_array, l_array, endface, tanpoint, r_start,
                                lat_start, lon_start, za_start, aa_start, 
                                rte_pressure, rte_temperature, rte_vmr_list, 
                                refr_index, refr_index_agenda, lraytrace, 
                                lat1, lat3, lon5, lon6, 
                                r15a, r35a, r36a, r15a, r15b, r35b, r36b, r15b,
                                rground15, rground35, rground36, rground16,
                                p_grid, lat_grid, lon_grid, 
                                r_geoid, z_field, t_field, vmr_field );
    }
  else
    {
      bool   known_ray_trace_method = false;
      assert( known_ray_trace_method );
    }


  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. 
  //
  Vector    r_v, lat_v, lon_v, za_v, aa_v;
  double   lstep;
  //
  from_raytracingarrays_to_ppath_vectors_3d( r_v, lat_v, lon_v, za_v, aa_v, 
                                        lstep, r_array, lat_array, lon_array, 
                                           za_array, aa_array, l_array, lmax );


  // Fill *ppath*
  ppath_end_3d( ppath, r_v, lat_v, lon_v, za_v, aa_v, lstep, lat_grid, 
                lon_grid, z_field, r_geoid, ip, ilat, ilon, endface, tanpoint,
                                                               method, 1, -1 );


  // Make part after a tangent point.
  //
  if( endface == 0  &&  tanpoint )
    {
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      // Call this function recursively
      ppath_step_refr_3d( ppath2, rte_pressure, rte_temperature, rte_vmr_list,
                          refr_index, refr_index_agenda, p_grid, lat_grid, 
                          lon_grid, z_field, t_field, vmr_field, 
                          r_geoid, z_ground, rtrace_method, lraytrace, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }
}




/*===========================================================================
  === Core of ppathCalc and help function(s)
  ===========================================================================*/

//! ppath_start_stepping
/*!
   Initiates a Ppath structure for calculation of a path with *ppath_step*.

   The function performs two main tasks. As mentioned above, it initiates
   a Ppath structure (a), but it also checks that the end point of the path is
   at an allowed location (b).

   (a): The Ppath structure is set to hold the position and LOS of the last
   point of the path inside the atmosphere. This point is either the
   sensor position, or the point where the path leaves the model atmosphere.
   If the path is totally outside the atmosphere, no point is put into the
   structure. If the (practical) end and start points are identical, such
   as when the sensor is inside the cloud box, the background field is set.

   (b): If it is found that the end point of the path is at an illegal position
   a detailed error message is given. Not allowed cases are: <br>  
      1. The sensor is placed below ground level. <br> 
      2. For 2D and 3D, the path leaves the model atmosphere at a latitude or
         longitude end face. <br> 
      3. For 2D and 3D, the path is totally outside the atmosphere and the 
         latitude and longitude of the tangent point is outside the range of
         the corresponding grids. 

   All input variables are identical with the WSV with the same name.
   The output variable is here called ppath for simplicity, but is in
   fact *ppath_step*.

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   lat_grid          The latitude grid.
   \param   lon_grid          The longitude grid.
   \param   z_field           The field of geometrical altitudes.
   \param   r_geoid           The geoid radius.
   \param   z_ground          Ground altitude.
   \param   cloudbox_on       Flag to activate the cloud box.
   \param   cloudbox_limits   Index limits of the cloud box.
   \param   outside_cloudbox  Boolean to flag if the propagation path is 
                              (expected to be) outside the cloudbox. Ordinary
                              clerar sky calculations are selected by the value
                              1. The value 0 means tracking of a propagation 
                              path inside the cloudbox. The path is then
                              tracked to the cloudbox boundary.
   \param   rte_pos           The position of the sensor.
   \param   rte_los           The line-of-sight of the sensor.

   \author Patrick Eriksson
   \date   2002-05-17
*/
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
        const bool&           outside_cloudbox,
        ConstVectorView       rte_pos,
        ConstVectorView       rte_los )
{
  // This function contains no checks or asserts as it is only a sub-function
  // to ppathCalc where the input is checked carefully.

  // Allocate the ppath structure
  ppath_init_structure(  ppath, atmosphere_dim, 1 );

  // Number of pressure levels
  const Index np = p_grid.nelem();

  // The different atmospheric dimensionalities are handled seperately

  //-- 1D ---------------------------------------------------------------------
  if( atmosphere_dim == 1 )
    {
      // Radius for the ground
      const double r_ground = r_geoid(0,0) + z_ground(0,0);

      // Radius for the top of the atmosphere
      const double r_top = r_geoid(0,0) + z_field(np-1,0,0);

      // The only forbidden case here is that the sensor is below the ground
      if( rte_pos[0] < r_ground )
        {
          ostringstream os;
          os << "The sensor is placed " 
             << (r_ground - rte_pos[0])/1e3 << " km below ground level.\n"
             << "The sensor must be above the ground.";
          throw runtime_error(os.str());
        }

      out2 << "  sensor altitude        : " << (rte_pos[0]-r_geoid(0,0))/1e3 
           << " km\n";

      // If downwards, calculate geometrical tangent position
      Vector geom_tan_pos(0);
      if( rte_los[0] > 90 )
        {
          geom_tan_pos.resize(2);
          geom_tan_pos[0] = geometrical_ppc( rte_pos[0], rte_los[0] );
          geom_tan_pos[1] = geompath_lat_at_za( rte_los[0], 0, 90 );
          out2 << "  geom. tangent radius   : " << geom_tan_pos[0]/1e3 
               <<" km\n";
          out2 << "  geom. tangent latitude : " << geom_tan_pos[1] << "\n";
          out2 << "  geom. tangent altitude : " 
               << (geom_tan_pos[0]-r_geoid(0,0))/1e3 << " km\n";
        }

      // Put sensor position and LOS in ppath as first guess
      ppath.pos(0,0) = rte_pos[0];
      ppath.los(0,0) = rte_los[0];
      ppath.pos(0,1) = 0; 
      ppath.z[0]     = ppath.pos(0,0) - r_geoid(0,0);
      

      // The sensor is inside the model atmosphere, 1D ------------------------
      if( rte_pos[0] < r_top )
        {
          // Use below the values in ppath (instead of rte_pos and rte_los) as 
          // they can be modified on the way.
     
          // Is the sensor on the ground looking down?
          // If yes and the sensor is inside the cloudbox, the background will
          // be changed below.
          if( ppath.pos(0,0) == r_ground  &&  ppath.los(0,0) > 90 )
            { ppath_set_background( ppath, 2 ); }

          // Check sensor position with respect to cloud box.
          if( cloudbox_on )
            {
              if( outside_cloudbox )
                {
                  // Is the sensor inside the cloud box?
                  if( ppath.z[0] > z_field(cloudbox_limits[0],0,0)  && 
                                 ppath.z[0] < z_field(cloudbox_limits[1],0,0) )
                    { ppath_set_background( ppath, 4 ); }

                  else if( ( ppath.z[0] == z_field(cloudbox_limits[0],0,0)  && 
                                                        ppath.los(0,0) <= 90 )
                           || 
                           ( ppath.z[0] == z_field(cloudbox_limits[1],0,0)  && 
                                                        ppath.los(0,0) > 90 ) )
                    {
                      ppath_set_background( ppath, 3 );
                    }
                }
              else
                {
                  assert( ppath.z[0] >= z_field(cloudbox_limits[0],0,0)  && 
                          ppath.z[0] <= z_field(cloudbox_limits[1],0,0) );
                  
                }
            }
        }

      // The sensor is outside the model atmosphere, 1D -----------------------
      else
        {
          // Upward observations are not allowed here
          if( rte_los[0] <= 90 )
              throw runtime_error("When the sensor is placed outside the model"
                         " atmosphere, upward observations are not allowed." );

          // We can here set the path constant, that equals the radius of the
          // geometrical tangent point.
          ppath.constant = geom_tan_pos[0];
 
          // Path is above the atmosphere
          if( ppath.constant >= r_top )
            {
              ppath_set_background( ppath, 1 );
              out1 << "  --- WARNING ---, path is totally outside of the "
                   << "model atmosphere\n";
            }

          // Path enters the atmosphere
          else
            {
              ppath.z[0]     = z_field(np-1,0,0);
              ppath.pos(0,0) = r_top;
              ppath.los(0,0) = geompath_za_at_r( ppath.constant, rte_los[0], 
                                                                       r_top );
              ppath.pos(0,1) = geompath_lat_at_za( rte_los[0], 0, 
                                                              ppath.los(0,0) );
            }
        }

      // Get grid position for the end point, if it is inside the atmosphere.
      if( ppath.z[0] <= z_field(np-1,0,0) )
        { gridpos( ppath.gp_p, z_field(joker,0,0), ppath.z ); }

      // Handle possible numerical problems for grid positions
      gridpos_check_fd( ppath.gp_p[0] );

      // Set geometrical tangent point position
      if( geom_tan_pos.nelem() == 2 )
        {
          ppath.geom_tan_pos.resize(2);
          ppath.geom_tan_pos = geom_tan_pos;
        }

    }  // End 1D


  //-- 2D ---------------------------------------------------------------------
  else if( atmosphere_dim == 2 )
    {
      // Number of points in the latitude grids
      const Index   nlat = lat_grid.nelem();

      // Is the sensor inside the latitude range of the model atmosphere,
      // and below the top of the atmosphere? If yes, is_inside = true.
      // Store geoid and ground radii, grid position and interpolation weights
      // for later use.
      //
      bool      is_inside = false;   
      double    rv_geoid=-1, rv_ground=-1;  // -1 to avoid compiler warnings
      GridPos   gp_lat;
      Vector    itw(2);
      //
      if( rte_pos[1] >= lat_grid[0]  &&  rte_pos[1] <= lat_grid[nlat-1] )
        {
          gridpos( gp_lat, lat_grid, rte_pos[1] );
          interpweights( itw, gp_lat );

          rv_geoid  = interp( itw, r_geoid(joker,0), gp_lat );
          rv_ground = rv_geoid + interp( itw, z_ground(joker,0), gp_lat );

          out2 << "  sensor altitude        : " << (rte_pos[0]-rv_geoid)/1e3 
               << " km\n";

          if( rte_pos[0] < ( rv_geoid + 
                               interp( itw, z_field(np-1,joker,0), gp_lat ) ) )
            { is_inside = true; }
        }

      // If downwards, calculate geometrical tangent position. If the tangent
      // point is inside the covered latitude range, calculate also the 
      // geometrical altitude of the tangent point and the top of atmosphere.
      //
      Vector  geom_tan_pos(0);
      double geom_tan_z=-2, geom_tan_atmtop=-1;  // OK values if the variables
      //                                             are not set below
      if( abs(rte_los[0]) > 90 )
        {
          geom_tan_pos.resize(2);
          geom_tan_pos[0] = geometrical_ppc( rte_pos[0], rte_los[0] );
          if( rte_los[0] > 0 )
            { geom_tan_pos[1] = 
                            geompath_lat_at_za( rte_los[0], rte_pos[1], 90 ); }
          else
            { geom_tan_pos[1] = 
                           geompath_lat_at_za( rte_los[0], rte_pos[1], -90 ); }
          out2 << "  geom. tangent radius   : " << geom_tan_pos[0] / 1e3
               <<" km\n";
          out2 << "  geom. tangent latitude : " << geom_tan_pos[1] 
               << "\n";
          //
          if( geom_tan_pos[1] >= lat_grid[0]  &&  
                                          geom_tan_pos[1] <= lat_grid[nlat-1] )
            {
              GridPos   gp_tan;
              Vector    itw_tan(2);
              gridpos( gp_tan, lat_grid, geom_tan_pos[1] );
              interpweights( itw_tan, gp_tan );
              geom_tan_z = geom_tan_pos[0] - 
                                   interp( itw_tan, r_geoid(joker,0), gp_tan );
              geom_tan_atmtop = 
                              interp( itw_tan, z_field(np-1,joker,0), gp_tan );
              out2 << "  geom. tangent altitude : " << geom_tan_z/1e3 
                   << " km\n";
            }
        }

      // Put sensor position and LOS in ppath as first guess
      ppath.pos(0,0) = rte_pos[0];
      ppath.pos(0,1) = rte_pos[1];
      ppath.los(0,0) = rte_los[0];

      // The sensor is inside the model atmosphere, 2D ------------------------
      if( is_inside )
        {
          // Check that the sensor is above the ground
          if( rte_pos[0] < rv_ground )
            {
              ostringstream os;
              os << "The sensor is placed " 
                 << (rv_ground - rte_pos[0])/1e3 << " km below ground level.\n"
                 << "The sensor must be above the ground.";
              throw runtime_error(os.str());
            }

          // Check that not at latitude end point and looks out
          if( ( rte_pos[1] == lat_grid[0]  &&   rte_los[0] < 0 ) )
            throw runtime_error( "The sensor is at the lower latitude end "
                                           "point and the zenith angle < 0." );
          if( rte_pos[1] == lat_grid[nlat-1]  &&   rte_los[0] >= 0 ) 
            throw runtime_error( "The sensor is at the upper latitude end "
                                          "point and the zenith angle >= 0." );
          
          // Geometrical altitude
          ppath.z[0] = ppath.pos(0,0) - rv_geoid;

          // Use below the values in ppath (instead of rte_pos and rte_los) as 
          // they can be modified on the way.
     
          // Grid positions
          ppath.gp_lat[0].idx   = gp_lat.idx;
          ppath.gp_lat[0].fd[0] = gp_lat.fd[0];
          ppath.gp_lat[0].fd[1] = gp_lat.fd[1];

          // Create a vector with the geometrical altitude of the pressure 
          // surfaces for the sensor latitude and use it to get ppath.gp_p.
          Vector z_grid(np);
          z_at_lat_2d( z_grid, p_grid, lat_grid, 
                                              z_field(joker,joker,0), gp_lat );
          gridpos( ppath.gp_p, z_grid, ppath.z );

          // Is the sensor on the ground looking down?
          if( ppath.pos(0,0) == rv_ground )
            {
              // Calculate radial slope of the ground
              const double rslope = psurface_slope_2d( lat_grid, 
                        r_geoid(joker,0), z_ground(joker,0), 
                                                      gp_lat, ppath.los(0,0) );

              // Calculate angular tilt of the ground
              const double atilt = psurface_angletilt( rv_ground, rslope);

              // Are we looking down into the ground?
              // If yes and the sensor is inside the cloudbox, the background 
              // will be changed below.
              if( is_los_downwards( ppath.los(0,0), atilt ) )
                { ppath_set_background( ppath, 2 ); }
            }
        }

      // The sensor is outside the model atmosphere, 2D -----------------------
      else
        {
          // Upward observations are not allowed here
          if( abs(rte_los[0]) <= 90 )
            {
              ostringstream os;
              os << "When the sensor is placed outside the model atmosphere,\n"
                 << "upward observations are not allowed.";
              throw runtime_error( os.str() );
            }
          
          // We can here set the path constant, that equals the radius of the
          // geometrical tangent point.
          ppath.constant = geom_tan_pos[0];

          // Handle cases when the sensor appears to look the wrong way
          if( ( rte_pos[1] <= lat_grid[0]  &&  rte_los[0] <= 0 )  || 
                      ( rte_pos[1] >= lat_grid[nlat-1]  &&  rte_los[0] >= 0 ) )
            {
              ostringstream os;
              os << "The sensor is outside (or at the limit) of the model "
                 << "atmosphere but\nlooks in the wrong direction (wrong sign "
                 << "for the zenith angle?).\nThis case includes nadir "
                 << "looking exactly at the latitude end points.";
              throw runtime_error( os.str() );
            }

          // If the sensor is outside the latitude range, check that path is
          // above the closest corner of the model atmosphere
          if( rte_pos[1] < lat_grid[0]  ||  rte_pos[1] > lat_grid[nlat-1] )
            {
              Index   ic = 0;
              String  sc = "lower";
              if( rte_pos[1] > lat_grid[0] )
                { ic = nlat - 1;   sc = "upper"; }
              const double rv = geompath_r_at_lat( ppath.constant, rte_pos[1], 
                                                    rte_los[0], lat_grid[ic] );
              if( rv < ( r_geoid(ic,0) + z_field(np-1,ic,0) ) )
                {
                  ostringstream os;
                  os << "The sensor is outside of the model atmosphere and "
                     << "looks in the\n" << sc << " latitude end face.\n"
                     << "The geometrical altitude of the corner point is "
                     << z_field(np-1,ic,0)/1e3 << " km.\n"
                     << "The geometrical altitude of the entrance point is "
                     << (rv-r_geoid(ic,0))/1e3 << " km.";
                  throw runtime_error( os.str() );
                }
            }

          // If the tangent point is inside covered latitude range, everything
          // is OK. If not, the path must be below the corner of the model
          // atmosphere. 
          if( ( geom_tan_pos[1] < lat_grid[0]  ||  
                                         geom_tan_pos[1] > lat_grid[nlat-1] ) )
            {
              Index   ic = 0;
              String  sc = "lower";
              if( rte_los[0] >= 0 )
                { ic = nlat - 1;   sc = "upper"; }
              const double rv = geompath_r_at_lat( ppath.constant, rte_pos[1], 
                                                    rte_los[0], lat_grid[ic] );
              if( rv >= ( r_geoid(ic,0) + z_field(np-1,ic,0) ) )
                {
                  ostringstream os;
                  os << "The combination of sensor position and line-of-sight "
                     << "gives a\npropagation path that goes above the model "
                     << "atmosphere, with\na tangent point outside the covered"
                     << " latitude range.\nThe latitude of the tangent point "
                     << "is " << geom_tan_pos[1] << " degrees.";
                  throw runtime_error( os.str() );
                }
            }

          // That should be all needed checks. We know now that the path is
          // either totally outside the atmosphere, with a tangent point 
          // inside lat_grid, or enters the atmosphere from the top 
          // somewhere inside lat_grid. In the latter case we need to
          // determine the latitude of the entrance point.
          
          // Path is above the atmosphere:
          // Requieres that tangent point is inside lat_grid and above the
          // top of the atmosphere.
          if( geom_tan_pos[1] >= lat_grid[0]  &&  
                           geom_tan_pos[1] <= lat_grid[nlat-1]   &&  
                                                geom_tan_z >= geom_tan_atmtop )
            {
              ppath_set_background( ppath, 1 );
              out1 << "  ------- WARNING -------: path is totally outside of "
                   << "the model atmosphere\n";
            }

          // The path enters the atmosphere
          else
            {
              // Find the latitude where the path passes top of the atmosphere.

              // We are handling this in a rather dumb way. A test is performed
              // for each latitude range using psurface_crossing_2d.
              // A bit smarter algorithm was considered but that made the code 
              // more messy.
              // The case when the sensor is placed inside lat_grid must be
              // hanled seperetaly.

              // Determine first latitude range of interest, search direction
              // and first test latitude.
              double lat0;
              Index   ilat0, istep;
              //
              if( rte_pos[1] <= lat_grid[0] )
                {
                  lat0  = lat_grid[0]; 
                  ilat0 = 0;
                  istep = 1;
                }
              else if( rte_pos[1] >= lat_grid[nlat-1] )
                {
                  lat0  = lat_grid[nlat-1]; 
                  ilat0 = nlat-1;
                  istep = -1;
                }
              else
                {
                  lat0  = rte_pos[1]; 
                  if( rte_los[0] >= 0 )
                    {  
                      ilat0 = gridpos2gridrange( gp_lat, 1 );
                      istep = 1;
                    }
                  else
                    { 
                      ilat0 = gridpos2gridrange( gp_lat, 0 ) + 1;
                      istep = -1; 
                    }
                }


              // Loop until entrance point is found
              Index ready = 0;
              while( !ready )
                {
                  // The slope c is here calculated for the viewing direction
                  // and the zenith angle to apply shall then be positive.

                  // Calculate radius and zenith angle of path at lat0
                  double r0  = geompath_r_at_lat( ppath.constant, rte_pos[1], 
                                                            rte_los[0], lat0 );
                  double za0 = abs( geompath_za_at_r( ppath.constant, 
                                                            rte_los[0], r0 ) );

                  // Calculate radius and slope to use in psurface_crossing_2d
                  double rv1 = r_geoid(ilat0,0) + z_field(np-1,ilat0,0);
                  double rv2 = r_geoid(ilat0+istep,0) + 
                                                   z_field(np-1,ilat0+istep,0);
                  double latstep = abs( lat_grid[ilat0+istep] - 
                                                             lat_grid[ilat0] );
                  double c = ( rv2 - rv1 ) / latstep;

                  if( lat0 != lat_grid[ilat0] )
                    { rv1 = rv1 + c * abs( lat0 - lat_grid[ilat0] ); } 

                  double dlat = psurface_crossing_2d( r0, za0, rv1, c );

                  if( dlat <= latstep )
                    {
                      ready = 1;
                      ppath.pos(0,1) = lat0 + istep*dlat;
                      ppath.pos(0,0) = rv1 + c * dlat;
                      ppath.los(0,0) = geompath_za_at_r( ppath.constant, 
                                                  rte_los[0], ppath.pos(0,0) );
                      // Re-use some variables from above
                      rv1 = r_geoid(ilat0,0);
                      rv2 = r_geoid(ilat0+istep,0);
                      c   = ( rv2 - rv1 ) / latstep;
                      ppath.z[0] = ppath.pos(0,0) - ( rv1 + c * 
                                     abs( ppath.pos(0,1) - lat_grid[ilat0] ) );
                      ppath.gp_p[0].idx = np - 2;
                      ppath.gp_p[0].fd[0] = 1;
                      ppath.gp_p[0].fd[1] = 0;
                      gridpos( ppath.gp_lat[0], lat_grid, ppath.pos(0,1) );
                    } 
                  else
                    {
                      ilat0 += istep;
                      lat0   = lat_grid[ilat0];
                    }
                } 
            }
        }      

      // Handle possible numerical problems for grid positions
      gridpos_check_fd( ppath.gp_p[0] );
      gridpos_check_fd( ppath.gp_lat[0] );

      // Set geometrical tangent point position
      if( geom_tan_pos.nelem() == 2 )
        {
          ppath.geom_tan_pos.resize(2);
          ppath.geom_tan_pos = geom_tan_pos;
        }

    }  // End 2D


  //-- 3D ---------------------------------------------------------------------
  else
    {
      // Number of points in the latitude and longitude grids
      const Index   nlat = lat_grid.nelem();
      const Index   nlon = lon_grid.nelem();

      // Is the sensor inside the latitude and longitude ranges of the
      // model atmosphere, and below the top of the atmosphere? If
      // yes, is_inside = true. Store geoid and ground radii, grid
      // position and interpolation weights for later use.
      //
      bool      is_inside = false;   
      double   rv_geoid=-1, rv_ground=-1;  // -1 to avoid compiler warnings
      GridPos   gp_lat, gp_lon;
      Vector    itw(4);
      //
      if( rte_pos[1] >= lat_grid[0]  &&  rte_pos[1] <= lat_grid[nlat-1]  &&
                rte_pos[2] >= lon_grid[0]  &&  rte_pos[2] <= lon_grid[nlon-1] )
        {
          gridpos( gp_lat, lat_grid, rte_pos[1] );
          gridpos( gp_lon, lon_grid, rte_pos[2] );
          interpweights( itw, gp_lat, gp_lon );

          rv_geoid  = interp( itw, r_geoid, gp_lat, gp_lon );
          rv_ground = rv_geoid + interp( itw, z_ground, gp_lat, gp_lon );

          out2 << "  sensor altitude        : " << (rte_pos[0]-rv_geoid)/1e3 
               << " km\n";

          if( rte_pos[0] < ( rv_geoid + interp( itw, 
                                z_field(np-1,joker,joker), gp_lat, gp_lon ) ) )
            { is_inside = true; }
        }

      // If downwards, calculate geometrical tangent position. If the tangent
      // point is inside the covered latitude range, calculate also the 
      // geometrical altitude of the tangent point and the top of atmosphere.
      //
      Array<double>  geom_tan_pos(0);
      double geom_tan_z=-2, geom_tan_atmtop=-1;  // OK values if the variables
      //                                             are not set below
      if( rte_los[0] > 90 )
        {
          geom_tan_pos.resize(3);
          double    dummy;
          geompath_tanpos_3d( geom_tan_pos[0], geom_tan_pos[1], 
                   geom_tan_pos[2], dummy, rte_pos[0], rte_pos[1], rte_pos[2], 
                   rte_los[0], rte_los[1], 
                   geometrical_ppc( rte_pos[0], rte_los[0] ) );
          out2 << "  geom. tangent radius   : " << geom_tan_pos[0] / 1e3
               <<" km\n";
          out2 << "  geom. tangent latitude : " << geom_tan_pos[1] 
               << "\n";
          out2 << "  geom. tangent longitude: " << geom_tan_pos[2] 
               << "\n";
          //

          if( geom_tan_pos[1] >= lat_grid[0]       &&  
              geom_tan_pos[1] <= lat_grid[nlat-1]  &&
              geom_tan_pos[2] >= lon_grid[0]       &&
              geom_tan_pos[2] <= lon_grid[nlon-1] )
            {
              GridPos   gp_tanlat, gp_tanlon;
              Vector    itw_tan(4);
              gridpos( gp_tanlat, lat_grid, geom_tan_pos[1] );
              gridpos( gp_tanlon, lon_grid, geom_tan_pos[2] );
              interpweights( itw_tan, gp_tanlat, gp_tanlon );
              geom_tan_z = geom_tan_pos[0] - 
                              interp( itw_tan, r_geoid, gp_tanlat, gp_tanlon );
              geom_tan_atmtop = interp( itw_tan, 
               z_field(np-1,joker,joker), gp_tanlat, gp_tanlon );
              out2 << "  geom. tangent altitude : " << geom_tan_z/1e3 
                   << " km\n";
            }
        }

      // Put sensor position and LOS in ppath as first guess
      ppath.pos(0,0) = rte_pos[0];
      ppath.pos(0,1) = rte_pos[1];
      ppath.pos(0,2) = rte_pos[2];
      ppath.los(0,0) = rte_los[0];
      ppath.los(0,1) = rte_los[1];


      // The sensor is inside the model atmosphere, 3D ------------------------
      if( is_inside )
        {
          // Check that the sensor is above the ground
          if( rte_pos[0] < rv_ground )
            {
              ostringstream os;
              os << "The sensor is placed " 
                 << (rv_ground - rte_pos[0])/1e3 << " km below ground level.\n"
                 << "The sensor must be above the ground.";
              throw runtime_error(os.str());
            }

          // Check that not at latitude end point and looks out
          if( rte_pos[1] > -90  &&  rte_pos[1] == lat_grid[0]  &&  
                                                       abs( rte_los[1] > 90 ) )
            throw runtime_error( "The sensor is at the lower latitude end "
                   "point and the absolute value of the azimuth angle > 90." );
          if( rte_pos[1] < 90  &&  rte_pos[1] == lat_grid[nlat-1]  &&  
                                                     abs( rte_los[1] ) <= 90 ) 
            throw runtime_error( "The sensor is at the upper latitude end "
                  "point and the absolute value of the azimuth angle <= 90." );

          // Check that not at longitude end point and looks out
          if( rte_pos[2] == lon_grid[0]  &&  rte_los[1] < 0 )
            throw runtime_error( "The sensor is at the lower longitude end "
                                          "point and the azimuth angle < 0." );
          if( rte_pos[2] == lon_grid[nlon-1]  &&  rte_los[1] > 0 ) 
            throw runtime_error( "The sensor is at the upper longitude end "
                                          "point and the azimuth angle > 0." );
          
          // Geometrical altitude
          ppath.z[0] = ppath.pos(0,0) - rv_geoid;

          // Use below the values in ppath (instead of rte_pos and rte_los) as 
          // they can be modified on the way.
     
          // Grid positions
          ppath.gp_lat[0].idx   = gp_lat.idx;
          ppath.gp_lat[0].fd[0] = gp_lat.fd[0];
          ppath.gp_lat[0].fd[1] = gp_lat.fd[1];
          ppath.gp_lon[0].idx   = gp_lon.idx;
          ppath.gp_lon[0].fd[0] = gp_lon.fd[0];
          ppath.gp_lon[0].fd[1] = gp_lon.fd[1];

          // Create a vector with the geometrical altitude of the pressure 
          // surfaces for the sensor latitude and use it to get ppath.gp_p.
          Vector z_grid(np);
          z_at_latlon( z_grid, p_grid, lat_grid, lon_grid, z_field, 
                                                              gp_lat, gp_lon );
          gridpos( ppath.gp_p, z_grid, ppath.z );

          // Is the sensor on the ground looking down?
          if( ppath.pos(0,0) == rv_ground )
            {
              // Calculate radial slope of the ground
              const double rslope = psurface_slope_3d( lat_grid, lon_grid,
                          r_geoid, z_ground, gp_lat, gp_lon, ppath.los(0,1)  );

              // Calculate angular tilt of the ground
              const double atilt = psurface_angletilt( rv_ground, rslope);

              // Are we looking down into the ground?
              // If yes and the sensor is inside the cloudbox, the background 
              // will be changed below.
              if( is_los_downwards( ppath.los(0,0), atilt ) )
                { ppath_set_background( ppath, 2 ); }
            }

          // Check sensor position with respect to cloud box.
          if( cloudbox_on )
            {
              // To check all possible cases here when the sensor is at the
              // surface and can either look into or out from the box needs
              // a lot of coding.
              // So we are instead sloppy and set all cases when the sensor
              // is inside or at the surface to be inside the box.
              // The neglected cases should be very unlikely in for real
              // situations.

              if( ppath.pos(0,1) >= lat_grid[cloudbox_limits[2]]  &&
                  ppath.pos(0,1) <= lat_grid[cloudbox_limits[3]]  &&
                  ppath.pos(0,2) >= lon_grid[cloudbox_limits[4]]  &&
                  ppath.pos(0,2) <= lon_grid[cloudbox_limits[5]]  )
                {
                  // Calculate the lower and upper altitude radius limit for
                  // the cloud box at the latitude of the sensor

                  const double   rv_low = rv_geoid + interp( itw, 
                     z_field(cloudbox_limits[0],joker,joker), gp_lat, gp_lon );
                  const double   rv_upp = rv_geoid + interp( itw, 
                     z_field(cloudbox_limits[1],joker,joker), gp_lat, gp_lon );

                  if( ppath.pos(0,0) >= rv_low  &&  ppath.pos(0,0) <= rv_upp )
                    {
                      if( outside_cloudbox )
                        { ppath_set_background( ppath, 4 ); }       
                    }
                  else
                    { assert( outside_cloudbox ); }
                }
              else
                { assert( outside_cloudbox ); }
            }
        }

      // The sensor is outside the model atmosphere, 3D -----------------------
      else
        {
          // Upward observations are not allowed here
          if( abs(rte_los[0]) <= 90 )
            {
              ostringstream os;
              os << "When the sensor is placed outside the model atmosphere,\n"
                 << "upward observations are not allowed.";
              throw runtime_error( os.str() );
            }
          
          // We can here set the path constant, that equals the radius of the
          // geometrical tangent point.
          ppath.constant = geom_tan_pos[0];

          // Handle cases when the sensor appears to look the wrong way in
          // the north-south direction
          if( ( rte_pos[1] <= lat_grid[0]  &&  abs( rte_los[1] ) >= 90 )  || 
              ( rte_pos[1] >= lat_grid[nlat-1]  &&  abs( rte_los[1] ) <= 90 ) )
            {
              ostringstream os;
              os << "The sensor is north or south (or at the limit) of the "
                 << "model atmosphere but\nlooks in the wrong direction.\n";
              throw runtime_error( os.str() );
            }

          // Handle cases when the sensor appears to look the wrong way in
          // the west-east direction. We demand that the sensor is inside the
          // range of lon_grid even if all longitudes are covered.
          if( ( rte_pos[2] <= lon_grid[0]  &&  rte_los[1] < 0 )  || 
                       ( rte_pos[2] >= lon_grid[nlon-1]  &&  rte_los[1] > 0 ) )
            {
              ostringstream os;
              os << "The sensor is east or west (or at the limit) of the "
                 << "model atmosphere but\nlooks in the wrong direction.\n";
              throw runtime_error( os.str() );
            }

          // We either checks that:
          // 1. the tangent point is above the atmosphere, inside covered 
          //    lat and lon ranges
          // 2. We try to determine the entrance point, and issues an error
          //    if it is not at the top of the atmosphere.
          
          // Is tangent point above the model atmosphere
          if( geom_tan_pos[1] >= lat_grid[0]  &&  
                           geom_tan_pos[1] <= lat_grid[nlat-1]   &&  
                           geom_tan_pos[2] >= lon_grid[0]        &&  
                           geom_tan_pos[2] <= lon_grid[nlon-1]   &&  
                                                geom_tan_z >= geom_tan_atmtop )
            {
              ppath_set_background( ppath, 1 );
              out1 << "  ------- WARNING -------: path is totally outside of "
                   << "the model atmosphere\n";
            }

          // The path: does it enter the atmosphere, and then where?
          else
            {
              // Create a matrix for top of the atmosphere radii
              Matrix r_atmtop(nlat,nlon);
              //
              for( Index ilat=0; ilat<nlat; ilat++ )
                {
                  for( Index ilon=0; ilon<nlon; ilon++ )
                    { r_atmtop(ilat,ilon) = r_geoid(ilat,ilon) +
                                                     z_field(np-1,ilat,ilon); }
                }

              // Handle the case when the sensor radius is obviously too low
              double   rtopmin = min( r_atmtop );
              if( rte_pos[0] <= rtopmin )
                {
                  ostringstream os;
                  os << "The sensor radius is smaller than the minimum radius "
                     << "of the top of\nthe atmosphere. This gives no possible"
                     << " OK entrance point for the path.";
                  throw runtime_error( os.str() );
                }

              // Handle the case when the path clearly passes above the 
              // model atmosphere
              double   rtopmax = max( r_atmtop );
              if( geom_tan_pos[0] >= rtopmax )
                {
                  ostringstream os;
                  os << "The combination of sensor position and line-of-sight "
                     << "gives a\npropagation path that goes above the model "
                     << "atmosphere, with\na tangent point outside the covered"
                     << " latitude and longitude ranges.\nThe position of the "
                     << "tangent point is:\n   lat : " << geom_tan_pos[1] 
                     << "\n   lon : " << geom_tan_pos[2];
                  throw runtime_error( os.str() );
                }

              // Sensor pos and LOS in cartesian coordinates
              double   x, y, z, dx, dy, dz;
              poslos2cart( x, y, z, dx, dy, dz, rte_pos[0], 
                              rte_pos[1], rte_pos[2], rte_los[0], rte_los[1] );

              // Boolean for that entrance point CANNOT be found
              bool   failed = false;
             
              // Determine the entrance point for the minimum of *r_atmtop*
              double   r_top, lat_top, lon_top, l_top;
              psurface_crossing_3d( r_top, lat_top, lon_top, l_top, 
                                  rtopmin, rte_pos[0], rte_pos[1], rte_pos[2], 
                                  rte_los[0], x, y, z, dx, dy, dz );

              // If no crossing was found for min radius, or it is outside
              // covered lat and lon ranges, try max radius and make the same
              // check. If check not succesful, there is no OK entrance point.
              if( lat_top < lat_grid[0]  ||  lat_top > lat_grid[nlat-1] ||
                  lon_top < lat_grid[0]  ||  lon_top > lat_grid[nlat-1] )
                {
                  psurface_crossing_3d( r_top, lat_top, lon_top, l_top, 
                                   rtopmax, rte_pos[0], rte_pos[1], rte_pos[2],
                                   rte_los[0], x, y, z, dx, dy, dz );

                  if( lat_top < lat_grid[0]  ||  lat_top > lat_grid[nlat-1] ||
                      lon_top < lon_grid[0]  ||  lon_top > lon_grid[nlon-1] )
                    { failed = true; }
                }

              // Search iteratively for the entrance point until convergence
              // or moving out from covered lat and lon ranges
              //
              bool   ready = false;
              //
              while( !failed  &&  !ready )
                {
                  // Determine radius at found lat and lon
                  GridPos   gp_lattop, gp_lontop;
                  double   lat1, lat3, lon5, lon6, r15, r35, r36, r16;
                  Index     ilat, ilon;
                  gridpos( gp_lattop, lat_grid, lat_top );
                  ilat = gridpos2gridrange( gp_lattop, abs(rte_los[1]) >= 90 );
                  gridpos( gp_lontop, lon_grid, lon_top );
                  ilon  = gridpos2gridrange( gp_lontop, rte_los[1] > 0 );   
                  r15   = r_geoid(ilat,ilon) + z_field(np-1,ilat,ilon);
                  r35   = r_geoid(ilat+1,ilon) + z_field(np-1,ilat+1,ilon);
                  r36   = r_geoid(ilat+1,ilon+1) + z_field(np-1,ilat+1,ilon+1);
                  r16   = r_geoid(ilat,ilon+1) + z_field(np-1,ilat,ilon+1);
                  lat1  = lat_grid[ilat];
                  lat3  = lat_grid[ilat+1];
                  lon5  = lon_grid[ilon];
                  lon6  = lon_grid[ilon+1];
                  r_top = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                       r15, r35, r36, r16, lat_top, lon_top );

                  // Determine new entrance position for latest estimated
                  // entrance radius
                  double   lat_top2, lon_top2;
                  psurface_crossing_3d( r_top, lat_top2, lon_top2, l_top, 
                                    r_top, rte_pos[0], rte_pos[1], rte_pos[2], 
                                    rte_los[0], x, y, z, dx, dy, dz );
                  
                  // Check if new position is sufficiently close to old
                  if( abs( lat_top2 - lat_top ) < 1e-6  &&  
                                             abs( lon_top2 - lon_top ) < 1e-6 )
                    { ready = true; }

                  // Have we moved outside covered lat or lon range?
                  else if( lat_top < lat_grid[0]  ||  
                           lat_top > lat_grid[nlat-1] ||
                           lon_top < lon_grid[0]  ||  
                           lon_top > lon_grid[nlon-1] )
                    { failed = true; }
                  
                  lat_top = lat_top2;   
                  lon_top = lon_top2; 
                }

              if( failed )
                {
                  ostringstream os;
                  os << "The path does not enter the model atmosphere at the "
                     << "top.\nThe path reaches the top of the atmosphere "
                     << "altitude\napproximately at the position:\n"
                     << "   lat : " << lat_top << "\n   lon : " << lon_top;
                  throw runtime_error( os.str() );
                }

              // Correct found lat and lon for some special cases
              //
              if( rte_los[0] == 180 )
                {
                  lat_top = rte_pos[1];
                  lon_top = rte_pos[2];
                }
              else if( abs( rte_pos[1] ) < 90  && ( abs(rte_los[1]) < ANGTOL  
                                       ||  abs( rte_los[1] ) > 180-ANGTOL ) )
                { lon_top = rte_pos[2]; } 

              // Move found values to *ppath*
              //
              // Position
              ppath.pos(0,0) = r_top;
              ppath.pos(0,1) = lat_top;
              ppath.pos(0,2) = lon_top;
              //
              // Grid position
              ppath.gp_p[0].idx = np - 2;
              ppath.gp_p[0].fd[0] = 1;
              ppath.gp_p[0].fd[1] = 0;
              gridpos( ppath.gp_lat[0], lat_grid, lat_top ); 
              gridpos( ppath.gp_lon[0], lon_grid, lon_top ); 
              //
              // Geometrical altitude
              Vector   itw(4);
              interpweights( itw, ppath.gp_lat[0], ppath.gp_lon[0] );
              ppath.z[0] = ppath.pos(0,0) - interp(itw,  r_geoid,
                                            ppath.gp_lat[0], ppath.gp_lon[0] );
              //
              // LOS
              double   zatmp, aatmp;
              cart2poslos( r_top, lat_top, lon_top, zatmp, aatmp,
                           x+dx*l_top, y+dy*l_top, z+dz*l_top, dx, dy, dz );
              ppath.los(0,0) = zatmp;
              ppath.los(0,1) = aatmp;
              //
              // Correct found LOS for some special cases
              if( rte_los[0] == 180 )
                { ppath.los(0,0) = 180; }
              if( abs( rte_pos[1] ) < 90  &&  ( abs( rte_los[1] ) < ANGTOL  
                                         ||  abs( rte_los[1] ) > 180-ANGTOL ) )
                { ppath.los(0,1) = rte_los[1]; } 
            }

        } // else

      // Handle possible numerical problems for grid positions
      gridpos_check_fd( ppath.gp_p[0] );
      gridpos_check_fd( ppath.gp_lat[0] );
      gridpos_check_fd( ppath.gp_lon[0] );

      // Set geometrical tangent point position
      if( geom_tan_pos.nelem() == 3 )
        {
          ppath.geom_tan_pos.resize(3);
          ppath.geom_tan_pos[0] = geom_tan_pos[0];
          ppath.geom_tan_pos[1] = geom_tan_pos[1];          
          ppath.geom_tan_pos[2] = geom_tan_pos[2];        
        }
    }  // End 3D
}





//! ppath_calc
/*! 
   This is the core for the WSM ppathCalc.

   This function takes the same input as ppathCalc (that is, those
   input arguments are the WSV with the same name), but there are some
   additional argument(s):

   \param ppath              Output: A Ppath structure
   \param ppath_step         FIXME: Add documentation.
   \param ppath_step_agenda  FIXME: Add documentation.
   \param atmosphere_dim     The atmospheric dimensionality.
   \param p_grid             The pressure grid.
   \param lat_grid           The latitude grid.
   \param lon_grid           The longitude grid.
   \param z_field            The field of geometrical altitudes.
   \param r_geoid            The geoid radius.
   \param z_ground           Ground altitude.
   \param cloudbox_on        Flag to activate the cloud box.
   \param cloudbox_limits    Index limits of the cloud box.
   \param rte_pos            The position of the sensor.
   \param rte_los            The line-of-sight of the sensor.
   \param outside_cloudbox   Boolean to flag if the propagation path is 
                             (expected to be) outside the cloudbox. Ordinary
                             clerar sky calculations are selected by the value
                             1. The value 0 means tracking of a propagation 
                             path inside the cloudbox. The path is then
                             tracked to the cloudbox boundary.
   \param agenda_verb        This argument is given as input to agendas
                             to control the verbosity.

   \author Patrick Eriksson
   \date   2003-01-08
*/
void ppath_calc(
        // WS Output:
              Ppath&          ppath,
              Ppath&          ppath_step,
        // WS Input:
        const Agenda&         ppath_step_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         rte_pos,
        const Vector&         rte_los,
        const bool&           outside_cloudbox,
        const Index&          agenda_verb )
{
  // This function is a WSM but it is normally only called from RteCalc. 
  // For that reason, this function does not repeat input checks that are
  // performed in RteCalc, it only performs checks regarding the sensor 
  // position and LOS.

  //--- Check input -----------------------------------------------------------

  // Sensor position and LOS
  //
  chk_vector_length( "rte_pos", rte_pos, atmosphere_dim );
  chk_if_over_0( "sensor radius", rte_pos[0] );
  if( atmosphere_dim == 1 )
    {
      chk_vector_length( "rte_los", rte_los, 1 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], 0, 180 );
    }
  else if( atmosphere_dim == 2 )
    {
      chk_vector_length( "rte_los", rte_los, 1 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], -180, 180 );
      if( cloudbox_on )
        { throw runtime_error( "The cloud box is not defined for 2D." ); }
    }
  else
    {
      chk_if_in_range( "sensor latitude", rte_pos[1], -90, 90 );
      chk_if_in_range( "sensor longitude", rte_pos[2], -360, 360 );
      chk_vector_length( "rte_los", rte_los, 2 );
      chk_if_in_range( "sensor zenith angle", rte_los[0], 0, 180 );
      chk_if_in_range( "sensor azimuth angle", rte_los[1], -180, 180 );
    }
  
  //--- End: Check input ------------------------------------------------------


  // Some messages
  out2 << "  -------------------------------------\n";
  out2 << "  sensor radius          : " << rte_pos[0]/1e3 << " km\n";
  if( atmosphere_dim >= 2 )    
    out2 << "  sensor latitude        : " << rte_pos[1] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor longitude       : " << rte_pos[2] << "\n";
  out2 << "  sensor zenith angle    : " << rte_los[0] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor azimuth angle   : " << rte_los[1] << "\n";


  // Initiate the partial Ppath structure. 
  // The function doing the work sets ppath_step to the point of the path
  // inside the atmosphere closest to the sensor, if the path is at all inside
  // the atmosphere.
  // If the background field is set by the function this flags that there is no
  // path to follow (for example when the sensor is inside the cloud box).
  // The function checks also that the sensor position and LOS give an
  // allowed path.
  //
  ppath_start_stepping( ppath_step, atmosphere_dim, p_grid, lat_grid, 
            lon_grid, z_field, r_geoid, z_ground,
            cloudbox_on, cloudbox_limits, outside_cloudbox, rte_pos, rte_los );

  out2 << "  -------------------------------------\n";

  // Perform propagation path steps until the starting point is found, which
  // is flagged by ppath_step by setting the background field.
  //
  // The results of each step, returned by ppath_step_agenda as a new 
  // ppath_step, are stored as an array of Ppath structures.
  //
  Array<Ppath>   ppath_array;
  ppath_array.push_back( ppath_step );
  // 
  Index   np = ppath_step.np;   // Counter for number of points of the path
  Index   istep = 0;            // Counter for number of steps
  //
  const Index imax_p   = p_grid.nelem() - 1;
  const Index imax_lat = lat_grid.nelem() - 1;
  const Index imax_lon = lon_grid.nelem() - 1;

  
  while( !ppath_what_background( ppath_step ) )
    {

      // Call ppath_step agenda. 
      // The new path step is added to *ppath_array* last in the while block
      //
      istep++;
      //
      ppath_step_agenda.execute( agenda_verb + ( istep - 1 ) );

      // Before everything is tested carefully, we consider more than 5000
      // path points to be an indication on that the calcululations have
      // got stuck in an infinite loop.
      if( istep > 5000 )
        {
          throw logic_error(
             "5000 path points have been reached. Is this an infinite loop?" );
        }
      
      // Number of points in returned path step
      const Index n = ppath_step.np;

      // Increase the total number
      np += n - 1;

      // Check if the top of the atmosphere is reached
      if( is_gridpos_at_index_i( ppath_step.gp_p[n-1], imax_p ) )
        { ppath_set_background( ppath_step, 1 ); }

      // Check that path does not exit at a latitude or longitude end face
      if( atmosphere_dim == 2 )
        {
          // Latitude 
          if( is_gridpos_at_index_i( ppath_step.gp_lat[n-1], 0 ) )
            {
              ostringstream os;
              os << "The path exits the atmosphere through the lower latitude"
                 << " end face.\nThe exit point is at an altitude of " 
                 << ppath_step.z[n-1]/1e3 << " km.";
              throw runtime_error( os.str() );
            }
          if( is_gridpos_at_index_i( ppath_step.gp_lat[n-1], imax_lat ) )
            {
              ostringstream os;
              os << "The path exits the atmosphere through the upper latitude"
                 << " end face.\nThe exit point is at an altitude of " 
                 << ppath_step.z[n-1]/1e3 << " km.";
              throw runtime_error( os.str() );
            }
        }
      if( atmosphere_dim == 3 )
        {
          // Latitude 
          if( lat_grid[0] > -90  && 
                           is_gridpos_at_index_i( ppath_step.gp_lat[n-1], 0 ) )
            {
              ostringstream os;
              os << "The path exits the atmosphere through the lower latitude"
                 << " end face.\nThe exit point is at an altitude of " 
                 << ppath_step.z[n-1]/1e3 << " km.";
              throw runtime_error( os.str() );
            }
          if( lat_grid[imax_lat] < 90  && 
                    is_gridpos_at_index_i( ppath_step.gp_lat[n-1], imax_lat ) )
            {
              ostringstream os;
              os << "The path exits the atmosphere through the upper latitude"
                 << " end face.\nThe exit point is at an altitude of " 
                 << ppath_step.z[n-1]/1e3 << " km.";
              throw runtime_error( os.str() );
            }

          // Longitude 
          // Note that it must be if and else if here. Otherwise e.g. -180 
          // will be shifted to 180 and then later back to -180.
          if( is_gridpos_at_index_i( ppath_step.gp_lon[n-1], 0 )  &&
             ppath_step.los(n-1,1) < 0  &&  abs( ppath_step.pos(n-1,1) ) < 90 )
            {
              // Check if the longitude point can be shifted +360 degrees
              if( lon_grid[imax_lon] - lon_grid[0] >= 360 )
                {
                  ppath_step.pos(n-1,2) = ppath_step.pos(n-1,2) + 360;
                  gridpos( ppath_step.gp_lon[n-1], lon_grid, 
                                                       ppath_step.pos(n-1,2) );
                }
              else
                {
                  ostringstream os;
                  os << "The path exits the atmosphere through the lower " 
                     << "longitude end face.\nThe exit point is at an "
                     << "altitude of " << ppath_step.z[n-1]/1e3 << " km.";
                  throw runtime_error( os.str() );
                }
            }
          else if( is_gridpos_at_index_i( ppath_step.gp_lon[n-1], imax_lon ) &&
             ppath_step.los(n-1,1) > 0  &&  abs( ppath_step.pos(n-1,1) ) < 90 )
            {
              // Check if the longitude point can be shifted -360 degrees
              if( lon_grid[imax_lon] - lon_grid[0] >= 360 )
                {
                  ppath_step.pos(n-1,2) = ppath_step.pos(n-1,2) - 360;
                  gridpos( ppath_step.gp_lon[n-1], lon_grid, 
                                                       ppath_step.pos(n-1,2) );
                }
              else
                {
                  ostringstream os;
                  os << "The path exits the atmosphere through the upper "
                     << "longitude end face.\nThe exit point is at an "
                     << "altitude of " << ppath_step.z[n-1]/1e3 << " km.";
                  throw runtime_error( os.str() );
                }
            }
        }
      
    
      // Check if there is an intersection with an active cloud box
      if( cloudbox_on )
        {
          double ipos = double( ppath_step.gp_p[n-1].idx ) + 
                                                    ppath_step.gp_p[n-1].fd[0];
          if( ipos >= double( cloudbox_limits[0] )  && 
                                        ipos <= double( cloudbox_limits[1] ) )
            {
              if( atmosphere_dim == 1 )
                { ppath_set_background( ppath_step, 3 ); }
              else
                {
                  ipos = double( ppath_step.gp_lat[n-1].idx ) + 
                                                  ppath_step.gp_lat[n-1].fd[0];
                  if( ipos >= double( cloudbox_limits[2] )  && 
                                        ipos <= double( cloudbox_limits[3] ) )
                    {
                      ipos = double( ppath_step.gp_lon[n-1].idx ) + 
                                                  ppath_step.gp_lon[n-1].fd[0];
                      if( ipos >= double( cloudbox_limits[4] )  && 
                                        ipos <= double( cloudbox_limits[5] ) )
                        { ppath_set_background( ppath_step, 3 ); } 
                    }
                }
            }
        }

      // Put new ppath_step in ppath_array
      ppath_array.push_back( ppath_step );

    } // End path steps
  
 
  // Combine all structures in ppath_array to form the return Ppath structure.
  //
  ppath_init_structure( ppath, atmosphere_dim, np );
  //
  np = 0;   // Now used as counter for points moved to ppath
  //
  for( Index i=0; i<ppath_array.nelem(); i++ )
    {
      // For the first structure, the first point shall be included, but the
      // first structure can also be empty. 
      // For later structures, the first point shall not be included, but
      // there will always be at least two points.
      // Only the first structure can be empty.

      Index n = ppath_array[i].np;

      if( n )
        {
          // First index to include
          Index i1 = 1;
          if( i == 0 )
            { i1 = 0; }
          else
            { assert( n > 1 ); }

          // Vectors and matrices that can be handled by ranges.
          ppath.z[ Range(np,n-i1) ] = ppath_array[i].z[ Range(i1,n-i1) ];
          ppath.pos( Range(np,n-i1), joker ) = 
                                   ppath_array[i].pos( Range(i1,n-i1), joker );
          ppath.los( Range(np,n-i1), joker ) = 
                                   ppath_array[i].los( Range(i1,n-i1), joker );

          // For i==1, there is no defined l_step. For higher i, all 
          // values in l_step shall be copied.
          if( i > 0 )
            { ppath.l_step[ Range(np-1,n-1) ] = ppath_array[i].l_step; }

          // Grid positions must be handled by a loop
          for( Index j=i1; j<n; j++ )
            { ppath.gp_p[np+j-i1] = ppath_array[i].gp_p[j]; }
          if( atmosphere_dim >= 2 )
            {
              for( Index j=i1; j<n; j++ )
                { ppath.gp_lat[np+j-i1] = ppath_array[i].gp_lat[j]; }
            }
          if( atmosphere_dim == 3 )
            {
              for( Index j=i1; j<n; j++ )
                { ppath.gp_lon[np+j-i1] = ppath_array[i].gp_lon[j]; }
            }

          // Fields just set once
          if( ppath_array[i].tan_pos.nelem() )
            {
              ppath.tan_pos.resize( ppath_array[i].tan_pos.nelem() );
              ppath.tan_pos               = ppath_array[i].tan_pos; 
            }
          if( ppath_array[i].geom_tan_pos.nelem() )
            {
              ppath.geom_tan_pos.resize( ppath_array[i].tan_pos.nelem() );
              ppath.geom_tan_pos          = ppath_array[i].geom_tan_pos; 
            }

          // Increase number of points done
          np += n - i1;
         
        }
    }  
  ppath.method     = ppath_step.method;
  ppath.refraction = ppath_step.refraction;
  ppath.constant   = ppath_step.constant;
  ppath.background = ppath_step.background;

  out3 << "  number of path steps  : " << istep           << "\n";
  out3 << "  number of path points : " << ppath.z.nelem() << "\n";


  // If refraction has been considered, make a simple check that the
  // refraction at the top of the atmosphere is sufficiently close to 1.
  if( ppath.refraction  &&  min( z_field(z_field.npages()-1,0,0) ) < 60e3 )
    {
      out2 << "  *** WARNING****\n" 
           << "  The calculated propagation path can be inexact as the "
           << "atmosphere\n  only extends to " 
           <<  min( z_field(z_field.npages()-1,0,0) ) << " km. \n" 
           << "  The importance of this depends on the observation "
           << "geometry.\n  It is recommended that the top of the atmosphere "
           << "is not below 60 km.\n";
    }
}
