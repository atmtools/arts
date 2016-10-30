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
   \file   geodetic.cc
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date   2012-02-06 

   This file contains functions associated with the reference ellipsoid,
   conversion between latitudes and similar stuff.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "geodetic.h"
#include "math_funcs.h"
#include "ppath.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



/*===========================================================================
  === 2D functions
  ===========================================================================*/

// The 2D case is treated as being the 3D x/z-plane. That is, y-coordinate is
// skipped. For simplicity, the angle coordinate is denoted as latitude.
// However, the latitude is here not limited to [-90,90]. It is cyclic and can
// have any value. The input *lat0* is used to shift the output from atan2 with
// n*360 to return what should be the expected latitude. That is, it is assumed
// that no operation moves the latitude more than 180 degrees from the initial
// value *lat0*. Negative zeniath angles are handled, following ARTS definition
// of 2D geometry.


//! cart2pol
/*! 
   The inverse of *pol2cart*. 
   
   A 2D version of cart2sph 

   \param   r     Out: Radius of position.
   \param   lat   Out: Latitude of position.
   \param   x     x-coordinate of position.
   \param   z     z-coordinate of position.
   \param   lat0  Original latitude.
   \param   za0   Original zenith angle.

   \author Patrick Eriksson
   \date   2012-03-20
*/
void cart2pol(
            Numeric&   r,
            Numeric&   lat,
      const Numeric&   x,
      const Numeric&   z,
      const Numeric&   lat0,
      const Numeric&   za0 )
{
  r   = sqrt( x*x + z*z );

  // Zenith and nadir cases
  const Numeric absza = abs( za0 );
  if( absza < ANGTOL  ||  absza > 180-ANGTOL  )
    { lat = lat0; }

  else
    { // Latitude inside [0,360]
      lat = RAD2DEG * atan2( z, x );
      // Shift with n*360 to get as close as possible to lat0
      lat = lat - 360.0 * Numeric( round( ( lat - lat0 ) / 360.0 ) );
    }
}



//! cart2poslos
/*! 
   2D version of the 3D *cart2poslos*.

   \param   r     Out: Radius of observation position.
   \param   lat   Out: Latitude of observation position.
   \param   za    Out: LOS zenith angle at observation position.
   \param   x     x-coordinate of observation position.
   \param   z     z-coordinate of observation position.
   \param   dx    x-part of LOS unit vector.
   \param   dz    z-part of LOS unit vector.
   \param   ppc   Propagation path constant (r*sin(za))
   \param   lat0  Original latitude.
   \param   za0   Original zenith angle.

   \author Patrick Eriksson
   \date   2012-03-21
*/
void cart2poslos(
             Numeric&   r,
             Numeric&   lat,
             Numeric&   za,
       const Numeric&   x,
       const Numeric&   z,
       const Numeric&   dx,
       const Numeric&   dz,
       const Numeric&   ppc,
       const Numeric&   lat0,
       const Numeric&   za0 )
{
  r   = sqrt( x*x + z*z );

  // Zenith and nadir cases
  const Numeric absza = abs( za0 );
  if( absza < ANGTOL  ||  absza > 180-ANGTOL  )
    { 
      lat = lat0;
      za  = za0; 
    }

  else
    {
      lat = RAD2DEG * atan2( z, x );

      const Numeric   latrad = DEG2RAD * lat;
      const Numeric   coslat = cos( latrad );
      const Numeric   sinlat = sin( latrad );
      const Numeric   dr     = coslat*dx + sinlat*dz;

      // Use ppc for max accuracy, but dr required to resolve if up- 
      // and downward cases.

      // Another possibility to obtain (absolute value of) za is 
      // RAD2DEG*acos(dr).
      // It is checked that the two ways give consistent results, but
      // occasionally deviate with 1e-4 deg (due to numerical issues).

      za = RAD2DEG * asin( ppc / r );
      if( za0 > 0 )
        {
          if( isnan( za ) )
            { za = 90; }
          else if( dr < 0 )
            { za = 180.0 - za; }
        }
      else
        {
          if( isnan( za ) )
            { za = -90; }
          else if( dr < 0 )
            { za = -180.0 + za; }
          else 
            { za = -za; }
        }
    }
}



//! distance2D
/*! 
   The distance between two 2D points.
   
   The two latitudes can deviate with max 180 degrees.

   \param   l     Out: The distance
   \param   r1    Radius of position 1
   \param   lat1  Latitude of position 1
   \param   r2    Radius of position 2
   \param   lat2  Latitude of position 2

   \author Patrick Eriksson
   \date   2012-03-20
*/
void distance2D(
            Numeric&   l,
      const Numeric&   r1,
      const Numeric&   lat1,
      const Numeric&   r2,
      const Numeric&   lat2 )
{
  assert( abs( lat2 - lat1 ) <= 180 );

  Numeric x1, z1, x2, z2;
  pol2cart( x1, z1, r1, lat1 );
  pol2cart( x2, z2, r2, lat2 );

  const Numeric dx = x2 - x1; 
  const Numeric dz = z2 - z1; 
  l = sqrt( dx*dx + dz*dz ); 
}



//! geomtanpoint2d
/*! 
   Position of the tangent point for 3D cases.

   Calculates the 3D geometrical tangent point for arbitrary reference
   ellipsiod. That is, a spherical planet is not assumed. The tangent
   point is thus defined as the point closest to the ellipsoid (not as the
   ppoint with za=90).
  
   Geocentric coordinates are used for both sensor and tangent point
   positions.

   The algorithm used for non-spherical cases is derived by Nick Lloyd at
   University of Saskatchewan, Canada (nick.lloyd@usask.ca), and is part of
   the operational code for both OSIRIS and SMR on-board- the Odin
   satellite.

   The zenith angle must be >= 90.

   \param   r_tan       Out: Radius of tangent point.
   \param   lat_tan     Out: Latitude of tangent point.
   \param   lon_tan     Out: Longitude of tangent point.
   \param   r           Radius of observation position.
   \param   lat         Latitude of observation position.
   \param   lon         Longitude of observation position.
   \param   za          LOS zenith angle at observation position.
   \param   aa          LOS azimuth angle at observation position.

   \author Patrick Eriksson
   \date   2012-02-12
*/
/*
void geomtanpoint2d( 
             Numeric&    r_tan,
             Numeric&    lat_tan,
     ConstVectorView    refellipsoid,
       const Numeric&    r,
       const Numeric&    lat,
       const Numeric&    za )
{
  assert( refellipsoid.nelem() == 2 );
  assert( refellipsoid[0] > 0 );
  assert( refellipsoid[1] >= 0 );
  assert( refellipsoid[1] < 1 );
  assert( r > 0 );
  assert( za >= 90 );
                                // e=1e-7 corresponds to that polar radius
  if( refellipsoid[1] < 1e-7 )  // less than 1 um smaller than equatorial 
    {                           // one for the Earth
      r_tan = geometrical_ppc( r, za );
      if( za > 0 )
        { lat_tan = geompath_lat_at_za( za, lat, 90 ); }
      else
        { lat_tan = geompath_lat_at_za( za, lat, -90 ); }
    }

  else
    {
      assert( 0 );  // To be implemented
    }  
}  
*/


//! line_circle_intersect
/*! 
   Find the intersection between a line and a circle

   \param   x     Out: X-coordinate of intersection
   \param   z     Out: Z-coordinate of intersection
   \param   xl    A x-coordinate on the line
   \param   zl    A z-coordinate on the line
   \param   dx    X-component of line direction vector
   \param   dz    Z-component of line direction vector
   \param   xc    X-coordinate of sphere origo
   \param   zc    Z-coordinate of sphere origo
   \param   r     Radius of sphere

   \author Patrick Eriksson
   \date   2012-03-30
*/
void line_circle_intersect(
         Numeric&   x,
         Numeric&   z,
   const Numeric&   xl,
   const Numeric&   zl,
   const Numeric&   dx,
   const Numeric&   dz,
   const Numeric&   xc,
   const Numeric&   zc,
   const Numeric&   r )
{
  const Numeric a = dx*dx + dz*dz;
  const Numeric b = 2*( dx*(xl-xc) + dz*(zl-zc) );
  const Numeric c = xc*xc + zc*zc + xl*xl + zl*zl - 2*(xc*xl + zc*zl) - r*r;
  
  Numeric d = b*b - 4*a*c;
  assert( d > 0 );
  
  const Numeric a2 = 2*a;
  const Numeric b2 = -b / a2;
  const Numeric e  = sqrt( d ) / a2;

  const Numeric l1 = b2 + e;
  const Numeric l2 = b2 - e;

  Numeric l;
  if( l1 < 0 )
    { l = l2; }
  else  if( l2 < 0 )
    { l = l1; }
  else
    { l = min(l1,l2); assert( l>=0 ); }

  x = xl + l*dx;
  z = zl + l*dz;
}



//! pol2cart
/*! 
   Conversion from polar to cartesian coordinates.

   The cartesian coordinate system is defined such as the x-axis goes along
   lat=0 and lon=0 and z-axis goes along lat=90.

   \param   x     Out: x position.
   \param   z     Out: z position.
   \param   r     Radius.
   \param   lat   Latitude.

   \author Patrick Eriksson
   \date   2012-03-20
*/
void pol2cart(
            Numeric&   x,
            Numeric&   z,
      const Numeric&   r,
      const Numeric&   lat )
{
  assert( r > 0 );

  const Numeric   latrad = DEG2RAD * lat;

  x = r * cos( latrad );  
  z = r * sin( latrad );
}



//! poslos2cart
/*! 
   2D version of poslos2cart

   \param   x     Out: x-coordinate of observation position.
   \param   z     Out: z-coordinate of observation position.
   \param   dx    Out: x-part of LOS unit vector.
   \param   dz    Out: z-part of LOS unit vector.
   \param   r     Radius of observation position.
   \param   lat   Latitude of observation position.
   \param   za    LOS zenith angle at observation position.

   \author Patrick Eriksson
   \date   2012-03-20
*/
void poslos2cart(
             Numeric&   x,
             Numeric&   z,
             Numeric&   dx,
             Numeric&   dz,
       const Numeric&   r,
       const Numeric&   lat,
       const Numeric&   za )
{
  assert( r > 0 );
  assert( za >= -180 && za<=180 );

  const Numeric   latrad = DEG2RAD * lat;
  const Numeric   zarad  = DEG2RAD * za;

  const Numeric   coslat = cos( latrad );
  const Numeric   sinlat = sin( latrad );
  const Numeric   cosza  = cos( zarad );
  const Numeric   sinza  = sin( zarad );

  // This part as pol2cart but uses local variables
  x = r * coslat;  
  z = r * sinlat;

  const Numeric   dr   = cosza;
  const Numeric   dlat = sinza;         // r-term cancel out below

  dx = coslat * dr - sinlat * dlat;
  dz = sinlat * dr + coslat * dlat;
}





/*===========================================================================
  === 3D functions
  ===========================================================================*/

//! cart2poslos
/*! 
   The inverse of *poslos2cart*.

   Beside the cartesian coordinates (x,y,z,dx,dy,dz), the function takes as
   input information about the original position and LOS. The later data are
   used to improve the accuracy. For example, for zenith and nadir cases it is
   ensured that the latitude and longitude are not changed. This makes the
   function slower, but accuarcy is favoured as zenith, nadir, north and south
   line-of-sights are especially tricky as they can go along a grid box
   boundary and the smallest rounding error can move the path from one grid box
   to the neighbouring one.

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
   \param   ppc   Propagation path constant (r*sin(za))
   \param   lat0  Original latitude.
   \param   lon0  Original longitude.
   \param   za0   Original zenith angle.
   \param   aa0   Original azimuth angle.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void cart2poslos(
             Numeric&   r,
             Numeric&   lat,
             Numeric&   lon,
             Numeric&   za,
             Numeric&   aa,
       const Numeric&   x,
       const Numeric&   y,
       const Numeric&   z,
       const Numeric&   dx,
       const Numeric&   dy,
       const Numeric&   dz,
       const Numeric&   ppc,
       const Numeric&   x0,
       const Numeric&   y0,
       const Numeric&   z0,
       const Numeric&   lat0,
       const Numeric&   lon0,
       const Numeric&   za0,
       const Numeric&   aa0 )
{
  // Radius of new point
  r = sqrt( x*x + y*y + z*z );

  // Zenith and nadir cases
  if( za0 < ANGTOL  ||  za0 > 180-ANGTOL  )
    { 
      lat = lat0;
      lon = lon0;
      za  = za0; 
      aa  = aa0; 
    }

  else
    {
      lat = RAD2DEG * asin( z / r );
      lon = RAD2DEG * atan2( y, x );

      bool ns_case = false;
      bool lon_flip = false;

      // Make sure that lon is maintained for N-S cases (if not 
      // starting on a pole)
      if( ( abs(aa0) < ANGTOL  ||  abs(180-aa0) < ANGTOL )  && 
                                             abs( lat0 ) <= POLELAT )
        {
          ns_case = true;
          // Check that not lon changed with 180 deg
          if( abs(lon-lon0) < 1 )
            { lon = lon0; }
          else
            {
              lon_flip = true;
              if( lon0 > 0 )
                { lon = lon0 - 180; }
              else
                { lon = lon0 + 180; }
            }
        }

      const Numeric   latrad = DEG2RAD * lat;
      const Numeric   lonrad = DEG2RAD * lon;
      const Numeric   coslat = cos( latrad );
      const Numeric   sinlat = sin( latrad );
      const Numeric   coslon = cos( lonrad );
      const Numeric   sinlon = sin( lonrad );
      // dr only used to validate za
      const Numeric   dr     = coslat*coslon*dx + coslat*sinlon*dy + sinlat*dz;

      // Set za by ppc for max accuracy, but this does not resolve
      // za and 180-za. This was first resolved by dr, but using l and lmax was
      // foudn to be more stable.
      za = RAD2DEG * asin( ppc / r );

      // Correct and check za
      if( isnan( za ) )
        { za = 90; }
      // If za0 > 90, then correct za could be 180-za. Resolved by checking if
      // the tangent point is passed or not
      if( za0 > 90 )
        {
          const Numeric l = sqrt( pow(x-x0,2.0) + pow(y-y0,2.0) + pow(z-z0,2.0) );
          const Numeric r0 = sqrt( x0*x0 + y0*y0 + z0*z0 );
          const Numeric ltan = geompath_l_at_r( ppc, r0 );
          if( l < ltan )
            { za = 180.0 - za; }
        }
      // The difference below can at least be 1e-5 for tangent points
      if( abs( za - RAD2DEG*acos(dr) ) >= 1e-4 )
        {
          throw runtime_error(
            "Internal consistency check in *cart2poslos failed. If this "
            "happens to you and you critically need the ongoing calculations, "
            "try to change the observation LOS slightly. If you can reproduce "
            "this error, please contact Patrick in order to help tracking down "
            "the reason to this problem. If you see this message occasionally "
            "when doing MC calculations, it should not be critical. This path "
            "sampling will be rejected and replaced with a new one." );
        }

      // As last check of za, make sure that it has decreased
      //if( za >= za0 )
      //  {
      //    za = za0 - 1e-4;
      //    cout << "Increasing za! Reset to = " << za << " / ";          
      //  }

      // For lat = +- 90 the azimuth angle gives the longitude along which 
      // the LOS goes
      if( abs( lat ) >= POLELAT )      
        { aa = RAD2DEG * atan2( dy, dx ); }

      // N-S cases, not starting at a pole
      else if( ns_case )
        { 
          if( !lon_flip )
            { aa = aa0; }
          else
            {
              if( abs(aa0) < ANGTOL )
                { aa = 180; }
              else
                { aa = 0; }
            }
        }

      else
        {
          const Numeric   dlat = -sinlat*coslon/r*dx - sinlat*sinlon/r*dy + 
                                                             coslat/r*dz;
          const Numeric   dlon = -sinlon/coslat/r*dx + coslon/coslat/r*dy;

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
    }
}



//! cart2sph
/*! 
   The inverse of *sph2cart*.

   For definition of lat0, lon0, za0 and aa0, see *cart2poslos*.

   \param   r     Out: Radius of observation position.
   \param   lat   Out: Latitude of observation position.
   \param   lon   Out: Longitude of observation position.
   \param   x     x-coordinate of observation position.
   \param   y     y-coordinate of observation position.
   \param   z     z-coordinate of observation position.
   \param   lat0  Original latitude.
   \param   lon0  Original longitude.
   \param   za0   Original zenith angle.
   \param   aa0   Original azimuth angle.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void cart2sph(
             Numeric&   r,
             Numeric&   lat,
             Numeric&   lon,
       const Numeric&   x,
       const Numeric&   y,
       const Numeric&   z,
       const Numeric&   lat0,
       const Numeric&   lon0,
       const Numeric&   za0,
       const Numeric&   aa0 )
{
  r   = sqrt( x*x + y*y + z*z );

  // Zenith and nadir cases
  if( za0 < ANGTOL  ||  za0 > 180-ANGTOL  )
    { 
      lat = lat0;
      lon = lon0;
    }

  else
    {
      lat = RAD2DEG * asin( z / r );
      lon = RAD2DEG * atan2( y, x );

      // Make sure that lon is maintained for N-S cases (if not 
      // starting on a pole)
      if( ( abs(aa0) < ANGTOL  ||  abs(180-aa0) < ANGTOL )  && 
                                             abs( lat0 ) <= POLELAT )
        {
          // Check that not lon changed with 180 deg
          if( abs(lon-lon0) < 1 )
            { lon = lon0; }
          else
            {
              if( lon0 > 0 )
                { lon = lon0 - 180; }
              else
                { lon = lon0 + 180; }
            }
        }
    }
}



//! distance3D
/*! 
   The distance between two 3D points.
   
   \param   l     Out: The distance
   \param   r1    Radius of position 1
   \param   lat1  Latitude of position 1
   \param   lon1  Longitude of position 1
   \param   r2    Radius of position 2
   \param   lat2  Latitude of position 2
   \param   lon2  Longitude of position 2

   \author Patrick Eriksson
   \date   2012-03-20
*/
void distance3D(
            Numeric&   l,
      const Numeric&   r1,
      const Numeric&   lat1,
      const Numeric&   lon1,
      const Numeric&   r2,
      const Numeric&   lat2,
      const Numeric&   lon2 )
{
  Numeric x1, y1, z1, x2, y2, z2;
  sph2cart( x1, y1, z1, r1, lat1, lon1 );
  sph2cart( x2, y2, z2, r2, lat2, lon2 );

  const Numeric dx = x2 - x1; 
  const Numeric dy = y2 - y1; 
  const Numeric dz = z2 - z1; 
  l = sqrt( dx*dx + dy*dy + dz*dz ); 
}



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
   \param   ppc         Geometrical propagation path constant.

   \author Patrick Eriksson
   \date   2002-12-31
*/
void geompath_tanpos_3d( 
             Numeric&    r_tan,
             Numeric&    lat_tan,
             Numeric&    lon_tan,
             Numeric&    l_tan,
       const Numeric&    r,
       const Numeric&    lat,
       const Numeric&    lon,
       const Numeric&    za,
       const Numeric&    aa,
       const Numeric&    ppc )
{
  assert( za >= 90 );
  assert( r >= ppc );

  Numeric   x, y, z, dx, dy, dz; 

  poslos2cart( x, y, z, dx, dy, dz, r, lat, lon, za, aa );

  l_tan = sqrt( r*r - ppc*ppc );

  cart2sph( r_tan, lat_tan, lon_tan, x+dx*l_tan, y+dy*l_tan, z+dz*l_tan,
            lat, lon, za, aa );
}



//! geomtanpoint
/*! 
   Position of the tangent point for 3D cases.

   Calculates the 3D geometrical tangent point for arbitrary reference
   ellipsiod. That is, a spherical planet is not assumed. The tangent
   point is thus defined as the point closest to the ellipsoid (not as the
   ppoint with za=90).
  
   Geocentric coordinates are used for both sensor and tangent point
   positions.

   The algorithm used for non-spherical cases is derived by Nick Lloyd at
   University of Saskatchewan, Canada (nick.lloyd@usask.ca), and is part of
   the operational code for both OSIRIS and SMR on-board- the Odin
   satellite.

   The zenith angle must be >= 90.

   \param   r_tan       Out: Radius of tangent point.
   \param   lat_tan     Out: Latitude of tangent point.
   \param   lon_tan     Out: Longitude of tangent point.
   \param   r           Radius of observation position.
   \param   lat         Latitude of observation position.
   \param   lon         Longitude of observation position.
   \param   za          LOS zenith angle at observation position.
   \param   aa          LOS azimuth angle at observation position.

   \author Patrick Eriksson
   \date   2012-02-12
*/
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
       const Numeric&    aa )
{
  assert( refellipsoid.nelem() == 2 );
  assert( refellipsoid[0] > 0 );
  assert( refellipsoid[1] >= 0 );
  assert( refellipsoid[1] < 1 );
  assert( r > 0 );
  assert( za >= 90 );

  if( refellipsoid[1] < 1e-7 )        // e=1e-7 corresponds to that polar radius
    {                                 // less than 1 um smaller than equatorial 
      Numeric   x, y, z, dx, dy, dz;   // one for the Earth

      poslos2cart( x, y, z, dx, dy, dz, r, lat, lon, za, aa );
   
      const Numeric ppc   = r * sin( DEG2RAD * abs(za) );
      const Numeric l_tan = sqrt( r*r - ppc*ppc );
   
      cart2sph( r_tan, lat_tan, lon_tan, x+dx*l_tan, y+dy*l_tan, z+dz*l_tan );
    }

  else
    {
      // Equatorial and polar radii squared:
      const Numeric a2 = refellipsoid[0]*refellipsoid[0];
      const Numeric b2 = a2 * ( 1 - refellipsoid[1]*refellipsoid[1] ); 

      Vector X(3), xunit(3), yunit(3), zunit(3);

      poslos2cart( X[0], X[1], X[2], xunit[0], xunit[1], xunit[2], 
                                                         r, lat, lon, za, aa );
      cross( zunit, xunit, X );
      unitl( zunit );                // Normalisation of length to 1

      cross( yunit, zunit, xunit );
      unitl( yunit );                // Normalisation of length to 1

            Numeric x   = X[0];
            Numeric y   = X[1];
      const Numeric w11 = xunit[0];
      const Numeric w12 = yunit[0];
      const Numeric w21 = xunit[1];
      const Numeric w22 = yunit[1];
      const Numeric w31 = xunit[2];
      const Numeric w32 = yunit[2];

      const Numeric yr = X * yunit;
      const Numeric xr = X * xunit;

      const Numeric A = (w11*w11 + w21*w21)/a2 + w31*w31/b2;
      const Numeric B = 2.0*((w11*w12 + w21*w22)/a2 + (w31*w32)/b2);
      const Numeric C = (w12*w12 + w22*w22)/a2 + w32*w32/b2;

      if( B == 0.0 )
        { x = 0.0; }
      else 
        { 
          const Numeric K      = -2.0*A/B; 
          const Numeric factor = 1.0/(A+(B+C*K)*K);
          x = sqrt(factor);
          y = K*x;
        }

      const Numeric dist1 = (xr-X[0])*(xr-X[0]) + (yr-y)*(yr-y);
      const Numeric dist2 = (xr+X[0])*(xr+X[0]) + (yr+y)*(yr+y);
 	
      if( dist1 > dist2 )
        { x = -x; }

      cart2sph( r_tan, lat_tan, lon_tan, w11*x + w12*yr, w21*x + w22*yr,
                                                         w31*x + w32*yr );
    }
}
*/



//! line_sphere_intersect
/*! 
   Find the intersection between a line and a sphere

   \param   x     Out: X-coordinate of intersection
   \param   y     Out: Y-coordinate of intersection
   \param   z     Out: Z-coordinate of intersection
   \param   xl    A x-coordinate on the line
   \param   yl    A y-coordinate on the line
   \param   zl    A z-coordinate on the line
   \param   dx    X-component of line direction vector
   \param   dy    Y-component of line direction vector
   \param   dz    Z-component of line direction vector
   \param   xc    X-coordinate of sphere origo
   \param   yc    Y-coordinate of sphere origo
   \param   zc    Z-coordinate of sphere origo
   \param   r     Radius of sphere

   \author Patrick Eriksson
   \date   2012-03-30
*/
void line_sphere_intersect(
         Numeric&   x,
         Numeric&   y,
         Numeric&   z,
   const Numeric&   xl,
   const Numeric&   yl,
   const Numeric&   zl,
   const Numeric&   dx,
   const Numeric&   dy,
   const Numeric&   dz,
   const Numeric&   xc,
   const Numeric&   yc,
   const Numeric&   zc,
   const Numeric&   r )
{
  const Numeric a = dx*dx + dy*dy + dz*dz;
  const Numeric b = 2*( dx*(xl-xc) + dy*(yl-yc) + dz*(zl-zc) );
  const Numeric c = xc*xc + yc*yc + zc*zc + 
                    xl*xl + yl*yl + zl*zl - 2*(xc*xl + yc*yl + zc*zl) - r*r;
  
  Numeric d = b*b - 4*a*c;
  assert( d > 0 );
  
  const Numeric a2 = 2*a;
  const Numeric b2 = -b / a2;
  const Numeric e  = sqrt( d ) / a2;

  const Numeric l1 = b2 + e;
  const Numeric l2 = b2 - e;

  Numeric l;
  if( l1 < 0 )
    { l = l2; }
  else  if( l2 < 0 )
    { l = l1; }
  else
    { l = min(l1,l2); assert( l>=0 ); }

  x = xl + l*dx;
  y = yl + l*dy;
  z = zl + l*dz;
}



//! latlon_at_aa
/*! 
   Destination point given distance and bearing from start point

   Calculates the latitude and longitide, given a position, an azimuth angle
   and an angular distance.

   \param   lat2  Latitude of end position.
   \param   lon2  Longitude of end position.
   \param   lat1  Latitude of start position.
   \param   lon1  Longitude of start position.
   \param   aa    Azimuth/bearing at start point.
   \param   ddeg  Angular distance (in degrees).

   \author Patrick Eriksson
   \date   2012-04-23
*/
void latlon_at_aa(
         Numeric&   lat2,
         Numeric&   lon2,
   const Numeric&   lat1,
   const Numeric&   lon1,
   const Numeric&   aa,
   const Numeric&   ddeg )
{
  // Code from http://www.movable-type.co.uk/scripts/latlong.html
  // (but with short-cuts, such as asin(sin(lat2)) = lat2)
  // Note that lat1 here is another latitude
  
  const Numeric dang   = DEG2RAD * ddeg;
  const Numeric cosdang= cos( dang );
  const Numeric sindang= sin( dang );
  const Numeric latrad = DEG2RAD * lat1;
  const Numeric coslat = cos( latrad );
  const Numeric sinlat = sin( latrad );
  const Numeric aarad  = DEG2RAD * aa;

  lat2   = sinlat*cosdang + coslat*sindang*cos(aarad);
  lon2   = lon1 + RAD2DEG*( atan2( sin(aarad)*sindang*coslat,
                                   cosdang-sinlat*lat2 ) );
  lat2 = RAD2DEG * asin( lat2 );
}



//! los2xyz
/*! 
   Line-of-sight to another position given in cartesian coordinates.

   Calculates the zenith and azimuth angle for the geomrical path from 
   position 1 to position 2.

   \param   za    Out: LOS zenith angle at position 1.
   \param   aa    Out: LOS azimuth angle at position 1.
   \param   r     Radius of position 1.
   \param   lat   Latitude of position 1.
   \param   lon   Longitude of position 1.
   \param   x1    x-coordinate of position 1.
   \param   y1    y-coordinate of position 1.
   \param   z1    z-coordinate of position 1.
   \param   x2    x-coordinate of position 2.
   \param   y2    y-coordinate of position 2.
   \param   z2    z-coordinate of position 2.

   \author Patrick Eriksson
   \date   2012-03-26
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
   const Numeric&   z2 )
{
  Numeric dx = x2-x1, dy = y2-y1, dz = z2-z1;
  const Numeric ldxyz = sqrt( dx*dx + dy*dy + dz*dz );
  dx /= ldxyz;   
  dy /= ldxyz;   
  dz /= ldxyz;

  // All below extracted from 3D version of cart2poslos:
  const Numeric   latrad = DEG2RAD * lat1;
  const Numeric   lonrad = DEG2RAD * lon1;
  const Numeric   coslat = cos( latrad );
  const Numeric   sinlat = sin( latrad );
  const Numeric   coslon = cos( lonrad );
  const Numeric   sinlon = sin( lonrad );

  const Numeric   dr     = coslat*coslon*dx    + coslat*sinlon*dy    + sinlat*dz;
  const Numeric   dlat   = -sinlat*coslon/r1*dx - sinlat*sinlon/r1*dy + 
                                                                   coslat/r1*dz;
  const Numeric   dlon   = -sinlon/coslat/r1*dx + coslon/coslat/r1*dy;

  za = RAD2DEG * acos( dr );
  aa = RAD2DEG * acos( r1 * dlat / sin( DEG2RAD * za ) );
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
             Numeric&   x,
             Numeric&   y,
             Numeric&   z,
             Numeric&   dx,
             Numeric&   dy,
             Numeric&   dz,
       const Numeric&   r,
       const Numeric&   lat,
       const Numeric&   lon,
       const Numeric&   za,
       const Numeric&   aa )
{
  assert( r > 0 );
  assert( abs( lat ) <= 90 );
  assert( abs( lon ) <= 360 );
  assert( za >= 0 && za<=180 );

  // lat = +-90 
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  if( abs( lat ) > POLELAT )
    {
      const Numeric   s = sign( lat );

      x = 0;
      y = 0;
      z = s * r;

      dz = s * cos( DEG2RAD * za );
      dx = sin( DEG2RAD * za );
      dy = dx * sin( DEG2RAD * aa );
      dx = dx * cos( DEG2RAD * aa );
    }

  else
    {
      const Numeric   latrad = DEG2RAD * lat;
      const Numeric   lonrad = DEG2RAD * lon;
      const Numeric   zarad  = DEG2RAD * za;
      const Numeric   aarad  = DEG2RAD * aa;

      const Numeric   coslat = cos( latrad );
      const Numeric   sinlat = sin( latrad );
      const Numeric   coslon = cos( lonrad );
      const Numeric   sinlon = sin( lonrad );
      const Numeric   cosza  = cos( zarad );
      const Numeric   sinza  = sin( zarad );
      const Numeric   cosaa  = cos( aarad );
      const Numeric   sinaa  = sin( aarad );

      // This part as sph2cart but uses local variables
      x = r * coslat;   // Common term for x and y
      y = x * sinlon;
      x = x * coslon;
      z = r * sinlat;

      const Numeric   dr   = cosza;
      const Numeric   dlat = sinza * cosaa;         // r-term cancel out below
      const Numeric   dlon = sinza * sinaa / coslat; 

      dx = coslat*coslon * dr - sinlat*coslon * dlat - coslat*sinlon * dlon;
      dz =        sinlat * dr +        coslat * dlat;
      dy = coslat*sinlon * dr - sinlat*sinlon * dlat + coslat*coslon * dlon;
    }
}



//! pos2refell_r
/*!
    Extracts the reference ellipsoid for the specified position.

    The rference ellipsoid radius is defined slightly differently inside and
    outside of the model atmosphere. See *refell2r*. This function takes care
    of all those aspects. 

    If you know the latitude grid position of *rte_pos*, it is faster to use
    *refell2d*. 

    \return                 Ellispoid radius
    \param  atmosphere_dim  In: As the WSV with same name.
    \param  refellipsoid    In: As the WSV with same name.
    \param  lat_grid        In: As the WSV with same name.
    \param  lon_grid        In: As the WSV with same name.
    \param  rte_pos         In: As the WSV with same name.

    \author Patrick Eriksson 
    \date   2012-03-27
*/
Numeric pos2refell_r(
       const Index&     atmosphere_dim,
       ConstVectorView  refellipsoid,
       ConstVectorView  lat_grid,
       ConstVectorView  lon_grid,
       ConstVectorView  rte_pos )
{
  if( atmosphere_dim == 1 )
    { return refellipsoid[0]; }
  else
    {
      assert( rte_pos.nelem() > 1 );

      bool inside = true;

      if( rte_pos[1] < lat_grid[0] ||  rte_pos[1] > last(lat_grid) )
        { inside = false; }
      else if( atmosphere_dim == 3 )
        {
          assert( rte_pos.nelem() == 3 );
          if( rte_pos[2] < lon_grid[0] ||  rte_pos[2] > last(lon_grid) )
            { inside = false; }
        }
          
      if( inside )
        {
          GridPos gp_lat;
          gridpos( gp_lat, lat_grid, rte_pos[1] );
          return refell2d( refellipsoid, lat_grid, gp_lat );
        }
      else
        { return refell2r( refellipsoid, rte_pos[1] ); }
    }
}



//! refell2r
/*!
    Reference ellipsoid radius, directly from *refellipsoid*.

    Gives the distance from the Earth's centre and the reference ellipsoid as a
    function of geoCENTRIC latitude.

    For 1D, extract r directly as refellipsoid[0] to save time.

    This is the basic function to calculate the reference ellipsoid radius.
    However, inside the atmosphere this radius is just used at the positions of
    the lat_grid. A linear interpolation is applied between these points. This
    is handled by other functions. For 2D and 3D and the grid position is
    known, use *refell2d*. The function pos2refell_r handles all this in a
    general way (but not always the fastest option).

    \return                 Ellispoid radius
    \param  refellipsoid    In: As the WSV with same name.
    \param  latitude        In: A geoecentric latitude.

    \author Patrick Eriksson 
    \date   2012-02-07
*/
Numeric refell2r(
       ConstVectorView  refellipsoid,
       const Numeric&   lat )
{
  assert( refellipsoid.nelem() == 2 );
  assert( refellipsoid[0] > 0 );
  assert( refellipsoid[1] >= 0 );
  assert( refellipsoid[1] < 1 );

  if( refellipsoid[1] < 1e-7 )  // e=1e-7 corresponds to that polar radius  
    {                           // less than 1 um smaller than equatorial 
      return refellipsoid[0];   // one for the Earth
    }

  else
    {
      const Numeric   c = 1 - refellipsoid[1]*refellipsoid[1];
      const Numeric   b = refellipsoid[0] * sqrt( c );
      const Numeric   v = DEG2RAD * lat;
      const Numeric   ct = cos( v );
      const Numeric   st = sin( v );
      
      return b / sqrt( c*ct*ct + st*st );
    }
}



//! refell2d
/*!
    Reference ellipsoid radius for points inside 2D and 3D atmospheres.

    To be consistent with the ppath calculations, the ellipsoid radius shall be
    treated to vary linear between the latitude grid points. This function
    performs this operation. The latitude position is specified by its grid
    position (*gp*).

    \return                 Ellispoid radius
    \param  refellipsoid    In: As the WSV with same name.
    \param  lat_grid        In: As the WSV with same name.
    \param  gp              In: Latitude grid position.

    \author Patrick Eriksson 
    \date   2012-02-09
*/
Numeric refell2d(
       ConstVectorView  refellipsoid,
       ConstVectorView  lat_grid,
       const GridPos    gp )
{
  if( gp.fd[0] == 0 )
    return refell2r(refellipsoid,lat_grid[gp.idx]);
  else if( gp.fd[0] == 1 )
    return refell2r(refellipsoid,lat_grid[gp.idx+1]);
  else
    return gp.fd[1] * refell2r(refellipsoid,lat_grid[gp.idx]) +
           gp.fd[0] * refell2r(refellipsoid,lat_grid[gp.idx+1]);
}       



//! sphdist
/*!
    The distance between two geograpgical positions

    "As-the-crow-flies" distance between two points, specified by their
    latitude and longitude. 

    \return        Angular distance
    \param  lat1   Latitude of position 1.
    \param  lon1   Longitude of position 1.
    \param  lat2   Latitude of position 2.
    \param  lon2   Longitude of position 2.

    \author Patrick Eriksson 
    \date   2012-04-05
*/
Numeric sphdist(
   const Numeric&   lat1,
   const Numeric&   lon1,
   const Numeric&   lat2,
   const Numeric&   lon2 )
{
  // Equations taken from http://www.movable-type.co.uk/scripts/latlong.html
  const Numeric slat = sin( DEG2RAD*(lat2-lat1) / 2.0 );
  const Numeric slon = sin( DEG2RAD*(lon2-lon1) / 2.0 );
  const Numeric a = slat*slat + cos(DEG2RAD*lat1)*cos(DEG2RAD*lat2)*slon*slon;

  return RAD2DEG * 2 * atan2( sqrt(a), sqrt(1-a) );
}       





//! sph2cart
/*! 
   Conversion from spherical to cartesian coordinates.

   The cartesian coordinate system is defined such as the x-axis goes along
   lat=0 and lon=0, the z-axis goes along lat=0 and lon=90, and z-axis goes
   along lat=90. 

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
            Numeric&   x,
            Numeric&   y,
            Numeric&   z,
      const Numeric&   r,
      const Numeric&   lat,
      const Numeric&   lon )
{
  assert( r > 0 );
  assert( abs( lat ) <= 90 );
  assert( abs( lon ) <= 360 );

  const Numeric   latrad = DEG2RAD * lat;
  const Numeric   lonrad = DEG2RAD * lon;

  x = r * cos( latrad );   // Common term for x and z
  y = x * sin( lonrad );
  x = x * cos( lonrad );
  z = r * sin( latrad );
}





/*===========================================================================
  === coordinate transformations
  ===========================================================================*/

//! lon_shiftgrid
/*! 
   Shifting longitude grid by +/- 360 to the same region as a given longitude.

   Longitudes are allowed to be in [-360,360], but only covering at maximum
   360 degrees at once. Different variables might be specified in different
   spaces of the allowed range. lon_shiftgrid shifts the longitude grid to the
   same space as the lon value. However, no check is done that lon is indeed
   within the range of longrid_out (only same region).

   \param   longrid_out Out: shifted longitude grid.
   \param   longrid_in  Original longitude grid.
   \param   lon         A longitude value.

   \author Jana Mendrok
   \date   2012-06-28
*/
void lon_shiftgrid(
            Vector&    longrid_out,
      ConstVectorView  longrid_in,
      const Numeric    lon )
{
    longrid_out = longrid_in;
    if (longrid_in[longrid_in.nelem()-1] >= lon+360.)
      longrid_out += -360.;
    else if (longrid_in[0] <= lon-360.)
      longrid_out += 360.;
}

