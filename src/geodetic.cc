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



//! cart2poslos
/*! 
   The inverse of *poslos2cart*.

   The azimuth angle is set to: <br> 
      0 when the zenith angle is 0 or 180.
      atan2(dy,dx) at the poles (lat = +- 90).

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
  const double   dr   = coslat*coslon*dx + sinlat*dz + coslat*sinlon*dy;
  const double   dlat = -sinlat*coslon/r*dx + coslat/r*dz -sinlat*sinlon/r*dy;
  const double   dlon = -sinlon/coslat/r*dx + coslon/coslat/r*dy;

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
    { aa = RAD2DEG * atan2( dy, dx ); }
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
  lat = RAD2DEG * asin( z / r );
  lon = RAD2DEG * atan2( y, x ); 
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
void geomtanpoint2d( 
             double&    r_tan,
             double&    lat_tan,
     ConstVectorView    refellipsoid,
       const double&    r,
       const double&    lat,
       const double&    za )
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
void geomtanpoint( 
             double&    r_tan,
             double&    lat_tan,
             double&    lon_tan,
     ConstVectorView    refellipsoid,
       const double&    r,
       const double&    lat,
       const double&    lon,
       const double&    za,
       const double&    aa )
{
  assert( refellipsoid.nelem() == 2 );
  assert( refellipsoid[0] > 0 );
  assert( refellipsoid[1] >= 0 );
  assert( refellipsoid[1] < 1 );
  assert( r > 0 );
  assert( za >= 90 );

  if( refellipsoid[1] < 1e-7 )        // e=1e-7 corresponds to that polar radius
    {                                 // less than 1 um smaller than equatorial 
      double   x, y, z, dx, dy, dz;   // one for the Earth

      poslos2cart( x, y, z, dx, dy, dz, r, lat, lon, za, aa );
   
      const double ppc   = r * sin( DEG2RAD * abs(za) );
      const double l_tan = sqrt( r*r - ppc*ppc );
   
      cart2sph( r_tan, lat_tan, lon_tan, x+dx*l_tan, y+dy*l_tan, z+dz*l_tan );
    }

  else
    {
      // Equatorial and polar radii squared:
      const double a2 = refellipsoid[0]*refellipsoid[0];
      const double b2 = a2 * ( 1 - refellipsoid[1]*refellipsoid[1] ); 

      Vector X(3), xunit(3), yunit(3), zunit(3);

      poslos2cart( X[0], X[1], X[2], xunit[0], xunit[1], xunit[2], 
                                                         r, lat, lon, za, aa );
      cross( zunit, xunit, X );
      unitl( zunit );                // Normalisation of length to 1

      cross( yunit, zunit, xunit );
      unitl( yunit );                // Normalisation of length to 1

            double x   = X[0];
            double y   = X[1];
      const double w11 = xunit[0];
      const double w12 = yunit[0];
      const double w21 = xunit[1];
      const double w22 = yunit[1];
      const double w31 = xunit[2];
      const double w32 = yunit[2];

      const double yr = X * yunit;
      const double xr = X * xunit;

      const double A = (w11*w11 + w21*w21)/a2 + w31*w31/b2;
      const double B = 2.0*((w11*w12 + w21*w22)/a2 + (w31*w32)/b2);
      const double C = (w12*w12 + w22*w22)/a2 + w32*w32/b2;

      if( B == 0.0 )
        { x = 0.0; }
      else 
        { 
          const double K      = -2.0*A/B; 
          const double factor = 1.0/(A+(B+C*K)*K);
          x = sqrt(factor);
          y = K*x;
        }

      const double dist1 = (xr-X[0])*(xr-X[0]) + (yr-y)*(yr-y);
      const double dist2 = (xr+X[0])*(xr+X[0]) + (yr+y)*(yr+y);
 	
      if( dist1 > dist2 )
        { x = -x; }

      cart2sph( r_tan, lat_tan, lon_tan, w11*x + w12*yr, w21*x + w22*yr,
                                                         w31*x + w32*yr );
    }
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
  assert( r > 0 );
  assert( abs( lat ) <= 90 );
  assert( abs( lon ) <= 360 );
  assert( za >= 0 && za<=180 );

  // lat = +-90 
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  if( abs( lat ) > POLELAT )
    {
      const double   s = sign( lat );

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
      const double   latrad = DEG2RAD * lat;
      const double   lonrad = DEG2RAD * lon;
      const double   zarad  = DEG2RAD * za;
      const double   aarad  = DEG2RAD * aa;

      const double   coslat = cos( latrad );
      const double   sinlat = sin( latrad );
      const double   coslon = cos( lonrad );
      const double   sinlon = sin( lonrad );
      const double   cosza  = cos( zarad );
      const double   sinza  = sin( zarad );
      const double   cosaa  = cos( aarad );
      const double   sinaa  = sin( aarad );

      // This part as sph2cart but uses local variables
      x = r * coslat;   // Common term for x and z
      y = x * sinlon;
      x = x * coslon;
      z = r * sinlat;

      const double   dr   = cosza;
      const double   dlat = sinza * cosaa;         // r-terms cancel out below
      const double   dlon = sinza * sinaa / coslat; 

      dx = coslat*coslon * dr - sinlat*coslon * dlat - coslat*sinlon * dlon;
      dz =        sinlat * dr +        coslat * dlat;
      dy = coslat*sinlon * dr - sinlat*sinlon * dlat + coslat*coslon * dlon;
    }
}



//! refell2r
/*!
    Reference geoid radius, directly from *refellipsoid*.

    Gives the distance from the Earth's centre and the reference ellipsoid
    as a function of geoCENTRIC latitude. 

    For 1D, extract r directly as refellipsoid[0] to save time.

    For 2D and 3D and the position is inside the atmosphere, use *refell2d* and
    *refell3d* for highest internal consistency.

    \return                 Ellispoid radius
    \param  refellipsoid    In: As the WSV with same name.
    \param  latitude        In: A geoecentric latitude.

    \author Patrick Eriksson 
    \date   2012-02-07
*/
double refell2r(
       ConstVectorView  refellipsoid,
       const double&   lat )
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
      const double   c = 1 - refellipsoid[1]*refellipsoid[1];
      const double   b = refellipsoid[0] * sqrt( c );
      const double   v = DEG2RAD * lat;
      const double   ct = cos( v );
      const double   st = sin( v );
      
      return b / sqrt( c*ct*ct + st*st );
    }
}



//! refell2d
/*!
    Reference ellipsoid radius for points inside 2D atmospheres.

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
double refell2d(
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
  y = x * sin( lonrad );
  x = x * cos( lonrad );
  z = r * sin( latrad );
}






