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

#include "geodetic.h"
#include <cmath>
#include <stdexcept>
#include "arts_conversions.h"
#include "math_funcs.h"
#include "ppath.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);

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
void cart2pol(Numeric& r,
              Numeric& lat,
              const Numeric& x,
              const Numeric& z,
              const Numeric& lat0,
              const Numeric& za0) {
  r = sqrt(x * x + z * z);

  // Zenith and nadir cases
  const Numeric absza = abs(za0);
  if (absza < ANGTOL || absza > 180 - ANGTOL) {
    lat = lat0;
  }

  else {  // Latitude inside [0,360]
    lat = RAD2DEG * atan2(z, x);
    // Shift with n*360 to get as close as possible to lat0
    lat = lat - 360.0 * Numeric(round((lat - lat0) / 360.0));
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
void cart2poslos(Numeric& r,
                 Numeric& lat,
                 Numeric& za,
                 const Numeric& x,
                 const Numeric& z,
                 const Numeric& dx,
                 const Numeric& dz,
                 const Numeric& ppc,
                 const Numeric& lat0,
                 const Numeric& za0) {
  r = sqrt(x * x + z * z);

  // Zenith and nadir cases
  const Numeric absza = abs(za0);
  if (absza < ANGTOL || absza > 180 - ANGTOL) {
    lat = lat0;
    za = za0;
  }

  else {
    lat = RAD2DEG * atan2(z, x);

    const Numeric latrad = DEG2RAD * lat;
    const Numeric coslat = cos(latrad);
    const Numeric sinlat = sin(latrad);
    const Numeric dr = coslat * dx + sinlat * dz;

    // Use ppc for max accuracy, but dr required to resolve if up-
    // and downward cases.

    // Another possibility to obtain (absolute value of) za is
    // RAD2DEG*acos(dr).
    // It is checked that the two ways give consistent results, but
    // occasionally deviate with 1e-4 deg (due to numerical issues).

    za = RAD2DEG * asin(ppc / r);
    if (za0 > 0) {
      if (std::isnan(za)) {
        za = 90;
      } else if (dr < 0) {
        za = 180.0 - za;
      }
    } else {
      if (std::isnan(za)) {
        za = -90;
      } else if (dr < 0) {
        za = -180.0 + za;
      } else {
        za = -za;
      }
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
void distance2D(Numeric& l,
                const Numeric& r1,
                const Numeric& lat1,
                const Numeric& r2,
                const Numeric& lat2) {
  ARTS_ASSERT(abs(lat2 - lat1) <= 180);

  Numeric x1, z1, x2, z2;
  pol2cart(x1, z1, r1, lat1);
  pol2cart(x2, z2, r2, lat2);

  const Numeric dx = x2 - x1;
  const Numeric dz = z2 - z1;
  l = sqrt(dx * dx + dz * dz);
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
  ARTS_ASSERT( refellipsoid.nelem() == 2 );
  ARTS_ASSERT( refellipsoid[0] > 0 );
  ARTS_ASSERT( refellipsoid[1] >= 0 );
  ARTS_ASSERT( refellipsoid[1] < 1 );
  ARTS_ASSERT( r > 0 );
  ARTS_ASSERT( za >= 90 );
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
      ARTS_ASSERT( 0 );  // To be implemented
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
void line_circle_intersect(Numeric& x,
                           Numeric& z,
                           const Numeric& xl,
                           const Numeric& zl,
                           const Numeric& dx,
                           const Numeric& dz,
                           const Numeric& xc,
                           const Numeric& zc,
                           const Numeric& r) {
  const Numeric a = dx * dx + dz * dz;
  const Numeric b = 2 * (dx * (xl - xc) + dz * (zl - zc));
  const Numeric c =
      xc * xc + zc * zc + xl * xl + zl * zl - 2 * (xc * xl + zc * zl) - r * r;

  Numeric d = b * b - 4 * a * c;
  ARTS_ASSERT(d > 0);

  const Numeric a2 = 2 * a;
  const Numeric b2 = -b / a2;
  const Numeric e = sqrt(d) / a2;

  const Numeric l1 = b2 + e;
  const Numeric l2 = b2 - e;

  Numeric l;
  if (l1 < 0) {
    l = l2;
  } else if (l2 < 0) {
    l = l1;
  } else {
    l = min(l1, l2);
    ARTS_ASSERT(l >= 0);
  }

  x = xl + l * dx;
  z = zl + l * dz;
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
void pol2cart(Numeric& x, Numeric& z, const Numeric& r, const Numeric& lat) {
  ARTS_ASSERT(r > 0);

  const Numeric latrad = DEG2RAD * lat;

  x = r * cos(latrad);
  z = r * sin(latrad);
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
void poslos2cart(Numeric& x,
                 Numeric& z,
                 Numeric& dx,
                 Numeric& dz,
                 const Numeric& r,
                 const Numeric& lat,
                 const Numeric& za) {
  ARTS_ASSERT(r > 0);
  ARTS_ASSERT(za >= -180 && za <= 180);

  const Numeric latrad = DEG2RAD * lat;
  const Numeric zarad = DEG2RAD * za;

  const Numeric coslat = cos(latrad);
  const Numeric sinlat = sin(latrad);
  const Numeric cosza = cos(zarad);
  const Numeric sinza = sin(zarad);

  // This part as pol2cart but uses local variables
  x = r * coslat;
  z = r * sinlat;

  const Numeric dr = cosza;
  const Numeric dlat = sinza;  // r-term cancel out below

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
                 const Numeric& aa0) {
  // Radius of new point
  r = sqrt(x * x + y * y + z * z);

  // Zenith and nadir cases
  if (za0 < ANGTOL || za0 > 180 - ANGTOL) {
    lat = lat0;
    lon = lon0;
    za = za0;
    aa = aa0;
  }

  else {
    lat = RAD2DEG * asin(z / r);
    lon = RAD2DEG * atan2(y, x);

    bool ns_case = false;
    bool lon_flip = false;

    // Make sure that lon is maintained for N-S cases (if not
    // starting on a pole)
    if ((abs(aa0) < ANGTOL || abs(180 - aa0) < ANGTOL) &&
        abs(lat0) <= POLELAT) {
      ns_case = true;
      // Check that not lon changed with 180 deg
      if (abs(abs(lon - lon0) - 180) < 5) {
        lon_flip = true;
        if (lon0 > 0) {
          lon = lon0 - 180;
        } else {
          lon = lon0 + 180;
        }
      } else {
        lon = lon0;
      }
    }

    const Numeric latrad = DEG2RAD * lat;
    const Numeric lonrad = DEG2RAD * lon;
    const Numeric coslat = cos(latrad);
    const Numeric sinlat = sin(latrad);
    const Numeric coslon = cos(lonrad);
    const Numeric sinlon = sin(lonrad);

    // Set za by ppc for max accuracy, but this does not resolve
    // za and 180-za. This was first resolved by dr, but using l and lmax was
    // found to be more stable.
    za = RAD2DEG * asin(ppc / r);

    // Correct and check za
    if (std::isnan(za)) {
      za = 90;
    }
    // If za0 > 90, then correct za could be 180-za. Resolved by checking if
    // the tangent point is passed or not
    if (za0 > 90) {
      const Numeric l =
          sqrt(pow(x - x0, 2.0) + pow(y - y0, 2.0) + pow(z - z0, 2.0));
      const Numeric r0 = sqrt(x0 * x0 + y0 * y0 + z0 * z0);
      const Numeric ltan = geompath_l_at_r(ppc, r0);
      if (l < ltan) {
        za = 180.0 - za;
      }
    }

    // For lat = +- 90 the azimuth angle gives the longitude along which
    // the LOS goes
    if (abs(lat) >= POLELAT) {
      aa = RAD2DEG * atan2(dy, dx);
    }

    // N-S cases, not starting at a pole
    else if (ns_case) {
      if (!lon_flip) {
        aa = aa0;
      } else {
        if (abs(aa0) < ANGTOL) {
          aa = 180;
        } else {
          aa = 0;
        }
      }
    }

    else {
      const Numeric dlat = -sinlat * coslon / r * dx -
                           sinlat * sinlon / r * dy + coslat / r * dz;
      const Numeric dlon = -sinlon / coslat / r * dx + coslon / coslat / r * dy;

      aa = RAD2DEG * acos(r * dlat / sin(DEG2RAD * za));

      if (std::isnan(aa)) {
        if (dlat >= 0) {
          aa = 0;
        } else {
          aa = 180;
        }
      } else if (dlon < 0) {
        aa = -aa;
      }
    }
  }
}
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
                       const Numeric& dz) {
  cart2sph_plain(r, lat, lon, x, y, z);


  const Numeric latrad = DEG2RAD * lat;
  const Numeric lonrad = DEG2RAD * lon;
  const Numeric coslat = cos(latrad);
  const Numeric sinlat = sin(latrad);
  const Numeric coslon = cos(lonrad);
  const Numeric sinlon = sin(lonrad);

  const Numeric dr = coslat*coslon*dx + sinlat*dz + coslat*sinlon*dy;
  const Numeric dlat = -sinlat*coslon/r*dx + coslat/r*dz - sinlat*sinlon/r*dy;
  const Numeric dlon = -sinlon/coslat/r*dx + coslon/coslat/r*dy;

  za = acos( dr );
  aa = RAD2DEG * acos( r * dlat / sin( za ) );
  za *= RAD2DEG;

  // Corrections of aa
  if (std::isnan(aa)) {
    if (dlat >= 0)
      aa = 0;
    else
      aa = 180;
  } else if (dlon < 0) {
      aa = -aa;
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
void cart2sph(Numeric& r,
              Numeric& lat,
              Numeric& lon,
              const Numeric& x,
              const Numeric& y,
              const Numeric& z,
              const Numeric& lat0,
              const Numeric& lon0,
              const Numeric& za0,
              const Numeric& aa0) {
  r = sqrt(x * x + y * y + z * z);

  // Zenith and nadir cases
  if (za0 < ANGTOL || za0 > 180 - ANGTOL) {
    lat = lat0;
    lon = lon0;
  }

  else {
    lat = RAD2DEG * asin(z / r);
    lon = RAD2DEG * atan2(y, x);

    // Make sure that lon is maintained for N-S cases (if not
    // starting on a pole)
    if ((abs(aa0) < ANGTOL || abs(180 - aa0) < ANGTOL) &&
        abs(lat0) <= POLELAT) {
      // Check that not lon changed with 180 deg
      if (abs(lon - lon0) < 1) {
        lon = lon0;
      } else {
        if (lon0 > 0) {
          lon = lon0 - 180;
        } else {
          lon = lon0 + 180;
        }
      }
    }
  }
}
void cart2sph_plain(Numeric& r,
                    Numeric& lat,
                    Numeric& lon,
                    const Numeric& x,
                    const Numeric& y,
                    const Numeric& z) {
  r = sqrt(x * x + y * y + z * z);
  lat = RAD2DEG * asin(z / r);
  lon = RAD2DEG * atan2(y, x);
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
void distance3D(Numeric& l,
                const Numeric& r1,
                const Numeric& lat1,
                const Numeric& lon1,
                const Numeric& r2,
                const Numeric& lat2,
                const Numeric& lon2) {
  Numeric x1, y1, z1, x2, y2, z2;
  sph2cart(x1, y1, z1, r1, lat1, lon1);
  sph2cart(x2, y2, z2, r2, lat2, lon2);

  const Numeric dx = x2 - x1;
  const Numeric dy = y2 - y1;
  const Numeric dz = z2 - z1;
  l = sqrt(dx * dx + dy * dy + dz * dz);
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
void geompath_tanpos_3d(Numeric& r_tan,
                        Numeric& lat_tan,
                        Numeric& lon_tan,
                        Numeric& l_tan,
                        const Numeric& r,
                        const Numeric& lat,
                        const Numeric& lon,
                        const Numeric& za,
                        const Numeric& aa,
                        const Numeric& ppc) {
  ARTS_ASSERT(za >= 90);
  ARTS_ASSERT(r >= ppc);

  Numeric x, y, z, dx, dy, dz;

  poslos2cart(x, y, z, dx, dy, dz, r, lat, lon, za, aa);

  l_tan = sqrt(r * r - ppc * ppc);

  cart2sph(r_tan,
           lat_tan,
           lon_tan,
           x + dx * l_tan,
           y + dy * l_tan,
           z + dz * l_tan,
           lat,
           lon,
           za,
           aa);
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
  ARTS_ASSERT( refellipsoid.nelem() == 2 );
  ARTS_ASSERT( refellipsoid[0] > 0 );
  ARTS_ASSERT( refellipsoid[1] >= 0 );
  ARTS_ASSERT( refellipsoid[1] < 1 );
  ARTS_ASSERT( r > 0 );
  ARTS_ASSERT( za >= 90 );

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

//! line_refellipsoid_intersect
/*! 
   Find the intersection between a line and the reference ellipsoid

   A negative distance is returned if there is no intersection.

   If you the zenith angle and the reference ellipsoid is spherical, use the
   function geompath_l_at_r instead.

   \param   l     Out: Distance to intersection
   \param   refell Vector defining reference ellipsoid
   \param   xl    A x-coordinate on the line
   \param   yl    A y-coordinate on the line
   \param   zl    A z-coordinate on the line
   \param   dx    X-component of line direction vector
   \param   dy    Y-component of line direction vector
   \param   dz    Z-component of line direction vector

   \author Patrick Eriksson
   \date   2012-03-30
*/
void line_refellipsoid_intersect(Numeric& l,
                                 const Vector& refellipsoid,
                                 const Numeric& x,
                                 const Numeric& y,
                                 const Numeric& z,
                                 const Numeric& dx,
                                 const Numeric& dy,
                                 const Numeric& dz) {
  // Code taken from Atmlab's ellipsoid_intersection

  // Spherical case
  if (refellipsoid[1] < 1e-7) {
    const Numeric p  = x*dx + y*dy + z*dz;
    const Numeric pp = p*p;
    const Numeric q = x*x + y*y + z*z - refellipsoid[0]*refellipsoid[0];
    if (q>pp)
      l = -1.0;
    else {
      const Numeric sq = sqrt(pp - q);
      if (-p > sq)
        l = -p - sq;
      else
        l = -p + sq;
    }
  }

  // Ellipsoid case
  else {
    // Based on https://medium.com/@stephenhartzell/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6
    const Numeric a = refellipsoid[0];
    const Numeric b = refellipsoid[0];
    const Numeric c = refellipsoid[0] * sqrt(1-refellipsoid[1]*refellipsoid[1]);
    const Numeric a2 = a*a;
    const Numeric b2 = b*b;
    const Numeric c2 = c*c;
    const Numeric x2 = x*x;
    const Numeric y2 = y*y;
    const Numeric z2 = z*z;
    const Numeric dx2 = dx*dx;
    const Numeric dy2 = dy*dy;
    const Numeric dz2 = dz*dz;
    const Numeric rad = a2*b2*dz2 + a2*c2*dy2 - a2*dy2*z2 + 2*a2*dy*dz*y*z -
                        a2*dz2*y2 + b2*c2*dx2 - b2*dx2*z2 + 2*b2*dx*dz*x*z -
                        b2*dz2*x2 - c2*dx2*y2 + 2*c2*dx*dy*x*y - c2*dy2*x2;
    if (rad<0)
      l = -1.0;
    else {
      const Numeric val = -a2*b2*dz*z - a2*c2*dy*y - b2*c2*dx*x;
      const Numeric mag = a2*b2*dz2 + a2*c2*dy2 + b2*c2*dx2;
      const Numeric abc = a*b*c*sqrt(rad);
      if (val > abc)
        l = (val - abc) / mag;
      else
        l = (val + abc) / mag;
    }
  }
}

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
                           const Numeric& r) {
  const Numeric a = dx * dx + dy * dy + dz * dz;
  const Numeric b = 2 * (dx * (xl - xc) + dy * (yl - yc) + dz * (zl - zc));
  const Numeric c = xc * xc + yc * yc + zc * zc + xl * xl + yl * yl + zl * zl -
                    2 * (xc * xl + yc * yl + zc * zl) - r * r;

  Numeric d = b * b - 4 * a * c;
  ARTS_ASSERT(d > 0);

  const Numeric a2 = 2 * a;
  const Numeric b2 = -b / a2;
  const Numeric e = sqrt(d) / a2;

  const Numeric l1 = b2 + e;
  const Numeric l2 = b2 - e;

  Numeric l;
  if (l1 < 0) {
    l = l2;
  } else if (l2 < 0) {
    l = l1;
  } else {
    l = min(l1, l2);
    ARTS_ASSERT(l >= 0);
  }

  x = xl + l * dx;
  y = yl + l * dy;
  z = zl + l * dz;
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
void latlon_at_aa(Numeric& lat2,
                  Numeric& lon2,
                  const Numeric& lat1,
                  const Numeric& lon1,
                  const Numeric& aa,
                  const Numeric& ddeg) {
  // Code from http://www.movable-type.co.uk/scripts/latlong.html
  // (but with short-cuts, such as asin(sin(lat2)) = lat2)
  // Note that lat1 here is another latitude

  const Numeric dang = DEG2RAD * ddeg;
  const Numeric cosdang = cos(dang);
  const Numeric sindang = sin(dang);
  const Numeric latrad = DEG2RAD * lat1;
  const Numeric coslat = cos(latrad);
  const Numeric sinlat = sin(latrad);
  const Numeric aarad = DEG2RAD * aa;

  lat2 = sinlat * cosdang + coslat * sindang * cos(aarad);
  lon2 = lon1 + RAD2DEG * (atan2(sin(aarad) * sindang * coslat,
                                 cosdang - sinlat * lat2));
  lat2 = RAD2DEG * asin(lat2);
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
             const Numeric& z2) {
  Numeric dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
  const Numeric ldxyz = sqrt(dx * dx + dy * dy + dz * dz);
  dx /= ldxyz;
  dy /= ldxyz;
  dz /= ldxyz;

  // All below extracted from 3D version of cart2poslos:
  const Numeric latrad = DEG2RAD * lat1;
  const Numeric lonrad = DEG2RAD * lon1;
  const Numeric coslat = cos(latrad);
  const Numeric sinlat = sin(latrad);
  const Numeric coslon = cos(lonrad);
  const Numeric sinlon = sin(lonrad);

  const Numeric dr = coslat * coslon * dx + coslat * sinlon * dy + sinlat * dz;
  const Numeric dlat =
      -sinlat * coslon / r1 * dx - sinlat * sinlon / r1 * dy + coslat / r1 * dz;
  const Numeric dlon = -sinlon / coslat / r1 * dx + coslon / coslat / r1 * dy;

  za = RAD2DEG * acos(dr);
  aa = RAD2DEG * acos(r1 * dlat / sin(DEG2RAD * za));
  if (std::isnan(aa)) {
    if (dlat >= 0) {
      aa = 0;
    } else {
      aa = 180;
    }
  } else if (dlon < 0) {
    aa = -aa;
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
                 const Numeric& aa) {
  ARTS_ASSERT(r > 0);
  ARTS_ASSERT(abs(lat) <= 90);
  //ARTS_ASSERT( abs( lon ) <= 360 );
  ARTS_ASSERT(za >= 0 && za <= 180);

  // lat = +-90
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  if (abs(lat) > POLELAT) {
    const Numeric s = sign(lat);

    x = 0;
    y = 0;
    z = s * r;

    dz = s * cos(DEG2RAD * za);
    dx = sin(DEG2RAD * za);
    dy = dx * sin(DEG2RAD * aa);
    dx = dx * cos(DEG2RAD * aa);
  }

  else {
    const Numeric latrad = DEG2RAD * lat;
    const Numeric lonrad = DEG2RAD * lon;
    const Numeric zarad = DEG2RAD * za;
    const Numeric aarad = DEG2RAD * aa;

    const Numeric coslat = cos(latrad);
    const Numeric sinlat = sin(latrad);
    const Numeric coslon = cos(lonrad);
    const Numeric sinlon = sin(lonrad);
    const Numeric cosza = cos(zarad);
    const Numeric sinza = sin(zarad);
    const Numeric cosaa = cos(aarad);
    const Numeric sinaa = sin(aarad);

    // This part as sph2cart but uses local variables
    x = r * coslat;  // Common term for x and y
    y = x * sinlon;
    x = x * coslon;
    z = r * sinlat;

    const Numeric dr = cosza;
    const Numeric dlat = sinza * cosaa;  // r-term cancel out below
    const Numeric dlon = sinza * sinaa / coslat;

    dx = coslat * coslon * dr - sinlat * coslon * dlat - coslat * sinlon * dlon;
    dz = sinlat * dr + coslat * dlat;
    dy = coslat * sinlon * dr - sinlat * sinlon * dlat + coslat * coslon * dlon;
  }
}

//! pos2refell_r
/*!
    Extracts the reference ellipsoid for the specified position.

    The reference ellipsoid radius is defined slightly differently inside and
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
Numeric pos2refell_r(const Index& atmosphere_dim,
                     ConstVectorView refellipsoid,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstVectorView rte_pos) {
  if (atmosphere_dim == 1) {
    return refellipsoid[0];
  } else {
    ARTS_ASSERT(rte_pos.nelem() > 1);

    bool inside = true;

    if (rte_pos[1] < lat_grid[0] || rte_pos[1] > last(lat_grid)) {
      inside = false;
    } else if (atmosphere_dim == 3) {
      ARTS_ASSERT(rte_pos.nelem() == 3);
      if (rte_pos[2] < lon_grid[0] || rte_pos[2] > last(lon_grid)) {
        inside = false;
      }
    }

    if (inside) {
      GridPos gp_lat;
      gridpos(gp_lat, lat_grid, rte_pos[1]);
      return refell2d(refellipsoid, lat_grid, gp_lat);
    } else {
      return refell2r(refellipsoid, rte_pos[1]);
    }
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
Numeric refell2r(ConstVectorView refellipsoid, const Numeric& lat) {
  ARTS_ASSERT(refellipsoid.nelem() == 2);
  ARTS_ASSERT(refellipsoid[0] > 0);
  ARTS_ASSERT(refellipsoid[1] >= 0);
  ARTS_ASSERT(refellipsoid[1] < 1);

  if (refellipsoid[1] < 1e-7)  // e=1e-7 corresponds to that polar radius
  {                            // less than 1 um smaller than equatorial
    return refellipsoid[0];    // one for the Earth
  }

  else {
    const Numeric c = 1 - refellipsoid[1] * refellipsoid[1];
    const Numeric b = refellipsoid[0] * sqrt(c);
    const Numeric v = DEG2RAD * lat;
    const Numeric ct = cos(v);
    const Numeric st = sin(v);

    return b / sqrt(c * ct * ct + st * st);
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
Numeric refell2d(ConstVectorView refellipsoid,
                 ConstVectorView lat_grid,
                 const GridPos gp) {
  if (gp.fd[0] == 0)
    return refell2r(refellipsoid, lat_grid[gp.idx]);
  else if (gp.fd[0] == 1)
    return refell2r(refellipsoid, lat_grid[gp.idx + 1]);
  else
    return gp.fd[1] * refell2r(refellipsoid, lat_grid[gp.idx]) +
           gp.fd[0] * refell2r(refellipsoid, lat_grid[gp.idx + 1]);
}

//! sphdist
/*!
    The distance between two geograpgical positions

    "As-the-crow-flies" angular distance between two points, specified by their
    latitude and longitude. 

    Note that angular distance is returned. The distance in length varies with
    altitude.

    \return        Angular distance (in degrees)
    \param  lat1   Latitude of position 1.
    \param  lon1   Longitude of position 1.
    \param  lat2   Latitude of position 2.
    \param  lon2   Longitude of position 2.

    \author Patrick Eriksson 
    \date   2012-04-05
*/
Numeric sphdist(const Numeric& lat1,
                const Numeric& lon1,
                const Numeric& lat2,
                const Numeric& lon2) {
  // Equations taken from http://www.movable-type.co.uk/scripts/latlong.html
  const Numeric slat = sin(DEG2RAD * (lat2 - lat1) / 2.0);
  const Numeric slon = sin(DEG2RAD * (lon2 - lon1) / 2.0);
  const Numeric a =
      slat * slat + cos(DEG2RAD * lat1) * cos(DEG2RAD * lat2) * slon * slon;

  return RAD2DEG * 2 * atan2(sqrt(a), sqrt(1 - a));
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
void sph2cart(Numeric& x,
              Numeric& y,
              Numeric& z,
              const Numeric& r,
              const Numeric& lat,
              const Numeric& lon) {
  ARTS_ASSERT(r > 0);
  ARTS_ASSERT(abs(lat) <= 90);
  ARTS_ASSERT(abs(lon) <= 360);

  const Numeric latrad = DEG2RAD * lat;
  const Numeric lonrad = DEG2RAD * lon;

  x = r * cos(latrad);  // Common term for x and z
  y = x * sin(lonrad);
  x = x * cos(lonrad);
  z = r * sin(latrad);
}



/*===========================================================================
  === Fixes for longitudes
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
void lon_shiftgrid(Vector& longrid_out,
                   ConstVectorView longrid_in,
                   const Numeric lon) {
  longrid_out = longrid_in;
  if (longrid_in[longrid_in.nelem() - 1] >= lon + 360.)
    longrid_out += -360.;
  else if (longrid_in[0] <= lon - 360.)
    longrid_out += 360.;
}

//! Cyclic latitude longitude coordinates.
/* This function wraps around the given latitude and longitude coordinates
 *  to the ranges
 * - lat in [-90.0, 90.0]
 * - lon in [0.0, 360.0]
 * Only in the range [-180.0, 180.0] are accepted otherwise an error will
 * be thrown
 *
 * \param[in, out] lat The latitude coordinate.
 * \param[in, out] lon The longitude coordiante.
 *
 * \author Simon Pfreundschuh
 * \date   2018-05-22
*/
void cycle_lat_lon(Numeric& lat, Numeric& lon) {
  ARTS_USER_ERROR_IF (lat < -180.0,
                      "Latitude values < -180.0 are not supported.");
  ARTS_USER_ERROR_IF (lat > 180.0,
                      "Latitude values > 180.0 are not supported.");

  if (lat < -90.0) {
    lat = -180.0 - lat;
    lon += 180.0;
  }
  if (lat > 90.0) {
    lat = 180.0 - lat;
    lon += 180.0;
  }

  while (lon < 0.0) {
    lon += 360.0;
  }
  while (lon > 360.0) {
    lon -= 360.0;
  }
}



/*===========================================================================
  === Functions involving geodetic latitudes (based on functions in Atmlab)
  ===========================================================================*/

//! Conversion from cartesian to geodetic coordinates
/* 
 * \param[out]  h   Geodetic altitude
 * \param[out]  la  Geodetic latitude
 * \param[out]  lon Longitude
 * \param[in] x   x-coordinate (ECEF)
 * \param[in] y   y-coordinate (ECEF)
 * \param[in] z   z-coordinate (ECEF)
 * \param[in]  refellipsoid As the WSV with the same name.
 *
 * \author Patrick Eriksson
 * \date   2020-09-17
*/
void cart2geodetic(Numeric& h,
                   Numeric& lat,
                   Numeric& lon,
                   const Numeric& x,
                   const Numeric& y,
                   const Numeric& z,
                   const Vector& refellipsoid ) {
  // Use geocentric function if geoid is spherical
  if (refellipsoid[1] < 1e-7)
    { Numeric r;
      cart2sph_plain(r, lat, lon, x, y, z);
      h = r - refellipsoid[0];
    }
  else
    {
      lon = RAD2DEG * atan2(y,x);

      const Numeric sq = sqrt(x*x+y*y);
      Numeric B0 = atan2(z,sq);
      Numeric B = B0-1, N;
      const Numeric e2 = refellipsoid[1]*refellipsoid[1];
      
      while (abs(B-B0)>1e-10) {
        N = refellipsoid[0] / sqrt(1-e2*sin(B0)*sin(B0));
        h = sq / cos(B0) - N;
        B = B0;
        B0 = atan((z/sq) * 1/(1-e2*N/(N+h)));
      }
      lat = RAD2DEG * B;
    }
}



//! Conversion from geodetic to cartesian coordinates
/* 
 * \param[out] x   x-coordinate (ECEF)
 * \param[out] y   y-coordinate (ECEF)
 * \param[out] z   z-coordinate (ECEF)
 * \param[in]  h   Geodetic altitude
 * \param[in]  lat Geodetic latitude
 * \param[in]  lon Longitude
 * \param[in]  refellipsoid As the WSV with the same name.
 *
 * \author Patrick Eriksson
 * \date   2020-09-17
*/
void geodetic2cart(Numeric& x,
                   Numeric& y,
                   Numeric& z,
                   const Numeric& h,
                   const Numeric& lat,
                   const Numeric& lon,
                   const Vector& refellipsoid ) {
  // Use geocentric function if geoid is spherical
  if (refellipsoid[1] < 1e-7)
    { sph2cart(x, y, z, h+refellipsoid[0], lat, lon); }
  else
    {
      const Numeric a = refellipsoid[0];
      const Numeric e2 = refellipsoid[1]*refellipsoid[1];
      const Numeric sinlat = sin(DEG2RAD*lat);
      const Numeric coslat = cos(DEG2RAD*lat);
      const Numeric N = a / sqrt(1 - e2*sinlat*sinlat);

      x = (N + h) * coslat * cos(DEG2RAD*lon);
      y = (N + h) * coslat * sin(DEG2RAD*lon);
      z = (N*(1 - e2) + h) * sinlat;
    }
}



//! geodeticposlos2cart
/*! 
   As *poslos2cart* but starts with geodetic position and LOS.

   \param   x     Out: x-coordinate of observation position.
   \param   y     Out: y-coordinate of observation position.
   \param   z     Out: z-coordinate of observation position.
   \param   dx    Out: x-part of LOS unit vector.
   \param   dy    Out: y-part of LOS unit vector.
   \param   dz    Out: z-part of LOS unit vector.
   \param   h     Geodetic altitude of observation position.
   \param   lat   Geodetic latitude of observation position.
   \param   lon   Longitude of observation position.
   \param   za    LOS zenith angle at observation position.
   \param   aa    LOS azimuth angle at observation position.

   \author Patrick Eriksson
   \date   2020-09-17
*/
void geodeticposlos2cart(Numeric& x,
                         Numeric& y,
                         Numeric& z,
                         Numeric& dx,
                         Numeric& dy,
                         Numeric& dz,
                         const Numeric& h,
                         const Numeric& lat,
                         const Numeric& lon,
                         const Numeric& za,
                         const Numeric& aa,
                         const Vector& refellipsoid ) {
  ARTS_ASSERT(abs(lat) <= 90);
  ARTS_ASSERT(za >= 0 && za <= 180);

  // lat = +-90
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  // At the poles, no difference between geocentric and geodetic zenith
  if (abs(lat) > POLELAT) {
    const Numeric s = sign(lat);

    x = 0;
    y = 0;
    z = s * (h + refellipsoid[0]*sqrt(1-refellipsoid[1]*refellipsoid[1]));

    dz = s * cos(DEG2RAD * za);
    dx = sin(DEG2RAD * za);
    dy = dx * sin(DEG2RAD * aa);
    dx = dx * cos(DEG2RAD * aa);
  }

  else {      
    const Numeric coslat = cos(DEG2RAD * lat);
    const Numeric sinlat = sin(DEG2RAD * lat);
    const Numeric coslon = cos(DEG2RAD * lon);
    const Numeric sinlon = sin(DEG2RAD * lon);

    geodetic2cart(x, y, z, h, lat, lon, refellipsoid);

    Numeric de, dn, du;
    zaaa2enu(de, dn, du, za, aa);  

    // See
    // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
    dx = -sinlon*de - sinlat*coslon*dn + coslat*coslon*du;  
    dy =  coslon*de - sinlat*sinlon*dn + coslat*sinlon*du;  
    dz =              coslat*       dn + sinlat*       du;  
  }
}



//! cart2geodeticposlos
/*! 
   The inverse of *geodeticposlos2cart*.

   \param   h     Out: Geodetic altitude.
   \param   lat   Out: Geodetic latitude
   \param   lon   Out: Longitude
   \param   za    Out: LOS zenith angle
   \param   aa    Out: LOS azimuth angle
   \param   x     x-coordinate
   \param   y     y-coordinate
   \param   z     z-coordinate
   \param   dx    x-part of LOS unit vector.
   \param   dy    y-part of LOS unit vector.
   \param   dz    z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2020-09-17
*/
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
                         const Vector& refellipsoid ) {

  cart2geodetic(h, lat, lon, x, y, z, refellipsoid );
  
  const Numeric latrad = DEG2RAD * lat;
  const Numeric lonrad = DEG2RAD * lon;
  const Numeric coslat = cos(latrad);
  const Numeric sinlat = sin(latrad);
  const Numeric coslon = cos(lonrad);
  const Numeric sinlon = sin(lonrad);

  // See
  // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
  const Numeric de =        -sinlon*dx +        coslon*dy;  
  const Numeric dn = -sinlat*coslon*dx - sinlat*sinlon*dy + coslat*dz;  
  const Numeric du =  coslat*coslon*dx + coslat*sinlon*dy + sinlat*dz;  

  enu2zaaa(za, aa, de, dn, du);
}
