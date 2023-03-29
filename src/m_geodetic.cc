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
  \file   m_geodetic.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2012-02-06

  \brief  Workspace functions of geodetic character.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "arts_conversions.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "matpack_data.h"
#include "messages.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidForAzimuth(Vector& refellipsoid,
                            const Numeric& latitude,
                            const Numeric& azimuth,
                            const Verbosity&) {
  ARTS_USER_ERROR_IF (refellipsoid.nelem() != 2,
                      "Input *refellispoid must be a vector of length 2*.");

  if (refellipsoid[1] > 0) {
    const Numeric e2 = refellipsoid[1] * refellipsoid[1];
    const Numeric a = 1 - e2 * pow(sin(DEG2RAD * latitude), 2.0);

    const Numeric rn = 1 / sqrt(a);
    const Numeric rm = (1 - e2) * (rn / a);

    const Numeric v = DEG2RAD * azimuth;

    refellipsoid[0] =
        refellipsoid[0] / (pow(cos(v), 2.0) / rm + pow(sin(v), 2.0) / rn);
    refellipsoid[1] = 0;
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidOrbitPlane(Vector& refellipsoid,
                            const Numeric& orbitinc,
                            const Verbosity&) {
  ARTS_USER_ERROR_IF (refellipsoid.nelem() != 2,
                      "Input *refellispoid must be a vector of length 2*.");
  chk_if_in_range("orbitinc", orbitinc, 0, 180);

  // Radius at maximum latitude
  const Numeric rp = refell2r(refellipsoid, orbitinc);

  // New eccentricity
  refellipsoid[1] = sqrt(1 - pow(rp / refellipsoid[0], 2.0));
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidSet(Vector& refellipsoid,
                     const Numeric& re,
                     const Numeric& e,
                     const Verbosity&) {
  refellipsoid.resize(2);

  refellipsoid[0] = re;
  refellipsoid[1] = e;
}



/* Workspace method: Doxygen documentation will be auto-generated */
void rte_poslosFromECEF(Vector& rte_pos,
                        Vector& rte_los,
                        const Matrix& sensor_pos_ecef,
                        const Matrix& sensor_los_ecef,
                        const Vector& refellipsoid,
                        const Verbosity& v) {
  ARTS_USER_ERROR_IF (sensor_pos_ecef.nrows() != 1,
                      "For this WSM, *sensor_pos_ecef* can only have one row.");

  Matrix sensor_pos, sensor_los;
  sensor_poslosFromECEF(sensor_pos,
                        sensor_los,
                        sensor_pos_ecef,
                        sensor_los_ecef,
                        refellipsoid,
                        v );
  rte_pos = sensor_pos(0,joker);
  if(sensor_los.empty())
    rte_los.resize(0);
  else
    rte_los = sensor_los(0,joker);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void rte_poslosFromGeodetic(Vector& rte_pos,
                            Vector& rte_los,
                            const Matrix& sensor_pos_geodetic,
                            const Matrix& sensor_los_geodetic,
                            const Vector& refellipsoid,
                            const Verbosity& v) {
  ARTS_USER_ERROR_IF (sensor_pos_geodetic.nrows() != 1,
                      "For this WSM, *sensor_pos_geodetic* can only have one row.");

  Matrix sensor_pos, sensor_los;
  sensor_poslosFromGeodetic(sensor_pos,
                            sensor_los,
                            sensor_pos_geodetic,
                            sensor_los_geodetic,
                            refellipsoid,
                            v );
  rte_pos = sensor_pos(0,joker);
  if(sensor_los.empty())
    rte_los.resize(0);
  else
    rte_los = sensor_los(0,joker);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_poslosFromECEF(Matrix& sensor_pos,
                           Matrix& sensor_los,
                           const Matrix& sensor_pos_ecef,
                           const Matrix& sensor_los_ecef,
                           const Vector& refellipsoid,
                           const Verbosity&) {
  const Index ncols = sensor_pos_ecef.ncols();
  const Index nrows = sensor_pos_ecef.nrows();
  ARTS_USER_ERROR_IF (ncols != 3,
                      "*sensor_pos_geodetic* must have three columns.");
  
  sensor_pos.resize( nrows, ncols );

  if (sensor_los_ecef.empty()) {
    sensor_los = sensor_los_ecef;

    for (Index i=0; i<nrows; i++) {
      cart2sph_plain(sensor_pos(i,0),
                     sensor_pos(i,1),
                     sensor_pos(i,2),
                     sensor_pos_ecef(i,0),
                     sensor_pos_ecef(i,1),
                     sensor_pos_ecef(i,2));
      sensor_pos(i,0) -= refell2r(refellipsoid,sensor_pos(i,1));
    }
  }
  else {
    const Index ncols2 = sensor_los_ecef.ncols();
    ARTS_USER_ERROR_IF (ncols2 != 3,
                        "*sensor_los_ecef* must be empty or have "
                        "three columns.");
    ARTS_USER_ERROR_IF (nrows != sensor_los_ecef.nrows(),
                        "*sensor_los_ecef* must be empty or have the "
                        "same number of rows as *sensor_pos_ecef*");

    sensor_pos.resize( nrows, ncols );
    sensor_los.resize( nrows, 2 );
      
    for (Index i=0; i<nrows; i++) {
      cart2poslos_plain(sensor_pos(i,0),
                        sensor_pos(i,1),
                        sensor_pos(i,2),
                        sensor_los(i,0),
                        sensor_los(i,1),
                        sensor_pos_ecef(i,0),
                        sensor_pos_ecef(i,1),
                        sensor_pos_ecef(i,2),
                        sensor_los_ecef(i,0),
                        sensor_los_ecef(i,1),
                        sensor_los_ecef(i,2));                          
      sensor_pos(i,0) -= refell2r(refellipsoid,sensor_pos(i,1));
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_poslosFromGeodetic(Matrix& sensor_pos,
                               Matrix& sensor_los,
                               const Matrix& sensor_pos_geodetic,
                               const Matrix& sensor_los_geodetic,
                               const Vector& refellipsoid,
                               const Verbosity&) {
  const Index ncols = sensor_pos_geodetic.ncols();
  const Index nrows = sensor_pos_geodetic.nrows();
  ARTS_USER_ERROR_IF (ncols != 3,
                      "*sensor_pos_geodetic* must have three columns.");
  
  // No conversion to do if geoid is spherical
  if (refellipsoid[1] < 1e-7) {
    sensor_pos = sensor_pos_geodetic;
    sensor_los = sensor_los_geodetic;
  }
  else {
    sensor_pos.resize( nrows, ncols );

    if (sensor_los_geodetic.empty()) {
      sensor_los = sensor_los_geodetic;
      Numeric x, y, z;
      
      for (Index i=0; i<nrows; i++) {
        geodetic2cart(x, y, z,
                      sensor_pos_geodetic(i,0),
                      sensor_pos_geodetic(i,1),
                      sensor_pos_geodetic(i,2),
                      refellipsoid );
        cart2sph_plain(sensor_pos(i,0),
                       sensor_pos(i,1),
                       sensor_pos(i,2),
                       x, y, z);
        sensor_pos(i,0) -= refell2r(refellipsoid,sensor_pos(i,1));
      }
    } else {
      const Index ncols2 = sensor_los_geodetic.ncols();
      ARTS_USER_ERROR_IF (ncols2 != 2,
                          "*sensor_los_geodetic* must be empty or have "
                          "two columns.");
      ARTS_USER_ERROR_IF (nrows != sensor_los_geodetic.nrows(),
                          "*sensor_los_geodetic* must be empty or have the "
                          "same number of rows as *sensor_pos_geodetic*");

      sensor_pos.resize( nrows, ncols );
      sensor_los.resize( nrows, ncols2 );
      Numeric x, y, z, dx, dy, dz;
      
      for (Index i=0; i<nrows; i++) {
        geodeticposlos2cart(x, y, z, dx,dy, dz,
                            sensor_pos_geodetic(i,0),
                            sensor_pos_geodetic(i,1),
                            sensor_pos_geodetic(i,2),
                            sensor_los_geodetic(i,0),
                            sensor_los_geodetic(i,1),
                            refellipsoid );
        cart2poslos_plain(sensor_pos(i,0),
                          sensor_pos(i,1),
                          sensor_pos(i,2),
                          sensor_los(i,0),
                          sensor_los(i,1),
                          x, y, z, dx, dy, dz);
        sensor_pos(i,0) -= refell2r(refellipsoid,sensor_pos(i,1));
      }
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void y_geoToGeodetic(Matrix& y_geo,
                     const Vector& refellipsoid,
                     const Verbosity&) {
  if (y_geo.empty() ||  refellipsoid[1] < 1e-7 || std::isnan(y_geo(0,0))) {
    // Do nothing
  } else {
    for (Index i=0; i<y_geo.nrows(); i++) {
      Numeric x, y, z, dx, dy, dz;
      Numeric r = y_geo(i,0) + refell2r(refellipsoid,y_geo(i,1));
      poslos2cart(x, y, z, dx, dy, dz, r,
                  y_geo(i,1), y_geo(i,2), y_geo(i,3), y_geo(i,4));
      Numeric h, lat, lon, za, aa;
      cart2geodeticposlos(h, lat, lon, za, aa, x, y, z, dx, dy, dz,refellipsoid);
      y_geo(i,0) = h;
      y_geo(i,1) = lat;
      y_geo(i,2) = lon;
      y_geo(i,3) = za;
      y_geo(i,4) = aa;
    }
  }
}
  

/* Workspace method: Doxygen documentation will be auto-generated */
void IntersectionGeometricalWithAltitude(Matrix& pos,
                                         Matrix& los,
                                         const Matrix& sensor_pos,
                                         const Matrix& sensor_los,
                                         const Vector& refellipsoid,
                                         const Vector& lat_grid,
                                         const Vector& lon_grid,
                                         const Numeric& altitude,
                                         const Verbosity&) {
  //
  
  ARTS_USER_ERROR_IF (sensor_pos.ncols() != 3,
                      "*sensor_pos* must have three columns.");
  ARTS_USER_ERROR_IF (sensor_pos.nrows() == 0,
                      "*sensor_pos* must have at least one row.");
  for (Index i=0; i<sensor_pos.nrows(); i++) {
    ARTS_USER_ERROR_IF (sensor_pos(i,1) < -90 || sensor_pos(i,1) > 90,
                        "Unvalid latitude in *sensor_pos*.\n"
                        "Latitudes must be inside [-90,90],\n"
                        "but sensor_pos(",i,",1) is ", sensor_pos(i,1));
    ARTS_USER_ERROR_IF (sensor_pos(i,2) < -180 || sensor_pos(i,1) > 360,
                        "Unvalid longitude in *sensor_pos*.\n"
                        "Longitudes must be inside [-180,360],\n"
                        "but sensor_pos(",i,",2) is ", sensor_pos(i,2));
  }
  ARTS_USER_ERROR_IF (sensor_los.ncols() != 2,
                      "*sensor_los* must have two columns.");
  ARTS_USER_ERROR_IF (sensor_los.nrows() == 0,
                      "*sensor_los* must have at least one row.");
  for (Index i=0; i<sensor_los.nrows(); i++) {
    ARTS_USER_ERROR_IF (sensor_los(i,0) < 0 || sensor_los(i,0) > 180,
                        "Unvalid zenith angle in *sensor_los*.\n"
                        "Latitudes must be inside [0,180],\n"
                        "but sensor_los(",i,",0) is ", sensor_los(i,0));
    ARTS_USER_ERROR_IF (sensor_los(i,1) < -180 || sensor_los(i,1) > 180,
                        "Unvalid azimuth angle in *sensor_los*.\n"
                        "Latitudes must be inside [-180,180],\n"
                        "but sensor_los(",i,",1) is ", sensor_los(i,1));
  }
  ARTS_USER_ERROR_IF (sensor_los.nrows() != sensor_pos.nrows(),
      "*sensor_los* and *sensor_pos* must have the same number of rows.");

  const Index n = sensor_pos.nrows();
  pos.resize(n, 3);
  los.resize(n, 2);

  Numeric r1 = refellipsoid[0];
  Numeric r2 = r1 * sqrt(1 - refellipsoid[1]*refellipsoid[1]); 
  r1 += altitude;
  r2 += altitude;
  Vector refell_target(2);
  refell_target[0] = r1;
  refell_target[1] = sqrt(1-pow(r2/r1,2.0));
  
  
  const Index atmosphere_dim = 3;
  for (Index i=0; i<n; i++) {
    Numeric l, x, y, z, dx, dy, dz, r;
    r = sensor_pos(i,0) + pos2refell_r(atmosphere_dim,
                                       refellipsoid,
                                       lat_grid,
                                       lon_grid,
                                       sensor_pos(i,joker));
    Numeric za0 = sensor_los(i,0); // For fix below
    poslos2cart(x, y, z, dx, dy, dz,
                r, sensor_pos(i,1), sensor_pos(i,2),
                sensor_los(i,0), sensor_los(i,1));
    line_refellipsoid_intersect(l, refell_target, x, y, z, dx, dy, dz);
    
    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      cart2poslos_plain(pos(i,0),
                        pos(i,1),
                        pos(i,2),
                        los(i,0),
                        los(i,1),
                        x+l*dx,
                        y+l*dy,
                        z+l*dz,
                        dx,
                        dy,
                        dz);
      pos(i,0) -= pos2refell_r(atmosphere_dim,
                               refellipsoid,
                               lat_grid,
                               lon_grid,
                               pos(i,joker));
      // Check that pos matches expected altitude
      //ARTS_ASSERT(abs(pos(i,0)-altitude) < 0.1);

      // za=0 can give NaN for los. A simple fix:
      if (za0 == 0 || za0 == 180) {
        los(i,0) = za0;
        los(i,1) = 0;
      }
    }
  }
}
