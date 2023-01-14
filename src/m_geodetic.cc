/* Copyright (C) 2023
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

/**
    @file    m_geodetic.cc
    @author  Patrick Eriksson <patrick.eriksson@chalmers.se>
    @date    2023-01-14
 
    @brief   Workspace methods performing geodetic calculations of
             relatively general character.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "matpack_data.h"
#include "messages.h"
#include "variousZZZ.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void IntersectionGeometricAltitude(Matrix& pos,
                                   Matrix& los,
                                   const Matrix& sensor_pos,
                                   const Matrix& sensor_los,
                                   const Vector& refellipsoid,
                                   const Numeric& altitude,
                                   const Verbosity&)
{
  chk_sensor_poslos("sensor_pos", sensor_pos, "sensor_los", sensor_los);
  chk_refellipsoidZZZ(refellipsoid);

  const Index n = sensor_pos.nrows();
  pos.resize(n, 3);
  los.resize(n, 2);

  for (Index i=0; i<n; i++) {
    Numeric l;
    Vector ecef(3), decef(3);

    // Intersection not possible of position below altitude or view upwards
    // (the first part allowed by *intersection_altitude*)
    if (sensor_pos(i,0) < altitude || sensor_pos(i,0)<=90)
      { l = -1; }
    // Calculate distance to intersection 
    else
      {   
        geodetic_los2ecef(ecef,
                          decef,
                          sensor_pos(i,joker),
                          sensor_los(i,joker), 
                          refellipsoid);
        l = intersection_altitude(ecef,
                                  decef,
                                  refellipsoid,
                                  altitude);
      }

    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      poslos_at_distance(pos(i,joker), los(i,joker), ecef, decef, refellipsoid, l);
      // Check that pos matches expected altitude
      ARTS_ASSERT(abs(pos(i,0)-altitude) < 0.1);
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IntersectionGeometricLatitude(Matrix& pos,
                                   Matrix& los,
                                   const Matrix& sensor_pos,
                                   const Matrix& sensor_los,
                                   const Vector& refellipsoid,
                                   const Numeric& latitude,
                                   const Verbosity&)
{
  chk_sensor_poslos("sensor_pos", sensor_pos, "sensor_los", sensor_los);
  chk_refellipsoidZZZ(refellipsoid);
  chk_if_in_range("latitude", latitude, -90, 90);

  const Index n = sensor_pos.nrows();
  pos.resize(n, 3);
  los.resize(n, 2);

  for (Index i=0; i<n; i++) {
    Numeric l;
    Vector ecef(3), decef(3);

    geodetic_los2ecef(ecef,
                      decef,
                      sensor_pos(i,joker),
                      sensor_los(i,joker), 
                      refellipsoid);
    l = intersection_latitude(ecef,
                              decef,
                              sensor_pos(i,joker),
                              sensor_los(i,joker), 
                              refellipsoid,
                              latitude);
    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      poslos_at_distance(pos(i,joker), los(i,joker), ecef, decef, refellipsoid, l);
      // Check that pos matches expected longitude
      ARTS_ASSERT(abs(pos(i,1)-latitude) < 0.01);
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IntersectionGeometricLongitude(Matrix& pos,
                                    Matrix& los,
                                    const Matrix& sensor_pos,
                                    const Matrix& sensor_los,
                                    const Vector& refellipsoid,
                                    const Numeric& longitude,
                                    const Verbosity&)
{
  chk_sensor_poslos("sensor_pos", sensor_pos, "sensor_los", sensor_los);
  chk_refellipsoidZZZ(refellipsoid);
  chk_if_in_range("longitude", longitude, -180, 360);

  const Index n = sensor_pos.nrows();
  pos.resize(n, 3);
  los.resize(n, 2);

  for (Index i=0; i<n; i++) {
    Numeric l;
    Vector ecef(3), decef(3);

    geodetic_los2ecef(ecef,
                      decef,
                      sensor_pos(i,joker),
                      sensor_los(i,joker), 
                      refellipsoid);
    l = intersection_longitude(ecef,
                               decef,
                               sensor_pos(i,joker),
                               sensor_los(i,joker), 
                               longitude);

    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      poslos_at_distance(pos(i,joker), los(i,joker), ecef, decef, refellipsoid, l);
      // Check that pos matches expected longitude
      ARTS_ASSERT(abs(pos(i,2)-longitude) < 0.01);
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IntersectionGeometricSurface(Matrix& pos,
                                  Matrix& los,
                                  const Matrix& sensor_pos,
                                  const Matrix& sensor_los,
                                  const Vector& refellipsoid,
                                  const GriddedField2& surface_elevation,
                                  const Numeric& surface_search_accuracy,
                                  const Index& surface_search_safe,
                                  const Verbosity&) {
  chk_sensor_poslos("sensor_pos", sensor_pos, "sensor_los", sensor_los);
  chk_refellipsoid(refellipsoid);
  chk_surface_elevation(surface_elevation);
  chk_if_positive("surface_search_accuracy", surface_search_accuracy);
  chk_if_bool("surface_search_safe", surface_search_safe);

  const Index n = sensor_pos.nrows();
  pos.resize(n, 3);
  los.resize(n, 2);

  for (Index i=0; i<n; i++) {
    Vector ecef(3), decef(3);
    geodetic_los2ecef(ecef,
                      decef,
                      sensor_pos(i,joker),
                      sensor_los(i,joker), 
                      refellipsoid);
    const Numeric l = find_crossing_with_surface_z(sensor_pos(i,joker),
                                                   sensor_los(i,joker),
                                                   ecef,
                                                   decef,
                                                   refellipsoid,
                                                   surface_elevation,
                                                   surface_search_accuracy,
                                                   surface_search_safe);
    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      poslos_at_distance(pos(i,joker), los(i,joker), ecef, decef, refellipsoid, l);
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void TestBasicGeodeticAccuracy(Vector& rte_pos,
                               Numeric& max_dl,
                               Vector& max_dpos,
                               Vector& max_dlos,
                               const Vector& refellipsoid,
                               const Index& ntests,
                               const Numeric& max_allowed_dl,
                               const Verbosity& verbosity) {
  // Seed random generator
  Index seed;
  MCSetSeedFromTime(seed, verbosity);
  Rng rng;
  rng.seed(seed, verbosity);

  // Init GOUTs
  max_dl = 0;
  max_dpos.resize(3);
  max_dpos = 0;
  max_dlos.resize(2);
  max_dlos = 0;

  for (Index i=0; i<ntests; i++) {
    // Create a random pos and los
    Vector pos(3), los(2);
    pos[0] = 100e3*rng.draw() - 100.0;
    pos[1] = 180*rng.draw() - 90.0;
    pos[2] = 360*rng.draw() - 180.0;
    los[0] = 180*rng.draw();    
    los[1] = 360*rng.draw() - 180.0;

    // Convert to ECEF and back, and then back to ECEF to calculate distance
    Vector ecef(3), decef(3), pos_new(3), los_new(2), ecef_new(3), decef_new(3);
    geodetic_los2ecef(ecef, decef, pos, los, refellipsoid);
    ecef2geodetic_los(pos_new, los_new,ecef, decef, refellipsoid);
    geodetic_los2ecef(ecef_new, decef_new, pos_new, los_new, refellipsoid);
    const Numeric dl = ecef_distance(ecef, ecef_new);

    if (dl > max_allowed_dl) {
      ARTS_USER_ERROR("Maximum allowed error exceeded!\n"
                      "  Set distance error limit: ", max_allowed_dl, "\n"
                      "         Calculated  error: ", max_allowed_dl, "\n");
    }

    if (dl > max_dl) {
      max_dl = dl;
      rte_pos = pos;
    }

    Numeric d;
    for (Index j=0; j<3; j++) {
      d = abs(pos_new[j] - pos[j]);
      if (d > max_dpos[j])
        max_dpos[j] = d;
      if (j<2) {
        d = abs(los_new[j] - los[j]);
        if (d > max_dlos[j])
          max_dlos[j] = d;
      }
    }
  }
}
