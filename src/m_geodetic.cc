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

#include <workspace.h>

#include "check_input.h"
#include "geodetic.h"
#include "matpack_data.h"
#include "surf.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void IntersectionGeometricAltitude(Matrix& pos,
                                   Matrix& los,
                                   const Matrix& sensor_pos,
                                   const Matrix& sensor_los,
                                   const SurfaceField& surface_field,
                                   const Numeric& altitude)
{
  chk_sensor_poslos("sensor_pos", sensor_pos, "sensor_los", sensor_los);
  chk_refellipsoid(surface_field.ellipsoid);

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
                          surface_field.ellipsoid);
        l = intersection_altitude(ecef,
                                  decef,
                                  surface_field.ellipsoid,
                                  altitude);
      }

    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      poslos_at_distance(pos(i,joker), los(i,joker), ecef, decef, surface_field.ellipsoid, l);
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
                                   const SurfaceField& surface_field,
                                   const Numeric& latitude)
{
  chk_sensor_poslos("sensor_pos", sensor_pos, "sensor_los", sensor_los);
  chk_refellipsoid(surface_field.ellipsoid);
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
                      surface_field.ellipsoid);
    l = intersection_latitude(ecef,
                              decef,
                              sensor_pos(i,joker),
                              sensor_los(i,joker), 
                              surface_field.ellipsoid,
                              latitude);
    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      poslos_at_distance(pos(i,joker), los(i,joker), ecef, decef, surface_field.ellipsoid, l);
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
                                    const SurfaceField& surface_field,
                                    const Numeric& longitude)
{
  chk_sensor_poslos("sensor_pos", sensor_pos, "sensor_los", sensor_los);
  chk_refellipsoid(surface_field.ellipsoid);
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
                      surface_field.ellipsoid);
    l = intersection_longitude(ecef,
                               decef,
                               sensor_pos(i,joker),
                               sensor_los(i,joker), 
                               longitude);

    if (l<0) {
      pos(i,joker) = NAN;
      los(i,joker) = NAN;
    } else {
      poslos_at_distance(pos(i,joker), los(i,joker), ecef, decef, surface_field.ellipsoid, l);
      // Check that pos matches expected longitude
      ARTS_ASSERT(abs(pos(i,2)-longitude) < 0.01);
    }
  }
}
