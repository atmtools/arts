/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   m_planets.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2012-03-16

  \brief  Planet specific workspace methods.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>
#include <enumsEarthEllipsoid.h>
#include <enumsEuropaEllipsoid.h>
#include <enumsGanymedeEllipsoid.h>
#include <enumsIoEllipsoid.h>
#include <enumsJupiterEllipsoid.h>
#include <enumsMarsEllipsoid.h>
#include <enumsMoonEllipsoid.h>
#include <enumsPlanetOrMoonType.h>
#include <enumsVenusEllipsoid.h>
#include <operators.h>
#include <planet_data.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <cmath>

inline constexpr Numeric EARTH_RADIUS = Constant::earth_radius;

// Ref. 1:
// Seidelmann, P. Kenneth; Archinal, B. A.; A'hearn, M. F. et al (2007).
// "Report of the IAU/IAG Working Group on cartographic coordinates and
// rotational elements: 2006". Celestial Mechanics and Dynamical Astronomy 98
// (3): 155–180. Bibcode 2007CeMDA..98..155S. doi:10.1007/s10569-007-9072-y

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldEarth(SurfaceField &surf_field,
                     const String &model,
                     const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<EarthEllipsoid>(model)) {
    case EarthEllipsoid::WGS84:
      // https://en.wikipedia.org/wiki/World_Geodetic_System#1984_version
      surf_field.ellipsoid[0] = Body::Earth::a;
      surf_field.ellipsoid[1] = Body::Earth::b;
      break;
    case EarthEllipsoid::Sphere:
      surf_field.ellipsoid[0] = EARTH_RADIUS;
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldJupiter(SurfaceField &surf_field,
                       const String &model,
                       const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<JupiterEllipsoid>(model)) {
    case JupiterEllipsoid::Sphere:
      surf_field.ellipsoid[0] = 69911e3;  // From Ref. 1 (see above)
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
    case JupiterEllipsoid::Ellipsoid:
      surf_field.ellipsoid[0] = Body::Jupiter::a;  // From Ref. 1
      surf_field.ellipsoid[1] = Body::Jupiter::b;
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldMars(SurfaceField &surf_field,
                    const String &model,
                    const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<MarsEllipsoid>(model)) {
    case MarsEllipsoid::Sphere:
      surf_field.ellipsoid[0] = 3389.5e3;  // From Ref. 1 (see above)
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
    case MarsEllipsoid::Ellipsoid:
      surf_field.ellipsoid[0] = Body::Mars::a;  // From Ref. 1
      surf_field.ellipsoid[1] = Body::Mars::b;
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldMoon(SurfaceField &surf_field,
                    const String &model,
                    const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<MoonEllipsoid>(model)) {
    case MoonEllipsoid::Sphere:
      surf_field.ellipsoid[0] = 1737.4e3;  // From Ref. 1 (see above)
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
    case MoonEllipsoid::Ellipsoid:
      // https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
      surf_field.ellipsoid[0] = Body::Moon::a;
      surf_field.ellipsoid[1] = Body::Moon::b;
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldIo(SurfaceField &surf_field,
                  const String &model,
                  const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<IoEllipsoid>(model)) {
    case IoEllipsoid::Sphere:
      surf_field.ellipsoid[0] =
          1821.6e3;  // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldEuropa(SurfaceField &surf_field,
                      const String &model,
                      const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<EuropaEllipsoid>(model)) {
    case EuropaEllipsoid::Sphere:
      surf_field.ellipsoid[0] =
          1560.8e3;  // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldGanymede(SurfaceField &surf_field,
                        const String &model,
                        const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<GanymedeEllipsoid>(model)) {
    case GanymedeEllipsoid::Sphere:
      surf_field.ellipsoid[0] =
          2631e3;  // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldVenus(SurfaceField &surf_field,
                     const String &model,
                     const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  switch (to<VenusEllipsoid>(model)) {
    case VenusEllipsoid::Sphere:
      surf_field.ellipsoid[0] = 6051.8e3;  // From Ref. 1 (see above)
      surf_field.ellipsoid[1] = surf_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldInit(SurfaceField &surf_field,
                    const Numeric &r_equatorial,
                    const Numeric &r_polar,
                    const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(r_equatorial <= 0 || r_polar <= 0,
                     R"(Ellipsoid with a: {}, b: {} is invalid.

We do not support negative "a" or "b".

A planetary ellipsoid is defined by x^2/a^2 + y^2/a^2 + z^2/b^2 = 1 in ARTS.
)",
                     r_equatorial,
                     r_polar);
  ARTS_USER_ERROR_IF(std::abs(r_polar / r_equatorial - 1) > 0.5,
                     R"(Ellipsoid with a: {}, b: {} is invalid.

The ellipsoid is too flat.

A planetary ellipsoid is defined by x^2/a^2 + y^2/a^2 + z^2/b^2 = 1 in ARTS.
)",
                     r_equatorial,
                     r_polar);

  surf_field                = {};
  surf_field[SurfaceKey::h] = surf_elevation;

  surf_field.ellipsoid[0] = r_equatorial;
  surf_field.ellipsoid[1] = r_polar;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surf_fieldPlanet(SurfaceField &surf_field,
                      const String &option,
                      const Numeric &surf_elevation) {
  ARTS_TIME_REPORT

  using enum PlanetOrMoonType;
  switch (to<PlanetOrMoonType>(option)) {
    case Earth:
      surf_fieldEarth(surf_field, "WGS84", surf_elevation);
      // molarmass_dry_air = 28.966;
      // planet_rotation_period = 86164.1;
      break;
    case Io:
      surf_fieldIo(surf_field, "Sphere", surf_elevation);
      // molarmass_dry_air = 63.110068828000003;
      // planet_rotation_period = 152853;
      break;
    case Jupiter:
      surf_fieldJupiter(surf_field, "Sphere", surf_elevation);
      // molarmass_dry_air = 2.22;
      // planet_rotation_period = 35730;
      break;
    case Mars:
      surf_fieldMars(surf_field, "Sphere", surf_elevation);
      // molarmass_dry_air = 43.34;
      // planet_rotation_period = 88643;
      break;
    case Venus:
      surf_fieldVenus(surf_field, "Sphere", surf_elevation);
      // molarmass_dry_air = 43.45;
      // planet_rotation_period = -2.0997e7;
      break;
  }
}

void gravity_operatorCentralMass(NumericTernaryOperator &gravity_operator,
                                 const SurfaceField &surf_field,
                                 const Numeric &mass) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(surf_field.ellipsoid[0] <= 0,
                     "Ellipsoid has bad semi-major axis {:B,}",
                     surf_field.ellipsoid)
  ARTS_USER_ERROR_IF(surf_field.ellipsoid[1] <= 0 or
                         surf_field.ellipsoid[1] > surf_field.ellipsoid[0],
                     "Ellipsoid has bad semi-minor axis {:B,}",
                     surf_field.ellipsoid)

  gravity_operator = NumericTernaryOperator{EllipsoidGravity{
      .GM = Constant::G * mass,
      .a  = surf_field.ellipsoid[0],
      .e  = std::sqrt(
          1 - Math::pow2(surf_field.ellipsoid[1] / surf_field.ellipsoid[0]))}};
}
