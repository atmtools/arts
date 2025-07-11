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

#include <enumsEarthEllipsoid.h>
#include <enumsEuropaEllipsoid.h>
#include <enumsGanymedeEllipsoid.h>
#include <enumsIoEllipsoid.h>
#include <enumsJupiterEllipsoid.h>
#include <enumsMarsEllipsoid.h>
#include <enumsMoonEllipsoid.h>
#include <enumsPlanetOrMoonType.h>
#include <enumsVenusEllipsoid.h>
#include <planet_data.h>
#include <workspace.h>

#include <cmath>
#include <stdexcept>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "debug.h"
#include "operators.h"
#include "rtepack.h"
#include "surf.h"

inline constexpr Numeric EARTH_RADIUS = Constant::earth_radius;
inline constexpr Numeric DEG2RAD      = Conversion::deg2rad(1);

// Ref. 1:
// Seidelmann, P. Kenneth; Archinal, B. A.; A'hearn, M. F. et al (2007).
// "Report of the IAU/IAG Working Group on cartographic coordinates and
// rotational elements: 2006". Celestial Mechanics and Dynamical Astronomy 98
// (3): 155–180. Bibcode 2007CeMDA..98..155S. doi:10.1007/s10569-007-9072-y

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Earth(Numeric &g0, const Numeric &lat) {
  ARTS_TIME_REPORT

  // "Small g" at altitude=0, g0:
  // Expression for g0 taken from Wikipedia page "Gravity of Earth", that
  // is stated to be: International Gravity Formula 1967, the 1967 Geodetic
  // Reference System Formula, Helmert's equation or Clairault's formula.

  const Numeric x = DEG2RAD * fabs(lat);

  g0 = 9.780327 *
       (1 + 5.3024e-3 * pow(sin(x), 2.0) + 5.8e-6 * pow(sin(2 * x), 2.0));

  // Move to apparent gravity, i.e. include effect of the centrifugal force.
  // See: A first course in Atmospheric Thermodynamics by G. Petty (page 89) As
  // well as https://glossary.ametsoc.org/wiki/Apparent_gravity 0.033895 =
  // (7.29e-5)^2 * 6378e3
  g0 -= 0.033895 * pow(cos(x), 2.0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Jupiter(Numeric &g0) {
  ARTS_TIME_REPORT

  // value from MPS, ESA-planetary
  g0 = 23.12;
  // value (1bar level) from
  // http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html g0 = 24.79;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Mars(Numeric &g0) {
  ARTS_TIME_REPORT

  // value from MPS, ESA-planetary
  g0 = 3.690;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Venus(Numeric &g0) {
  ARTS_TIME_REPORT

  // value via MPS, ESA-planetary from Ahrens, 1995
  g0 = 8.870;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Io(Numeric &g0) {
  ARTS_TIME_REPORT

  // value via Wikipedia
  g0 = 1.796;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldEarth(SurfaceField &surface_field,
                        const String &model,
                        const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<EarthEllipsoid>(model)) {
    case EarthEllipsoid::WGS84:
      // https://en.wikipedia.org/wiki/World_Geodetic_System#1984_version
      surface_field.ellipsoid[0] = Body::Earth::a;
      surface_field.ellipsoid[1] = Body::Earth::b;
      break;
    case EarthEllipsoid::Sphere:
      surface_field.ellipsoid[0] = EARTH_RADIUS;
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldJupiter(SurfaceField &surface_field,
                          const String &model,
                          const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<JupiterEllipsoid>(model)) {
    case JupiterEllipsoid::Sphere:
      surface_field.ellipsoid[0] = 69911e3;  // From Ref. 1 (see above)
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
    case JupiterEllipsoid::Ellipsoid:
      surface_field.ellipsoid[0] = Body::Jupiter::a;  // From Ref. 1
      surface_field.ellipsoid[1] = Body::Jupiter::b;
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldMars(SurfaceField &surface_field,
                       const String &model,
                       const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<MarsEllipsoid>(model)) {
    case MarsEllipsoid::Sphere:
      surface_field.ellipsoid[0] = 3389.5e3;  // From Ref. 1 (see above)
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
    case MarsEllipsoid::Ellipsoid:
      surface_field.ellipsoid[0] = Body::Mars::a;  // From Ref. 1
      surface_field.ellipsoid[1] = Body::Mars::b;
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldMoon(SurfaceField &surface_field,
                       const String &model,
                       const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<MoonEllipsoid>(model)) {
    case MoonEllipsoid::Sphere:
      surface_field.ellipsoid[0] = 1737.4e3;  // From Ref. 1 (see above)
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
    case MoonEllipsoid::Ellipsoid:
      // https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
      surface_field.ellipsoid[0] = Body::Moon::a;
      surface_field.ellipsoid[1] = Body::Moon::b;
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldIo(SurfaceField &surface_field,
                     const String &model,
                     const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<IoEllipsoid>(model)) {
    case IoEllipsoid::Sphere:
      surface_field.ellipsoid[0] =
          1821.6e3;  // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldEuropa(SurfaceField &surface_field,
                         const String &model,
                         const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<EuropaEllipsoid>(model)) {
    case EuropaEllipsoid::Sphere:
      surface_field.ellipsoid[0] =
          1560.8e3;  // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldGanymede(SurfaceField &surface_field,
                           const String &model,
                           const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<GanymedeEllipsoid>(model)) {
    case GanymedeEllipsoid::Sphere:
      surface_field.ellipsoid[0] =
          2631e3;  // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldVenus(SurfaceField &surface_field,
                        const String &model,
                        const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  switch (to<VenusEllipsoid>(model)) {
    case VenusEllipsoid::Sphere:
      surface_field.ellipsoid[0] = 6051.8e3;  // From Ref. 1 (see above)
      surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
      break;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldInit(SurfaceField &surface_field,
                       const Numeric &r_equatorial,
                       const Numeric &r_polar,
                       const Numeric &surface_elevation) {
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

  surface_field                = {};
  surface_field[SurfaceKey::h] = surface_elevation;

  surface_field.ellipsoid[0] = r_equatorial;
  surface_field.ellipsoid[1] = r_polar;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldPlanet(SurfaceField &surface_field,
                         const String &option,
                         const Numeric &surface_elevation) {
  ARTS_TIME_REPORT

  using enum PlanetOrMoonType;
  switch (to<PlanetOrMoonType>(option)) {
    case Earth:
      surface_fieldEarth(surface_field, "WGS84", surface_elevation);
      // molarmass_dry_air = 28.966;
      // planet_rotation_period = 86164.1;
      break;
    case Io:
      surface_fieldIo(surface_field, "Sphere", surface_elevation);
      // molarmass_dry_air = 63.110068828000003;
      // planet_rotation_period = 152853;
      break;
    case Jupiter:
      surface_fieldJupiter(surface_field, "Sphere", surface_elevation);
      // molarmass_dry_air = 2.22;
      // planet_rotation_period = 35730;
      break;
    case Mars:
      surface_fieldMars(surface_field, "Sphere", surface_elevation);
      // molarmass_dry_air = 43.34;
      // planet_rotation_period = 88643;
      break;
    case Venus:
      surface_fieldVenus(surface_field, "Sphere", surface_elevation);
      // molarmass_dry_air = 43.45;
      // planet_rotation_period = -2.0997e7;
      break;
  }
}

void gravity_operatorCentralMass(NumericTernaryOperator &gravity_operator,
                                 const SurfaceField &surface_field,
                                 const Numeric &mass) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(surface_field.ellipsoid[0] <= 0,
                     "Ellipsoid has bad semi-major axis {:B,}",
                     surface_field.ellipsoid)
  ARTS_USER_ERROR_IF(
      surface_field.ellipsoid[1] <= 0 or
          surface_field.ellipsoid[1] > surface_field.ellipsoid[0],
      "Ellipsoid has bad semi-minor axis {:B,}",
      surface_field.ellipsoid)

  gravity_operator = NumericTernaryOperator{
      Gravity{.GM = Constant::G * mass,
              .a  = surface_field.ellipsoid[0],
              .e  = std::sqrt(1 - Math::pow2(surface_field.ellipsoid[1] /
                                            surface_field.ellipsoid[0]))}};
}
