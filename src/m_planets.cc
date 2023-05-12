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

#include "agenda_set.h"
#include "arts.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "auto_md.h"
#include "check_input.h"
#include "matpack_data.h"
#include "surf.h"
#include <cmath>
#include <stdexcept>

inline constexpr Numeric EARTH_RADIUS = Constant::earth_radius;
inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);

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
  // value from MPS, ESA-planetary
  g0 = 23.12;
  // value (1bar level) from
  // http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html g0 = 24.79;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Mars(Numeric &g0) {
  // value from MPS, ESA-planetary
  g0 = 3.690;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Venus(Numeric &g0) {
  // value via MPS, ESA-planetary from Ahrens, 1995
  g0 = 8.870;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void g0Io(Numeric &g0) {
  // value via Wikipedia
  g0 = 1.796;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldEarth(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] = EARTH_RADIUS;
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else if (model == "WGS84") {
    // https://en.wikipedia.org/wiki/World_Geodetic_System#1984_version
    surface_field.ellipsoid[0] = 6378137.0;
    surface_field.ellipsoid[1] = 6356752.314245;
  }

  else
    throw runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldJupiter(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] = 69911e3; // From Ref. 1 (see above)
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else if (model == "Ellipsoid") {
    surface_field.ellipsoid[0] = 71492e3; // From Ref. 1
    surface_field.ellipsoid[1] = 66854e3;
  }

  else
    throw runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldMars(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] = 3389.5e3; // From Ref. 1 (see above)
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else if (model == "Ellipsoid") {
    surface_field.ellipsoid[0] = 3396.19e3; // From Ref. 1
    surface_field.ellipsoid[1] = 3376.20e3;
  }

  else
    throw runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldMoon(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] = 1737.4e3; // From Ref. 1 (see above)
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else if (model == "Ellipsoid") {
    // https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
    surface_field.ellipsoid[0] = 1738.1e3;
    surface_field.ellipsoid[1] = 1736.0e3;
  }

  else
    throw runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldIo(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] =
        1821.6e3; // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else
    throw std::runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldEuropa(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] =
        1560.8e3; // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else
    throw std::runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldGanymede(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] =
        2631e3; // From Wikipedia (and http://ssd.jpl.nasa.gov/?sat_phys_par)
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else
    throw std::runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldVenus(SurfaceField &surface_field, const String &model) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;

  if (model == "Sphere") {
    surface_field.ellipsoid[0] = 6051.8e3; // From Ref. 1 (see above)
    surface_field.ellipsoid[1] = surface_field.ellipsoid[0];
  }

  else
    throw runtime_error("Unknown selection for input argument *model*.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_fieldInit(SurfaceField &surface_field, const Numeric &r_equatorial,
                       const Numeric &r_polar) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;
  
  surface_field.ellipsoid[0] = r_equatorial;
  surface_field.ellipsoid[1] = r_polar;
  chk_refellipsoid(surface_field.ellipsoid);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void PlanetSet(Workspace &ws, Agenda &g0_agenda, SurfaceField &surface_field,
               Numeric &molarmass_dry_air, Numeric &planet_rotation_period,
               const String &option) {
  surface_field = {};
  surface_field[Surf::Key::h] = 0.0;
  molarmass_dry_air = 0.0;
  planet_rotation_period = 0.0;

  using enum Options::planetDefaultOptions;
  switch (Options::toplanetDefaultOptionsOrThrow(option)) {
  case Earth:
    surface_fieldEarth(surface_field, "Sphere");
    molarmass_dry_air = 28.966;
    planet_rotation_period = 86164.1;
    break;
  case Io:
    surface_fieldIo(surface_field, "Sphere");
    molarmass_dry_air = 63.110068828000003;
    planet_rotation_period = 152853;
    break;
  case Jupiter:
    surface_fieldJupiter(surface_field, "Sphere");
    molarmass_dry_air = 2.22;
    planet_rotation_period = 35730;
    break;
  case Mars:
    surface_fieldMars(surface_field, "Sphere");
    molarmass_dry_air = 43.34;
    planet_rotation_period = 88643;
    break;
  case Venus:
    surface_fieldVenus(surface_field, "Sphere");
    molarmass_dry_air = 43.45;
    planet_rotation_period = -2.0997e7;
    break;
  case FINAL:
    break;
  }

  g0_agenda = AgendaManip::get_g0_agenda(ws, option);
}
