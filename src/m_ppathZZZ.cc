/* Copyright (C) 2021 Patrick Eriksson <Patrick.Eriksson@chalmers.se>
                            
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
 * @file   m_ppath.cc
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2021-08-11 
 *
 * @brief  Workspace functions releated to calculation of propagation paths.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "geodeticZZZ.h"
#include "math_funcs.h"
#include "ppathZZZ.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathCheckEndPoint(const Ppath& ppath,
                        const Index& background,
                        const Index& np,
                        const Numeric& altitude,
                        const Numeric& daltitude,
                        const Numeric& latitude,
                        const Numeric& dlatitude,
                        const Numeric& longitude,
                        const Numeric& dlongitude,
                        const Numeric& zenith_angle,
                        const Numeric& dzenith_angle,
                        const Numeric& azimuth_angle,
                        const Numeric& dazimuth_angle,
                        const Verbosity&) {
  // pos and los to check
  ConstVectorView pos = ppath.end_pos, los = ppath.end_los;

  if (background >= 0 && ppath.backgroundZZZ != background) {
    ARTS_USER_ERROR(
        "Radiative background not as expected!\n"
        "  background in ppath: ",
        ppath.backgroundZZZ,
        "\n  background expected: ",
        background);
  }
  if (np >= 0 && ppath.np != np) {
    ARTS_USER_ERROR(
        "Number of ppath points not as expected!\n"
        "  number in ppath: ",
        ppath.np,
        "\n  number expected: ",
        np);
  }
  if (daltitude >= 0 && abs(pos[0] - altitude) > daltitude) {
    ARTS_USER_ERROR(
        "Start altitude not as expected!\n"
        "  altitude in ppath: ",
        pos[0],
        "\n  altitude expected: ",
        altitude,
        "\n         difference: ",
        abs(pos[0] - altitude),
        "\n      set tolarance: ",
        daltitude);
  }
  if (dlatitude >= 0 && abs(pos[1] - latitude) > dlatitude) {
    ARTS_USER_ERROR(
        "Start latitude not as expected!\n"
        "  latitude in ppath: ",
        pos[1],
        "\n  latitude expected: ",
        latitude,
        "\n         difference: ",
        abs(pos[1] - latitude),
        "\n      set tolarance: ",
        dlatitude);
  }
  if (dlongitude >= 0 && abs(pos[2] - longitude) > dlongitude) {
    ARTS_USER_ERROR(
        "Start longitude not as expected!\n"
        "  longitude in ppath: ",
        pos[2],
        "\n  longitude expected: ",
        longitude,
        "\n          difference: ",
        abs(pos[2] - longitude),
        "\n       set tolarance: ",
        dlongitude);
  }
  if (dzenith_angle >= 0 && abs(los[0] - zenith_angle) > dzenith_angle) {
    ARTS_USER_ERROR(
        "Start zenith angle not as expected!\n"
        "  zenith angle in ppath: ",
        los[0],
        "\n  zenith angle expected: ",
        zenith_angle,
        "\n             difference: ",
        abs(los[0] - zenith_angle),
        "\n          set tolarance: ",
        dzenith_angle);
  }
  if (dazimuth_angle >= 0 && abs(los[1] - azimuth_angle) > dazimuth_angle) {
    ARTS_USER_ERROR(
        "Start azimuth angle not as expected!\n"
        "  azimuth angle in ppath: ",
        los[1],
        "\n  azimuth angle expected: ",
        azimuth_angle,
        "\n              difference: ",
        abs(los[1] - azimuth_angle),
        "\n           set tolarance: ",
        dazimuth_angle);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ppathPassive(Ppath& ppath,
                  const Vector& rte_pos,
                  const Vector& rte_los,
                  const Index& atmosphere_dim,
                  const Vector& refellipsoid,
                  const Vector& z_grid,
                  const Vector& lat_grid,
                  const Vector& lon_grid,
                  const Index& cloudbox_on,
                  const ArrayOfIndex& cloudbox_limits,
                  const GriddedField2& surface_elevation,
                  const Index& refraction_do,
                  const Index& add_grid_crossings,
                  const Numeric& l_step_max,
                  const Numeric& l_total_max,
                  const Numeric& l_raytrace_geom,
                  const Numeric& l_raytrace_refr,
                  const Numeric& l_accuracy,
                  const Index& safe_surface_search,
                  const Index& do_not_calc_gps,
                  const Verbosity&) {
  // Check lengths always used
  if (l_step_max <= 0) {
    ARTS_USER_ERROR("GIN l_step_max must be > 0.\n");
  }
  if (l_accuracy <= 0) {
    ARTS_USER_ERROR("GIN l_accuracy must be > 0.\n");
  }

  // Step 1, with refraction
  if (refraction_do) {
    if (l_raytrace_refr <= 0 || l_raytrace_refr > l_step_max / 4) {
      ARTS_USER_ERROR("GIN l_raytrace_geom must be > 0 and < l_step_max/4.\n");
    }
    ARTS_USER_ERROR("Refraction not yet implemented\n", l_raytrace_refr);

    // Step 1, geometrical
  } else {
    if (l_raytrace_geom <= 0) {
      ARTS_USER_ERROR("GIN l_raytrace_geom must be > 0.\n");
    }
    ppath_geom_const_lstep(ppath,
                           rte_pos,
                           rte_los,
                           atmosphere_dim,
                           refellipsoid,
                           z_grid,
                           lat_grid,
                           lon_grid,
                           cloudbox_on,
                           cloudbox_limits,
                           surface_elevation,
                           add_grid_crossings ? l_raytrace_geom : l_step_max,
                           l_total_max,
                           l_accuracy,
                           safe_surface_search,
                           do_not_calc_gps && !add_grid_crossings);
  }

  // Add grid crossings?
  if (add_grid_crossings) {
    ppath_grid_crossings(ppath,
                         atmosphere_dim,
                         refellipsoid,
                         z_grid,
                         lat_grid,
                         lon_grid,
                         l_step_max,
                         l_accuracy,
                         do_not_calc_gps);
  }
}
