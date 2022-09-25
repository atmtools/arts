/* Copyright (C) 2021 Patrick Eriksson <patrick.eriksson@chalmers.se>

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
 * @file   ppath_funcs.cc
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2021-07-28
 * 
 * @brief  Higher level functions related to calculation of propagation paths.
 *
 * The term propagation path is here shortened to ppath.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "ppath_funcs.h"

#include "geodeticZZZ.h"
#include "math_funcs.h"
#include "ppath.h"
#include "ppathZZZ.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);


/*===========================================================================
  === The functions, in alphabetical order
  ===========================================================================*/

Numeric find_crossing_with_surface_z(const Vector rte_pos,
                                     const Vector rte_los,
                                     const Vector ecef,
                                     const Vector decef,
                                     const Index& atmosphere_dim,
                                     const Vector& refellipsoid,
                                     const GriddedField2& surface_elevation,
                                     const Numeric& l_accuracy,
                                     const Index& safe_search) {
  // Find min and max surface altitude
  Numeric z_min, z_max;
  if (atmosphere_dim == 1) {
    z_min = surface_elevation.data(0, 0);
    z_max = surface_elevation.data(0, 0);
  } else {
    z_min = min(surface_elevation.data);
    z_max = max(surface_elevation.data);
  }

  // Catch upward looking cases that can not have a surface intersection
  if (rte_pos[0] >= z_max && rte_los[0] <= 90) {
    return -1;
  }

  // Check that observation position is above ground
  if (rte_pos[0] < z_max) {
    Numeric z_surf =
        surface_z_at_pos(rte_pos, atmosphere_dim, surface_elevation);
    if (rte_pos[0] < z_surf)
      ARTS_USER_ERROR(
          "The sensor is below the surface. Not allowed!\n"
          "The sensor altitude is at ",
          rte_pos[0],
          " m\n"
          "The surface altitude is ",
          z_surf,
          " m\n"
          "The position is (lat,lon): (",
          rte_pos[1],
          ",",
          rte_pos[2],
          ")");
  }

  // Constant surface altitude (in comparison to *l_accuracy*)
  if (atmosphere_dim == 1 || z_max - z_min < l_accuracy / 100) {
    return intersection_altitude(ecef, decef, refellipsoid, z_min);

    // The general case
  } else {
    // Find max distance that is guaranteed above or at surface
    // This is the minimum distance to the surface
    // If below z_max, this distance is 0. Otherwise given by z_max
    Numeric l_min;
    if (rte_pos[0] <= z_max)
      l_min = 0;
    else {
      l_min = intersection_altitude(ecef, decef, refellipsoid, z_max);
      // No intersection if not even z_max is reached
      if (l_min < 0) return -1;
    }
    // Find max distance for search.
    // If below z_max and upward, given by z_max
    // Otherwise in general given by z_min. If z_min not reached, the distance
    // is instead given by tangent point
    Numeric l_max;
    bool l_max_could_be_above_surface = false;
    if (rte_pos[0] <= z_max && rte_los[0] <= 90) {
      l_max = intersection_altitude(ecef, decef, refellipsoid, z_max);
      l_max_could_be_above_surface = true;
    } else {
      l_max = intersection_altitude(ecef, decef, refellipsoid, z_min);
    }
    if (l_max < 0) {
      Vector ecef_tan(3);
      approx_geometrical_tangent_point(ecef_tan, ecef, decef, refellipsoid);
      l_max = ecef_distance(ecef, ecef_tan);
      // To not miss intersections just after the tangent point, we add a
      // a distance that depends om planet radius (for Earth 111 km).
      l_max += refellipsoid[0] * sin(DEG2RAD);
      l_max_could_be_above_surface = true;
    }

    // Safe but slow approach
    // ----------------------
    if (safe_search) {
      Numeric l_test = l_min - l_accuracy / 2;
      bool above_surface = true;  // if true l_test is 0
      while (above_surface && l_test < l_max) {
        l_test += l_accuracy;
        Vector ecef_test, pos(3);
        ecef_at_distance(ecef_test, ecef, decef, l_test);
        ecef2geodetic(pos, ecef_test, refellipsoid);
        Numeric z_surf =
            surface_z_at_pos(pos, atmosphere_dim, surface_elevation);
        if (pos[0] < z_surf) above_surface = false;
      }
      if (above_surface) {
        return -1;
      } else {
        return l_test - l_accuracy / 2;
      }

      // Bisection search
      // ----------------------
    } else {
      // If l_max matches a point above the surface, we have no intersection
      // according to this search algorithm. And the search fails. So we need
      // to check that point if status unclear
      if (l_max_could_be_above_surface) {
        Vector ecef_test, pos(3);
        ecef_at_distance(ecef_test, ecef, decef, l_max);
        ecef2geodetic(pos, ecef_test, refellipsoid);
        Numeric z_surf =
            surface_z_at_pos(pos, atmosphere_dim, surface_elevation);
        if (pos[0] > z_surf) return -1;
      }
      // Start bisection
      while (l_max - l_min > 2 * l_accuracy) {
        const Numeric l_test = (l_min + l_max) / 2;
        Vector ecef_test, pos(3);
        ecef_at_distance(ecef_test, ecef, decef, l_test);
        ecef2geodetic(pos, ecef_test, refellipsoid);
        Numeric z_surf =
            surface_z_at_pos(pos, atmosphere_dim, surface_elevation);
        if (pos[0] >= z_surf)
          l_min = l_test;
        else
          l_max = l_test;
      }
      return (l_min + l_max) / 2;
    }
  }
}

enum PpathBackground ppath_init_calc(Numeric& l2toa,
                                     VectorView pos_toa,
                                     const Vector rte_pos,
                                     const Vector rte_los,
                                     ConstVectorView ecef,
                                     ConstVectorView decef,
                                     const Index& atmosphere_dim,
                                     const Vector& refellipsoid,
                                     const Vector& z_grid,
                                     const Vector& lat_grid,
                                     const Vector& lon_grid,
                                     const Index& cloudbox_on,
                                     const ArrayOfIndex& cloudbox_limits) {
  // Default values
  l2toa = -1;
  pos_toa = 0;

  // Inside or outside of atmosphere?
  const Numeric z_toa = last(z_grid);
  const Index outside = is_pos_outside_atmosphere(
      rte_pos, atmosphere_dim, z_toa, lat_grid, lon_grid);

  // -------------------------------------------------------------------------
  // Outside of atmosphere
  // -------------------------------------------------------------------------
  if (outside) {
    // Outside but sensor below TOA
    if (outside == 2) {
      ARTS_USER_ERROR(
          "The sensor is below the top-of-atmosphere (TOA) altitude,\n"
          "but outside in the latitude range. This is not allowed.\n"
          " Sensor altitude: ",
          rte_pos[0],
          " m\n"
          "    TOA altitude: ",
          z_toa,
          "m\n"
          " Sensor latitude: ",
          rte_pos[1],
          "\n"
          "   Latitude grid: ",
          lat_grid[0],
          "-",
          last(lat_grid));
    }
    if (outside == 3) {
      ARTS_USER_ERROR(
          "The sensor is below the top-of-atmosphere (TOA) altitude,\n"
          "but outside in the longitude range. This is not allowed.\n"
          "  Sensor altitude: ",
          rte_pos[0],
          " m\n"
          "     TOA altitude: ",
          z_toa,
          "m\n"
          " Sensor longitude: ",
          rte_pos[2],
          "\n"
          "   Longitude grid: ",
          lon_grid[0],
          "-",
          last(lon_grid));
    }

    // If looking up, space is background
    if (rte_los[0] >= 90) {
      return PPATH_BACKGROUND_SPACE;
    }

    // Otherwise, do we reach TOA?
    l2toa = intersection_altitude(ecef, decef, refellipsoid, z_toa);
    if (l2toa < 00) {
      return PPATH_BACKGROUND_SPACE;
    }

    // If TOA reached, is that inside defined atmosphere?
    if (atmosphere_dim > 1) {
      Vector ecef_toa;
      ecef_at_distance(ecef_toa, ecef, decef, l2toa);
      ecef2geodetic(pos_toa, ecef, refellipsoid);
      if (atmosphere_dim == 2 &&
          (pos_toa[1] < lat_grid[0] || pos_toa[1] > last(lat_grid))) {
        ARTS_USER_ERROR(
            "The sensor is above the top-of-atmosphere (TOA) altitude. "
            "The propagation path reaches TOA outside of covered latitude "
            "range. This is not allowed.\n"
            "         Sensor altitude: ",
            rte_pos[0],
            " m\n"
            "            TOA altitude: ",
            z_toa,
            "m\n"
            " Latitude of path at TOA: ",
            pos_toa[1],
            "\n"
            "           Latitude grid: ",
            lat_grid[0],
            "-",
            last(lat_grid));
      } else if (pos_toa[1] < lat_grid[0] || pos_toa[1] > last(lat_grid) ||
                 pos_toa[2] < lon_grid[0] || pos_toa[2] > last(lon_grid)) {
        ARTS_USER_ERROR(
            "The sensor is above the top-of-atmosphere (TOA) altitude. "
            "The propagation path reaches TOA outside of covered latitude-"
            "longitude domain. This is not allowed.\n"
            "          Sensor altitude: ",
            rte_pos[0],
            " m\n"
            "             TOA altitude: ",
            z_toa,
            "m\n"
            "  Latitude of path at TOA: ",
            pos_toa[1],
            "\n"
            "            Latitude grid: ",
            lat_grid[0],
            "-",
            last(lat_grid),
            " Longitide of path at TOA: ",
            pos_toa[2],
            "\n"
            "           Longitude grid: ",
            lon_grid[0],
            "-",
            last(lon_grid));
      }
    }

    // Do we hit a cloudbox boundary at TOA?
    if (cloudbox_on && cloudbox_limits[1] == z_grid.nelem() - 1) {
      if (atmosphere_dim == 1) {
        return PPATH_BACKGROUND_CLOUDBOX;
      } else if (atmosphere_dim == 2) {
        if (rte_pos[1] >= lat_grid[cloudbox_limits[2]] &&
            rte_pos[1] <= lat_grid[cloudbox_limits[3]])
          return PPATH_BACKGROUND_CLOUDBOX;
      } else {
        if (rte_pos[1] >= lat_grid[cloudbox_limits[2]] &&
            rte_pos[1] <= lat_grid[cloudbox_limits[3]] &&
            rte_pos[2] >= lon_grid[cloudbox_limits[4]] &&
            rte_pos[2] <= lon_grid[cloudbox_limits[5]])
          return PPATH_BACKGROUND_CLOUDBOX;
      }
    }
  }

  // -------------------------------------------------------------------------
  // Inside of atmosphere
  // -------------------------------------------------------------------------
  else {
    // Check if at boundary or inside of cloudbox
    if (cloudbox_on) {
      if (atmosphere_dim == 1) {
        if (rte_pos[0] >= z_grid[cloudbox_limits[0]] &&
            rte_pos[0] <= z_grid[cloudbox_limits[1]])
          return PPATH_BACKGROUND_CLOUDBOX;
      } else if (atmosphere_dim == 2) {
        if (rte_pos[0] >= z_grid[cloudbox_limits[0]] &&
            rte_pos[0] <= z_grid[cloudbox_limits[1]] &&
            rte_pos[1] >= lat_grid[cloudbox_limits[2]] &&
            rte_pos[1] <= lat_grid[cloudbox_limits[3]])
          return PPATH_BACKGROUND_CLOUDBOX;
      } else {
        if (rte_pos[0] >= z_grid[cloudbox_limits[0]] &&
            rte_pos[0] <= z_grid[cloudbox_limits[1]] &&
            rte_pos[1] >= lat_grid[cloudbox_limits[2]] &&
            rte_pos[1] <= lat_grid[cloudbox_limits[3]] &&
            rte_pos[2] >= lon_grid[cloudbox_limits[4]] &&
            rte_pos[2] <= lon_grid[cloudbox_limits[5]])
          return PPATH_BACKGROUND_CLOUDBOX;
      }
    }
  }

  return PPATH_BACKGROUND_UNDEFINED;
}
