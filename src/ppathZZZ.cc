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
 * @file   ppath_basic.cc
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2021-07-28
 *
 * @brief  Basic functions releated to calculation of propagation paths.
 *
 * The term propagation path is here shortened to ppath.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/


#include <algorithm>

#include "geodeticZZZ.h"
#include "ppathZZZ.h"
#include "variousZZZ.h"


/*===========================================================================
  === Internal functions, in alphabetical order
  ===========================================================================*/

/** Checks that a ppath is fully inside the atmoshere

   If not, an error is issued.

   @param[in]  ppath           As the WSV with the same name.
   @param[in]  atmosphere_dim  As the WSV with the same name.
   @param[in]  z_grid          As the WSV with the same name.
   @param[in]  lat_grid        As the WSV with the same name.
   @param[in]  lon_grid        As the WSV with the same name.

   @author Patrick Eriksson
   @date   2022-09-28
 */
void is_ppath_inside_atmosphere(const Ppath& ppath,
                                const Index& atmosphere_dim,
                                const Vector& z_grid,
                                const Vector& lat_grid,
                                const Vector& lon_grid) {
  if (atmosphere_dim > 1 && ppath.np > 1 &&
      is_pos_outside_atmosphere(ppath.pos(ppath.np - 1, joker),
                                atmosphere_dim,
                                last(z_grid),
                                lat_grid,
                                lon_grid)) {
    ConstVectorView pos = ppath.pos(ppath.np - 1, joker);
    if (atmosphere_dim == 2) {
      ARTS_USER_ERROR(
          "The propagation path is not fully inside the atmosphere.\n"
          "The model atmosphere has this extend\n"
          "  altitude: ",
          z_grid[0],
          "-",
          last(z_grid),
          "\n"
          "  latitude: ",
          lat_grid[0],
          "-",
          last(lat_grid),
          "\n"
          "while the position of *ppath* furthest away from the sensor\n"
          "is at (h,lat): (",
          pos[0],
          ",",
          pos[1],
          ")\n"
          "You need to extend the model atmosphere to cover this point.");
    } else {
      ARTS_USER_ERROR(
          "The propagation path is not fully inside the atmosphere.\n"
          "The model atmosphere has this extend\n"
          "  altitude: ",
          z_grid[0],
          "-",
          last(z_grid),
          "\n"
          "  latitude: ",
          lat_grid[0],
          "-",
          last(lat_grid),
          "\n"
          " longitude: ",
          lon_grid[0],
          "-",
          last(lon_grid),
          "\n"
          "while the position of *ppath* furthest away from the sensor\n"
          "is at (h,lat,lon): (",
          pos[0],
          ",",
          pos[1],
          ",",
          pos[2],
          ")\n"
          "You need to extend the model atmosphere to cover this point.");
    }
  }
}

/** Calculates grid positions of a ppath

   @param[in,out]  ppath      ppath-structure with at least fields dim, np,
                              pos, end_pos and start_pos set.
   @param[in]   z_grid        As the WSV with the same name.
   @param[in]   lat_grid      As the WSV with the same name.
   @param[in]   lon_grid      As the WSV with the same name.

   @author Patrick Eriksson
   @date   2021-08-16
 */
void ppath_calc_gp(Ppath& ppath,
                   const Vector& z_grid,
                   const Vector& lat_grid,
                   const Vector& lon_grid) {
  // Calculate grid positions
  ppath.gp_p.resize(ppath.np);
  if (ppath.np) {
    gridpos(ppath.gp_p, z_grid, ppath.pos(joker, 0));
  }
  if (ppath.dim >= 2 && ppath.np) {
    ppath.gp_lat.resize(ppath.np);
    gridpos(ppath.gp_lat, lat_grid, ppath.pos(joker, 1));
  } else {
    ppath.gp_lat.resize(0);
  }
  if (ppath.dim == 3 && ppath.np) {
    ppath.gp_lon.resize(ppath.np);
    gridpos(ppath.gp_lon, lon_grid, ppath.pos(joker, 2));
  } else {
    ppath.gp_lon.resize(0);
  }
}

/** Adjusts longitudes of a ppath

   If 3D, adjusts the longitudes in pos, start_pos and end_pos so they match
   lon_grid

   @param[in,out]  ppath      ppath-structure with at least fields dim, np,
                              pos, end_pos and start_pos set.
   @param[in]   lon_grid      As the WSV with the same name.

   @author Patrick Eriksson
   @date   2021-08-16
 */
void ppath_fix_lon(Ppath& ppath, const Vector& lon_grid) {
  if (ppath.dim == 3) {
    const Numeric lon_min = lon_grid[0];
    const Numeric lon_max = last(lon_grid);
    for (Index i = 0; i < ppath.np; i++) {
      ppath.pos(i, 2) = move_lon_to_range(ppath.pos(i, 2), lon_min, lon_max);
    }
    ppath.end_pos[2] = move_lon_to_range(ppath.end_pos[2], lon_min, lon_max);
    ppath.start_pos[2] =
        move_lon_to_range(ppath.start_pos[2], lon_min, lon_max);
  }
}

/** Some initial steps to determine a propagation path

   The function performs initial checks that are independent if the geometrical
   or refracted path shall be calculated.

   If sensor is outside of the atmosphere, the function checks that the
   propagation path enters the atmosphere from the top (or is totally in
   space).

   The function also checks if rte_pos is inside the cloudbox. Inside here
   includes to be at the boundary.

   No check with respect to the surface is made.

   @param[out]  l2toa             Length to TOA. Set to be negative if sensor
                                  is inside  of atmosphere.
   @param[out]  pos_toa           If l2toa >= 0, holds the position where the
                                  path enters the atmosphere.
   @param[out]  los_toa           If l2toa >= 0, holds the LOS at pos_toa
   @param[in]   rte_pos           As the WSV with the same name.
   @param[in]   rte_los           As the WSV with the same name.
   @param[in]   ecef              rte_pos in ECEF.
   @param[in]   decef             rte_los in ECEF.
   @param[in]   atmosphere_dim    As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with same name.
   @param[in]   z_grid            As the WSV with same name.
   @param[in]   lat_grid          As the WSV with same name.
   @param[in]   lon_grid          As the WSV with same name.
   @param[in]   cloudbox_on       As the WSV with same name.
   @param[in]   cloudbox_limits   As the WSV with same name.

   @return  Radiative background. Set to be undefined if not found to be space
   or cloudbox.

   @author Patrick Eriksson
   @date   2021-08-11
 */
enum PpathBackground ppath_init_calc(Numeric& l2toa,
                                     VectorView pos_toa,
                                     VectorView los_toa,
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
    if (rte_los[0] <= 90) {
      return PPATH_BACKGROUND_SPACE;
    }

    // Otherwise, do we reach TOA?
    l2toa = intersection_altitude(ecef, decef, refellipsoid, z_toa);
    if (l2toa < 0) {
      return PPATH_BACKGROUND_SPACE;
    } else {
      poslos_at_distance(pos_toa, los_toa, ecef, decef, refellipsoid, l2toa);
    }

    // If TOA reached, is that inside defined 2D or 3D atmosphere?
    if (atmosphere_dim > 1) {
      if (atmosphere_dim == 2) {
        if (pos_toa[1] < lat_grid[0] || pos_toa[1] > last(lat_grid))
          ARTS_USER_ERROR(
              "The sensor is above the top-of-atmosphere (TOA) altitude.\n"
              "The propagation path reaches TOA outside of covered\n"
              "latitude range. This is not allowed.\n"
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
                 !is_lon_in_range(pos_toa[2], lon_grid[0], last(lon_grid))) {
        ARTS_USER_ERROR(
            "The sensor is above the top-of-atmosphere (TOA) altitude.\n"
            "The propagation path reaches TOA outside of covered"
            " latitude-longitude domain. This is not allowed.\n"
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
            "\n"
            " Longitude of path at TOA: ",
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
        if (pos_toa[1] >= lat_grid[cloudbox_limits[2]] &&
            pos_toa[1] <= lat_grid[cloudbox_limits[3]])
          return PPATH_BACKGROUND_CLOUDBOX;
      } else {
        if (pos_toa[1] >= lat_grid[cloudbox_limits[2]] &&
            pos_toa[1] <= lat_grid[cloudbox_limits[3]] &&
            is_lon_in_range(pos_toa[2],
                            lon_grid[cloudbox_limits[4]],
                            lon_grid[cloudbox_limits[5]]))
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
            is_lon_in_range(rte_pos[2],
                            lon_grid[cloudbox_limits[4]],
                            lon_grid[cloudbox_limits[5]]))
          return PPATH_BACKGROUND_CLOUDBOX;
      }
    }
  }

  return PPATH_BACKGROUND_UNDEFINED;
}

/*===========================================================================
  === External functions, in alphabetical order
  ===========================================================================*/

Numeric find_crossing_cloudbox(const Vector rte_pos,
                               const Vector rte_los,
                               const Vector ecef,
                               const Vector decef,
                               const Index& atmosphere_dim,
                               const Vector& refellipsoid,
                               const Vector& z_grid,
                               const Vector& lat_grid,
                               const Vector& lon_grid,
                               const Index& cloudbox_on,
                               const ArrayOfIndex& cloudbox_limits,
                               const bool& is_outside) {
  if (!cloudbox_on) return -1;

  Numeric l2cbox = -1;

  // Outside of cloudbox
  // ----------------------------------------------------------------------
  if (is_outside) {
    // 1D
    if (atmosphere_dim == 1) {
      // Below, looking up?
      if (rte_pos[0] < z_grid[cloudbox_limits[0]] && rte_los[0] <= 90) {
        l2cbox = intersection_altitude(
            ecef, decef, refellipsoid, z_grid[cloudbox_limits[0]]);
        // Above, looking down?
      } else if (rte_pos[0] > z_grid[cloudbox_limits[1]] && rte_los[0] > 90) {
        l2cbox = intersection_altitude(
            ecef, decef, refellipsoid, z_grid[cloudbox_limits[1]]);
      }

      // 2D
    } else if (atmosphere_dim == 2) {
      Numeric lt = -1;  // Test length
      Vector pt(3);     // Test position
      // Intersection with low altitude face?
      if (rte_pos[0] < z_grid[cloudbox_limits[0]] && rte_los[0] <= 90) {
        lt = intersection_altitude(
            ecef, decef, refellipsoid, z_grid[cloudbox_limits[0]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[0] - z_grid[cloudbox_limits[0]]) < 0.1);
          if (pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]])
            l2cbox = lt;
        }
        // Intersection with high altitude face?
      } else if (rte_pos[0] > z_grid[cloudbox_limits[1]] && rte_los[0] > 90) {
        lt = intersection_altitude(
            ecef, decef, refellipsoid, z_grid[cloudbox_limits[1]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[0] - z_grid[cloudbox_limits[1]]) < 0.1);
          if (pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]])
            l2cbox = lt;
        }
      }
      // Intersection with south face?
      if (rte_pos[1] < lat_grid[cloudbox_limits[2]] && abs(rte_los[1]) < 90) {
        lt = intersection_latitude(ecef,
                                   decef,
                                   rte_pos,
                                   rte_los,
                                   refellipsoid,
                                   lat_grid[cloudbox_limits[2]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[1] - lat_grid[cloudbox_limits[2]]) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]])
            l2cbox = min_geq(l2cbox, lt, 0);
        }
        // Intersection with north face?
      } else if (rte_pos[1] > lat_grid[cloudbox_limits[3]] &&
                 abs(rte_los[1]) >= 90) {
        lt = intersection_latitude(ecef,
                                   decef,
                                   rte_pos,
                                   rte_los,
                                   refellipsoid,
                                   lat_grid[cloudbox_limits[3]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[1] - lat_grid[cloudbox_limits[3]]) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]])
            l2cbox = min_geq(l2cbox, lt, 0);
        }
      }

      // 3D
    } else {
      Numeric lt = -1;  // Test length
      Vector pt(3);     // Test position
      // Intersection with low altitude face?
      if (rte_pos[0] < z_grid[cloudbox_limits[0]] && rte_los[0] <= 90) {
        lt = intersection_altitude(
            ecef, decef, refellipsoid, z_grid[cloudbox_limits[0]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[0] - z_grid[cloudbox_limits[0]]) < 0.1);
          if (pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]] &&
              is_lon_in_range(pt[2],
                              lon_grid[cloudbox_limits[4]],
                              lon_grid[cloudbox_limits[5]]))
            l2cbox = lt;
        }
        // Intersection with high altitude face?
      } else if (rte_pos[0] > z_grid[cloudbox_limits[1]] && rte_los[0] > 90) {
        lt = intersection_altitude(
            ecef, decef, refellipsoid, z_grid[cloudbox_limits[1]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[0] - z_grid[cloudbox_limits[1]]) < 0.1);
          if (pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]] &&
              is_lon_in_range(pt[2],
                              lon_grid[cloudbox_limits[4]],
                              lon_grid[cloudbox_limits[5]]))
            l2cbox = lt;
        }
      }
      // Intersection with south face?
      if (rte_pos[1] < lat_grid[cloudbox_limits[2]] && abs(rte_los[1]) < 90) {
        lt = intersection_latitude(ecef,
                                   decef,
                                   rte_pos,
                                   rte_los,
                                   refellipsoid,
                                   lat_grid[cloudbox_limits[2]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[1] - lat_grid[cloudbox_limits[2]]) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]] &&
              is_lon_in_range(pt[2],
                              lon_grid[cloudbox_limits[4]],
                              lon_grid[cloudbox_limits[5]]))
            l2cbox = min_geq(l2cbox, lt, 0);
        }
        // Intersection with north face?
      } else if (rte_pos[1] > lat_grid[cloudbox_limits[3]] &&
                 abs(rte_los[1]) >= 90) {
        lt = intersection_latitude(ecef,
                                   decef,
                                   rte_pos,
                                   rte_los,
                                   refellipsoid,
                                   lat_grid[cloudbox_limits[3]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[1] - lat_grid[cloudbox_limits[3]]) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]] &&
              is_lon_in_range(pt[2],
                              lon_grid[cloudbox_limits[4]],
                              lon_grid[cloudbox_limits[5]]))
            l2cbox = min_geq(l2cbox, lt, 0);
        }
      }
      // For longitude faces we need to handle [-180,180] vs [0,360]
      const Numeric rte_lon_adjusted =
          move_lon_to_range(rte_pos[2],
                            lon_grid[cloudbox_limits[4]],
                            lon_grid[cloudbox_limits[5]]);
      Vector rte_pos_adjusted = rte_pos;
      rte_pos_adjusted[2] = rte_lon_adjusted;
      // Intersection with west face?
      if (rte_lon_adjusted < lon_grid[cloudbox_limits[4]] && rte_los[1] > 0) {
        lt = intersection_longitude(ecef,
                                    decef,
                                    rte_pos_adjusted,
                                    rte_los,
                                    lon_grid[cloudbox_limits[4]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[2] - shift_lon_to_pm180(
                                      lon_grid[cloudbox_limits[4]])) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]] &&
              pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]])
            l2cbox = min_geq(l2cbox, lt, 0);
        }
        // Intersection with east face?
      } else if (rte_lon_adjusted > lon_grid[cloudbox_limits[5]] &&
                 rte_los[1] < 0) {
        lt = intersection_longitude(ecef,
                                    decef,
                                    rte_pos_adjusted,
                                    rte_los,
                                    lon_grid[cloudbox_limits[5]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT(abs(pt[2] - shift_lon_to_pm180(
                                      lon_grid[cloudbox_limits[5]])) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]] &&
              pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]])
            l2cbox = min_geq(l2cbox, lt, 0);
        }
      }
    }

    // Inside of cloudbox
    // ----------------------------------------------------------------------
  } else {
    // Not yet implemented
    ARTS_ASSERT(0);
  }

  return l2cbox;
}

Index is_pos_outside_atmosphere(const Vector pos,
                                const Index& atmosphere_dim,
                                const Numeric& z_toa,
                                const Vector& lat_grid,
                                const Vector& lon_grid) {
  if (pos[0] > z_toa) return 1;
  if (atmosphere_dim >= 2) {
    if (pos[1] < lat_grid[0] || pos[1] > last(lat_grid)) return 2;
    if (atmosphere_dim == 3) {
      if (!is_lon_in_range(pos[2], lon_grid[0], last(lon_grid))) return 3;
    }
  }
  return 0;
}

void ppath_geom_const_lstep(Ppath& ppath,
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
                            const Numeric& l_step_max,
                            const Numeric& l_total_max,
                            const Numeric& l_accuracy,
                            const Index& safe_surface_search,
                            const Index& do_not_calc_gps) {
  // Convert rte_pos/los to ECEF
  Vector ecef(3), decef(3);
  geodetic_los2ecef(ecef, decef, rte_pos, rte_los, refellipsoid);

  // Number of points of ppath (-1 signifies still not known)
  Index np = -1;

  // Some initial checks and calculations
  Numeric l2toa;
  Vector pos_toa(3), los_toa(2);
  enum PpathBackground background = ppath_init_calc(l2toa,
                                                    pos_toa,
                                                    los_toa,
                                                    rte_pos,
                                                    rte_los,
                                                    ecef,
                                                    decef,
                                                    atmosphere_dim,
                                                    refellipsoid,
                                                    z_grid,
                                                    lat_grid,
                                                    lon_grid,
                                                    cloudbox_on,
                                                    cloudbox_limits);
  if (background == PPATH_BACKGROUND_SPACE) np = 0;
  if (background == PPATH_BACKGROUND_CLOUDBOX) np = 1;
  // Length from sensor to start of ppath inside the atmosphere
  const Numeric l_end = l2toa > 0 ? l2toa : 0;

  // If np < 0 we have a path to determine
  Numeric l_start = -1;
  if (np < 0) {
    // Length from sensor to end of ppath inside the atmosphere or
    // to tangent point if it's closer
    Numeric l_test;

    // Distance along ppath to the surface (negative if no intersection)
    // The function also checks that rte_pos is above the surface
    l_test = find_crossing_with_surface_z(rte_pos,
                                          rte_los,
                                          ecef,
                                          decef,
                                          atmosphere_dim,
                                          refellipsoid,
                                          surface_elevation,
                                          l_accuracy,
                                          safe_surface_search);
    if (l_test >= 0) {
      l_start = l_test;
      background = PPATH_BACKGROUND_SURFACE;
    }

    // Determine length to cloudbox (negative if no intersection)
    l_test = find_crossing_cloudbox(rte_pos,
                                    rte_los,
                                    ecef,
                                    decef,
                                    atmosphere_dim,
                                    refellipsoid,
                                    z_grid,
                                    lat_grid,
                                    lon_grid,
                                    cloudbox_on,
                                    cloudbox_limits,
                                    true);
    if (l_test >= 0 && (l_start < 0 || l_test < l_start)) {
      l_start = l_test;
      background = PPATH_BACKGROUND_CLOUDBOX;
    }

    // If still no intersection, l_start must be at TOA
    if (l_start < 0) {
      if (l2toa < 0) {
        // Observation from within the atmosphere
        l_start =
            intersection_altitude(ecef, decef, refellipsoid, last(z_grid));
      } else {
        // Observation from outside the atmosphere
        Vector ecef_toa(3);
        ecef_at_distance(ecef_toa, ecef, decef, l2toa);
        // Ignore lengths < 1m to find exit point, and not the entrance point
        // from which we start
        l_start = l2toa + intersection_altitude(
                              ecef_toa, decef, refellipsoid, last(z_grid), 1);
      }
      background = PPATH_BACKGROUND_SPACE;
    }
    // Apply total length criterion?
    if (l_total_max > 0) {
      l_start = l2toa > 0 ? l2toa + l_total_max : l_total_max;
      background = PPATH_BACKGROUND_STOP_DISTANCE;
    }
    ARTS_ASSERT(background);

    // Determine np
    ARTS_ASSERT(l_start >= l_end);
    np = 1 + Index(ceil((l_start - l_end) / l_step_max));
  }

  // Fill ppath
  ppath.dim = atmosphere_dim;
  ppath.np = np;
  ppath.backgroundZZZ = background;
  ppath.start_lstep = 0;
  ppath.start_pos = rte_pos;
  ppath.start_los = rte_los;
  ppath.start_lstep = l_end;
  ppath.nreal = Vector(np, 1.0);
  ppath.ngroup = Vector(np, 1.0);
  ppath.pos.resize(np, 3);
  ppath.los.resize(np, 2);
  //
  if (np == 0) {
    ppath.lstep.resize(0);
  } else if (np == 1) {
    ppath.lstep.resize(0);
    if (l2toa >= 0) {
      ppath.pos(0, joker) = pos_toa;
      ppath.los(0, joker) = los_toa;
    } else {
      ppath.pos(0, joker) = rte_pos;
      ppath.los(0, joker) = rte_los;
    }
  } else {
    // Create equidistant length vector and fill pos and los
    Vector l;
    nlinspace(l, l_end, l_start, np);
    for (Index i = 0; i < np; i++) {
      poslos_at_distance(ppath.pos(i, joker),
                         ppath.los(i, joker),
                         ecef,
                         decef,
                         refellipsoid,
                         l[i]);
    }
    ppath.lstep.resize(np - 1);
    ppath.lstep = l[1] - l[0];
  }
  if (np == 0) {
    ppath.end_los = ppath.start_los;
  } else {
    ppath.end_pos = ppath.pos(ppath.np - 1, joker);
    ppath.end_los = ppath.los(ppath.np - 1, joker);

    // Adjust longitudes
    ppath_fix_lon(ppath, lon_grid);

    // Check that ppath is fully inside the atmosphere
    is_ppath_inside_atmosphere(
        ppath, atmosphere_dim, z_grid, lat_grid, lon_grid);

    // Calculate grid positions?
    if (!do_not_calc_gps) {
      ppath_calc_gp(ppath, z_grid, lat_grid, lon_grid);
    }
  }
  if (np == 0 || do_not_calc_gps) {
    ppath.gp_p.resize(0);
    ppath.gp_lat.resize(0);
    ppath.gp_lon.resize(0);
  }
}

void ppath_grid_crossings(Ppath& ppath,
                          const Index& atmosphere_dim,
                          const Vector& refellipsoid,
                          const Vector& z_grid,
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Numeric& l_step_max,
                          const Numeric& l_accuracy,
                          const Index& do_not_calc_gps) {
  // Nothing to do if there is no ppath step
  if (ppath.np < 2) {
    return;
  }

  // l means distance from ppath pos[0]
  // dl means distance from some other ppath point
  // Excpetion: l_step_max is still a local length

  // Process ppath to get some data in another form
  Vector l_acc_ppath(ppath.np);  // Accumulated length along ppath
  Vector gp_z(ppath.np);         // Grid positions as Numeric, i.e. idx+fd[0]
  Vector gp_lat(atmosphere_dim < 2 ? 0 : ppath.np);
  Vector gp_lon(atmosphere_dim < 3 ? 0 : ppath.np);
  for (Index ip = 0; ip < ppath.np; ++ip) {
    if (ip == 0) {
      l_acc_ppath[ip] = 0;
    } else {
      l_acc_ppath[ip] = l_acc_ppath[ip - 1] + ppath.lstep[ip - 1];
    }
    gp_z[ip] = Numeric(ppath.gp_p[ip].idx) + ppath.gp_p[ip].fd[0];
    if (atmosphere_dim >= 2) {
      gp_lat[ip] = Numeric(ppath.gp_lat[ip].idx) + ppath.gp_lat[ip].fd[0];
      if (atmosphere_dim >= 3) {
        gp_lon[ip] = Numeric(ppath.gp_lon[ip].idx) + ppath.gp_lon[ip].fd[0];
      }
    }
  }

  // Fixes to handle numerical inaccuracy:
  // Make sure that start gp_z is inside the atmosphere
  if (gp_z[0] >= Numeric(z_grid.nelem() - 1)) {
    gp_z[0] = Numeric(z_grid.nelem() - 1) - 1e-9;
  }
  // Total length of ppath, minus a small distance to avoid that end point gets repeated
  const Numeric l2end = l_acc_ppath[ppath.np - 1] - l_accuracy / 100.0;

  // Containers for new ppath points (excluding start and end points, that
  // always are taken from original ppath)
  ArrayOfIndex istart_array(0);
  ArrayOfNumeric l_array(0);

  // Loop ppath steps to find grid crossings
  Numeric l_last_inserted = 0;
  for (Index ip = 0; ip < ppath.np - 1; ++ip) {
    // Length to grid crossings inside ppath step
    ArrayOfNumeric dl_from_ip(0);

    // Change in integer grid position for each dimension
    const Index dgp_z = n_int_between(gp_z[ip], gp_z[ip + 1]);
    const Index dgp_lat =
        atmosphere_dim < 2 ? 0 : n_int_between(gp_lat[ip], gp_lat[ip + 1]);
    const Index dgp_lon =
        atmosphere_dim < 3 ? 0 : n_int_between(gp_lon[ip], gp_lon[ip + 1]);

    if (dgp_z || dgp_lat || dgp_lon) {
      // ECEF at start end of ppath step
      Vector ecef(3), decef(3);
      geodetic_los2ecef(ecef,
                        decef,
                        ppath.pos(ip, joker),
                        ppath.los(ip, joker),
                        refellipsoid);

      // Crossing(s) of z_grid
      for (Index i = 1; i <= abs(dgp_z); ++i) {
        const Numeric dl_test = intersection_altitude(
            ecef,
            decef,
            refellipsoid,
            z_grid[int_at_step(gp_z[ip], sign(dgp_z) * i)]);
        if (dl_test > 0 && l_acc_ppath[ip] + dl_test < l2end) {
          dl_from_ip.push_back(dl_test);
        }
      }

      // Crossing(s) of lat_grid
      for (Index i = 1; i <= abs(dgp_lat); ++i) {
        const Numeric dl_test = intersection_latitude(
            ecef,
            decef,
            ppath.pos(ip, joker),
            ppath.los(ip, joker),
            refellipsoid,
            lat_grid[int_at_step(gp_lat[ip], sign(dgp_lat) * i)]);
        if (dl_test > 0 && l_acc_ppath[ip] + dl_test < l2end) {
          dl_from_ip.push_back(dl_test);
        }
      }

      // Crossing(s) of lon_grid
      for (Index i = 1; i <= abs(dgp_lon); ++i) {
        const Numeric dl_test = intersection_longitude(
            ecef,
            decef,
            ppath.pos(ip, joker),
            ppath.los(ip, joker),
            lon_grid[int_at_step(gp_lon[ip], sign(dgp_lon) * i)]);
        if (dl_test > 0 && l_acc_ppath[ip] + dl_test < l2end) {
          dl_from_ip.push_back(dl_test);
        }
      }

      // Sort dl_from_ip
      std::sort(dl_from_ip.begin(), dl_from_ip.end());

      // Move to overall arrays and add points if l_step_max that requires
      for (Index i = 0; i < dl_from_ip.nelem(); ++i) {
        // Some useful lengths
        const Numeric l_next = l_acc_ppath[ip] + dl_from_ip[i];
        const Numeric dl = l_next - l_last_inserted;
        // Number of points needed to fulfill l_step_max
        if (dl > l_step_max) {
          const Index n_extra = Index(std::floor(dl / l_step_max));
          const Numeric dl_step = dl / Numeric(n_extra + 1);
          for (Index extra = 0; extra < n_extra; ++extra) {
            const Numeric l_extra = l_last_inserted + dl_step;
            istart_array.push_back(ip);
            l_array.push_back(l_extra);
            l_last_inserted = l_extra;
          }
        }
        // Add grid crossing point
        istart_array.push_back(ip);
        l_array.push_back(l_next);
        l_last_inserted = l_next;
      }
    }  // if dgp_p
  }    // ip loop

  // The distance between last grid crossing and end point can exceed l_step_max
  // Fix (largely same code as above)!
  const Numeric dl = l_acc_ppath[ppath.np - 1] - l_last_inserted;
  if (dl > l_step_max) {
    const Index n_extra = Index(std::floor(dl / l_step_max));
    const Numeric dl_step = dl / Numeric(n_extra + 1);
    for (Index extra = 0; extra < n_extra; ++extra) {
      const Numeric l_extra = l_last_inserted + dl_step;
      istart_array.push_back(ppath.np - 2);
      l_array.push_back(l_extra);
      l_last_inserted = l_extra;
    }
  }
  //-------------------------------------------------------------

  // Make copies of data in ppath that will change, but we need
  Index np = ppath.np;
  Vector nreal = ppath.nreal;
  Vector ngroup = ppath.ngroup;
  Matrix pos = ppath.pos;
  Matrix los = ppath.los;

  // New size of ppath
  const Index nl = l_array.nelem();

  ppath.np = nl + 2;
  ppath.nreal = Vector(ppath.np, 1.0);  // We guess on no refraction
  ppath.ngroup = Vector(ppath.np, 1.0);
  ppath.lstep.resize(ppath.np - 1);
  ppath.pos.resize(ppath.np, 3);
  ppath.los.resize(ppath.np, 2);

  // Pos and los at end points
  ppath.pos(0, joker) = pos(0, joker);
  ppath.los(0, joker) = los(0, joker);
  ppath.pos(ppath.np - 1, joker) = pos(np - 1, joker);
  ppath.los(ppath.np - 1, joker) = los(np - 1, joker);

  // Calculate and insert new pos and los, and do lstep in parallel
  Vector l_array_as_vector(nl);
  if (nl) {
    for (Index i = 0; i < nl; ++i) {
      l_array_as_vector[i] = l_array[i];
      Vector ecef(3), decef(3);
      geodetic_los2ecef(ecef,
                        decef,
                        pos(istart_array[i], joker),
                        los(istart_array[i], joker),
                        refellipsoid);
      poslos_at_distance(ppath.pos(i + 1, joker),
                         ppath.los(i + 1, joker),
                         ecef,
                         decef,
                         refellipsoid,
                         l_array[i] - l_acc_ppath[istart_array[i]]);
      if (i > 0) {
        ppath.lstep[i] = l_array[i] - l_array[i - 1];
        //ARTS_ASSERT(ppath.lstep[i] > 0)
      }
    }
    ppath.lstep[0] = l_array[0];
    ppath.lstep[nl] = l_acc_ppath[np - 1] - l_array[nl - 1];
    ARTS_ASSERT(ppath.lstep[0] > 0)
    ARTS_ASSERT(ppath.lstep[nl] > 0)
  } else {
    ppath.lstep[0] = l_acc_ppath[np - 1];
  }

  // New refractive indices, mainly set by interpolation
  if (max(nreal) > 1.0) {
    ppath.nreal[0] = nreal[0];
    ppath.ngroup[0] = ngroup[0];
    ppath.nreal[ppath.np - 1] = nreal[np - 1];
    ppath.ngroup[ppath.np - 1] = ngroup[np - 1];
    //
    if (nl) {
      ArrayOfGridPos gp(nl);
      gridpos(gp, l_acc_ppath, l_array_as_vector);
      Matrix itw(nl, 2);
      interpweights(itw, gp);
      interp(ppath.nreal[Range(1, nl)], itw, nreal, gp);
      interp(ppath.ngroup[Range(1, nl)], itw, ngroup, gp);
    }
  }

  // Adjust longitudes
  ppath_fix_lon(ppath, lon_grid);

  // Calculate grid positions?
  if (do_not_calc_gps) {
    ppath.gp_p.resize(0);
    ppath.gp_lat.resize(0);
    ppath.gp_lon.resize(0);
  } else {
    ppath_calc_gp(ppath, z_grid, lat_grid, lon_grid);
  }
}
