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

#include "geodeticZZZ.h"
#include "interpolation.h"
#include "jacobian.h"
#include "ppathZZZ.h"
#include "math_funcs.h"

using GriddedFieldGrids::GFIELD2_LAT_GRID;
using GriddedFieldGrids::GFIELD2_LON_GRID;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);


/*===========================================================================
  === The functions, in alphabetical order
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
  if (!cloudbox_on)
    return -1;

  Numeric l2cbox = -1;
  
  // Outside of cloudbox
  // ----------------------------------------------------------------------
  if (is_outside) {

    // 1D
    if (atmosphere_dim == 1) {
      // Below, looking up?
      if (rte_pos[0] < z_grid[cloudbox_limits[0]] && rte_los[0] <= 90) {
        l2cbox = intersection_altitude(ecef, decef, refellipsoid,
                                       z_grid[cloudbox_limits[0]]);
      // Above, looking down?
      } else if (rte_pos[0] > z_grid[cloudbox_limits[1]] && rte_los[0] > 90) {
        l2cbox = intersection_altitude(ecef, decef, refellipsoid,
                                       z_grid[cloudbox_limits[1]]);
      }
      
    // 2D
    } else if (atmosphere_dim == 2) {
      Numeric lt = -1;  // Test length
      Vector pt(3);     // Test position
      // Intersection with low altitude face?
      if (rte_pos[0] < z_grid[cloudbox_limits[0]] && rte_los[0] <= 90) {
        lt = intersection_altitude(ecef, decef, refellipsoid,
                                   z_grid[cloudbox_limits[0]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[0]-z_grid[cloudbox_limits[0]]) < 0.1);
          if (pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]])
            l2cbox = lt;
        }
      // Intersection with high altitude face?
      } else if (rte_pos[0] > z_grid[cloudbox_limits[1]] && rte_los[0] > 90) {
        lt = intersection_altitude(ecef, decef, refellipsoid,
                                   z_grid[cloudbox_limits[1]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[0]-z_grid[cloudbox_limits[1]]) < 0.1);
          if (pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]])
            l2cbox = lt;
        }
      }
      // Intersection with south face?
      if (rte_pos[1] < lat_grid[cloudbox_limits[2]] && abs(rte_los[1]) < 90) {
        lt = intersection_latitude(ecef, decef, rte_pos, rte_los, refellipsoid,
                                   lat_grid[cloudbox_limits[2]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[1]-lat_grid[cloudbox_limits[2]]) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]])
            l2cbox = min_geq(l2cbox, lt, 0);
        }
      // Intersection with north face?
      } else if (rte_pos[1] > lat_grid[cloudbox_limits[3]] && abs(rte_los[1]) >= 90) {
        lt = intersection_latitude(ecef, decef, rte_pos, rte_los, refellipsoid,
                                   lat_grid[cloudbox_limits[3]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[1]-lat_grid[cloudbox_limits[3]]) < 0.1);
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
        lt = intersection_altitude(ecef, decef, refellipsoid,
                                   z_grid[cloudbox_limits[0]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[0]-z_grid[cloudbox_limits[0]]) < 0.1);
          if (pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]] &&
              is_lon_in_range(pt[2],
                              lon_grid[cloudbox_limits[4]],
                              lon_grid[cloudbox_limits[5]]))
            l2cbox = lt;
        }
      // Intersection with high altitude face?
      } else if (rte_pos[0] > z_grid[cloudbox_limits[1]] && rte_los[0] > 90) {
        lt = intersection_altitude(ecef, decef, refellipsoid,
                                   z_grid[cloudbox_limits[1]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[0]-z_grid[cloudbox_limits[1]]) < 0.1);
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
        lt = intersection_latitude(ecef, decef, rte_pos, rte_los, refellipsoid,
                                   lat_grid[cloudbox_limits[2]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[1]-lat_grid[cloudbox_limits[2]]) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]] &&
              is_lon_in_range(pt[2],
                              lon_grid[cloudbox_limits[4]],
                              lon_grid[cloudbox_limits[5]]))
            l2cbox = min_geq(l2cbox, lt, 0);
        }
      // Intersection with north face?
      } else if (rte_pos[1] > lat_grid[cloudbox_limits[3]] && abs(rte_los[1]) >= 90) {
        lt = intersection_latitude(ecef, decef, rte_pos, rte_los, refellipsoid,
                                   lat_grid[cloudbox_limits[3]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[1]-lat_grid[cloudbox_limits[3]]) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]] &&
              is_lon_in_range(pt[2],
                              lon_grid[cloudbox_limits[4]],
                              lon_grid[cloudbox_limits[5]]))
            l2cbox = min_geq(l2cbox, lt, 0);
        }
      }
      // For longitude faces we need to handle [-180,180] vs [0,360]
      const Numeric rte_lon_adjusted = move_lon_to_range(rte_pos[2],
                                                         lon_grid[cloudbox_limits[4]],
                                                         lon_grid[cloudbox_limits[5]]);
      Vector rte_pos_adjusted = rte_pos;
      rte_pos_adjusted[2] = rte_lon_adjusted;
      // Intersection with west face?
      if (rte_lon_adjusted < lon_grid[cloudbox_limits[4]] && rte_los[1] > 0) {
        lt = intersection_longitude(ecef, decef, rte_pos_adjusted, rte_los,
                                    lon_grid[cloudbox_limits[4]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[2]-
                           shift_lon_to_pm180(lon_grid[cloudbox_limits[4]])) < 0.1);
          if (pt[0] >= z_grid[cloudbox_limits[0]] &&
              pt[0] <= z_grid[cloudbox_limits[1]] &&
              pt[1] >= lat_grid[cloudbox_limits[2]] &&
              pt[1] <= lat_grid[cloudbox_limits[3]])
            l2cbox = min_geq(l2cbox, lt, 0);
        }
      // Intersection with east face?
      } else if (rte_lon_adjusted > lon_grid[cloudbox_limits[5]] && rte_los[1] < 0) {
        lt = intersection_longitude(ecef, decef, rte_pos_adjusted, rte_los,
                                    lon_grid[cloudbox_limits[5]]);
        if (lt >= 0) {
          pos_at_distance(pt, ecef, decef, refellipsoid, lt);
          ARTS_ASSERT (abs(pt[2]-
                           shift_lon_to_pm180(lon_grid[cloudbox_limits[5]])) < 0.1);
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
    ARTS_ASSERT (0);
  }
  
  return l2cbox;
}

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
    z_min = surface_elevation.data(0,0);
    z_max = surface_elevation.data(0,0);    
  } else {
    z_min = min(surface_elevation.data);
    z_max = max(surface_elevation.data);
  }

  // Catch upward looking cases that can not have a surface intersection
  if (rte_pos[0]>=z_max && rte_los[0]<=90) {
    return -1;
  }

  // Check that observation position is above ground
  if (rte_pos[0] < z_max) {
    Numeric z_surf = surface_z_at_pos(rte_pos, atmosphere_dim, surface_elevation);
    if (rte_pos[0] < z_surf-l_accuracy)
      ARTS_USER_ERROR("The sensor is below the surface. Not allowed!\n"
                      "The sensor altitude is at ", rte_pos[0], " m\n"
                      "The surface altitude is ", z_surf, " m\n"
                      "The position is (lat,lon): (", rte_pos[1], ",",
                      rte_pos[2], ")");
  }

  // Constant surface altitude (in comparison to *l_accuracy*)
  if (atmosphere_dim == 1 || z_max-z_min < l_accuracy/100) {
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
      if (l_min < 0)  
        return -1;
    }
    // Find max distance for search. 
    // If below z_max and upward, given by z_max
    // Otherwise in general given by z_min. If z_min not reached, the distance
    // is instead given by tangent point
    Numeric l_max;
    bool l_max_could_be_above_surface = false;
    if (rte_pos[0]<=z_max && rte_los[0]<=90) {
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
      Numeric l_test = l_min - l_accuracy/2;  // Remove l/2 to get exact result 
      bool above_surface = true;              // if true l_test is 0
      while (above_surface && l_test<l_max) {
        l_test += l_accuracy;
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_test);
        Numeric z_surf = surface_z_at_pos(pos, atmosphere_dim, surface_elevation);
        if (pos[0] < z_surf)
          above_surface = false;
      }
      if (above_surface) {
        return -1;
      } else {
        return l_test - l_accuracy/2;
      }
      
    // Bisection search
    // ----------------------
    } else {
      // If l_max matches a point above the surface, we have no intersection
      // according to this search algorithm. And the search fails. So we need
      // to check that point if status unclear
      if (l_max_could_be_above_surface) {
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_max);
        Numeric z_surf = surface_z_at_pos(pos, atmosphere_dim, surface_elevation);
        if (pos[0] > z_surf)
          return -1;
      }
      // Start bisection
      while (l_max-l_min > 2*l_accuracy) {
        const Numeric l_test = (l_min+l_max)/2;
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_test);
        Numeric z_surf = surface_z_at_pos(pos, atmosphere_dim, surface_elevation);
        if (pos[0] >= z_surf)
          l_min = l_test;
        else
          l_max = l_test;
      }
      return (l_min+l_max)/2;
    }
  }
}



Index is_pos_outside_atmosphere(const Vector pos,
                                const Index& atmosphere_dim,
                                const Numeric& z_toa,
                                const Vector& lat_grid,
                                const Vector& lon_grid) {
  if (pos[0] > z_toa)
    return 1;   
  if (atmosphere_dim >= 2) {
    if (pos[1] < lat_grid[0] || pos[1] > last(lat_grid))
      return 2;  
    if (atmosphere_dim == 3) {
      if (!is_lon_in_range(pos[2], lon_grid[0], last(lon_grid)))
        return 3;
    }
  }
  return 0;
}



void ppath_fix_lon_and_gp(Ppath& ppath,
                          const Vector& z_grid,
                          const Vector& lat_grid,
                          const Vector& lon_grid) {
  // Adjust longitudes?
  if (ppath.dim == 3) {
    const Numeric lon_min = lon_grid[0];
    const Numeric lon_max = last(lon_grid);
    for (Index i=0; i<ppath.np; i++) {
      ppath.pos(i,2) = move_lon_to_range(ppath.pos(i,2), lon_min, lon_max);
    }
      ppath.end_pos[2] = move_lon_to_range(ppath.end_pos[2], lon_min, lon_max);
      ppath.start_pos[2] = move_lon_to_range(ppath.start_pos[2], lon_min, lon_max);
  }
  // Calculate grid positions
  ppath.gp_p.resize(ppath.np);
  if (ppath.np)
    gridpos(ppath.gp_p, z_grid, ppath.pos(joker,0));
  if (ppath.dim >= 2 && ppath.np) {
    ppath.gp_lat.resize(ppath.np);
    gridpos(ppath.gp_lat, lat_grid, ppath.pos(joker,1));
  } else {
    ppath.gp_lat.resize(0);
  }
  if (ppath.dim == 3  && ppath.np) {
    ppath.gp_lon.resize(ppath.np);
    gridpos(ppath.gp_lon, lon_grid, ppath.pos(joker,2));
  } else {
    ppath.gp_lon.resize(0);
  }
}



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
  const Index outside = is_pos_outside_atmosphere(rte_pos,
                                                 atmosphere_dim,
                                                 z_toa,
                                                 lat_grid,
                                                 lon_grid);

  // -------------------------------------------------------------------------
  // Outside of atmosphere
  // -------------------------------------------------------------------------
  if (outside) {
    // Outside but sensor below TOA
    if (outside == 2) {
      ARTS_USER_ERROR("The sensor is below the top-of-atmosphere (TOA) altitude,\n"
                      "but outside in the latitude range. This is not allowed.\n"
                      " Sensor altitude: ", rte_pos[0], " m\n"
                      "    TOA altitude: ", z_toa, "m\n"
                      " Sensor latitude: ", rte_pos[1], "\n"
                      "   Latitude grid: ", lat_grid[0], "-", last(lat_grid));
    }
    if (outside == 3) {
      ARTS_USER_ERROR("The sensor is below the top-of-atmosphere (TOA) altitude,\n"
                      "but outside in the longitude range. This is not allowed.\n"
                      "  Sensor altitude: ", rte_pos[0], " m\n"
                      "     TOA altitude: ", z_toa, "m\n"
                      " Sensor longitude: ", rte_pos[2], "\n"
                      "   Longitude grid: ", lon_grid[0], "-", last(lon_grid));
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
      if (atmosphere_dim==2) {
        if (pos_toa[1]<lat_grid[0] || pos_toa[1]>last(lat_grid))
          ARTS_USER_ERROR("The sensor is above the top-of-atmosphere (TOA) altitude.\n"
                          "The propagation path reaches TOA outside of covered\n"
                          "latitude range. This is not allowed.\n"
                          "         Sensor altitude: ", rte_pos[0], " m\n"
                          "            TOA altitude: ", z_toa, "m\n"
                          " Latitude of path at TOA: ", pos_toa[1], "\n"
                          "           Latitude grid: ", lat_grid[0], "-",
                          last(lat_grid));
      } else if (pos_toa[1]<lat_grid[0] || pos_toa[1]>last(lat_grid) ||
                 !is_lon_in_range(pos_toa[2], lon_grid[0], last(lon_grid))) {
        ARTS_USER_ERROR("The sensor is above the top-of-atmosphere (TOA) altitude.\n"
                        "The propagation path reaches TOA outside of covered"
                        " latitude-longitude domain. This is not allowed.\n"
                        "          Sensor altitude: ", rte_pos[0], " m\n"
                        "             TOA altitude: ", z_toa, "m\n"
                        "  Latitude of path at TOA: ", pos_toa[1], "\n"
                        "            Latitude grid: ", lat_grid[0], "-",
                        last(lat_grid), "\n"
                        " Longitude of path at TOA: ", pos_toa[2], "\n"
                        "           Longitude grid: ", lon_grid[0], "-", last(lon_grid));
      }
    }

    // Do we hit a cloudbox boundary at TOA?
    if (cloudbox_on && cloudbox_limits[1] == z_grid.nelem()-1) {
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



Numeric surface_z_at_pos(const Vector pos,
                         const Index& atmosphere_dim,
                         const GriddedField2& surface_elevation) {
  // We want here an extrapolation to infinity ->
  //                                        extremly high extrapolation factor
  const Numeric inf_proxy = 1.0e99;
  
  if (atmosphere_dim == 1) {
    return surface_elevation.data(0,0);
    
  } else if (atmosphere_dim == 2) {
    const Vector& lat_grid = surface_elevation.get_numeric_grid(GFIELD2_LAT_GRID);
    ArrayOfGridPos gp_lat(1);
    if (lat_grid.nelem() > 1) {
      gridpos(gp_lat, lat_grid, pos[1], inf_proxy);
      jacobian_type_extrapol(gp_lat);
    } else {
      gp4length1grid(gp_lat);
    }
    Vector itw(2);
    interpweights(itw, gp_lat);
    return interp(itw, surface_elevation.data(joker, 0), gp_lat[0]);
    
  } else if (atmosphere_dim == 3) {
    const Vector& lat_grid = surface_elevation.get_numeric_grid(GFIELD2_LAT_GRID);
    const Vector& lon_grid = surface_elevation.get_numeric_grid(GFIELD2_LON_GRID);
    ArrayOfGridPos gp_lat(1);
    if (lat_grid.nelem() > 1) {
      gridpos(gp_lat, lat_grid, pos[1], inf_proxy);
      jacobian_type_extrapol(gp_lat);
    } else {
      gp4length1grid(gp_lat);
    }
    ArrayOfGridPos gp_lon(1);
    if (lon_grid.nelem() > 1) {
      gridpos(gp_lon, lon_grid, pos[2], inf_proxy);
      jacobian_type_extrapol(gp_lon);
    } else {
      gp4length1grid(gp_lon);
    }
    Vector itw(4);
    interpweights(itw, gp_lat, gp_lon);
    return interp(itw, surface_elevation.data, gp_lat[0], gp_lon[0]);

  } else {
    ARTS_ASSERT(0);
  }
}
  
