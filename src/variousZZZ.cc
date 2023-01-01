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
 * @date   2023-01-01
 *
 * @brief  Functions to later be placed elsewhere
 */

#include "arts_conversions.h"
#include "geodeticZZZ.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "messages.h"
#include "variousZZZ.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);

using GriddedFieldGrids::GFIELD2_LAT_GRID;
using GriddedFieldGrids::GFIELD2_LON_GRID;


/* Workspace method: Doxygen documentation will be auto-generated */
void surface_elevationSet(GriddedField2& surface_elevation,
                          const Vector& latitude_grid,
                          const Vector& longitude_grid,
                          const Matrix& elevations,          
                          const Verbosity&) {
  surface_elevation.set_name("Surface elevation map");

  surface_elevation.set_grid_name(GFIELD2_LAT_GRID, "Latitude");
  surface_elevation.set_grid(GFIELD2_LAT_GRID, latitude_grid);

  surface_elevation.set_grid_name(GFIELD2_LON_GRID, "Longitude");
  surface_elevation.set_grid(GFIELD2_LON_GRID, longitude_grid);

  surface_elevation.data = elevations;

  // To not be too restrictive, we check assuming 3D
  chk_surface_elevation(surface_elevation);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_elevationSetConstant(GriddedField2& surface_elevation,
                                  const Numeric& elevation,
                                  const Verbosity&) {
  surface_elevation.set_name("Surface elevation map");

  surface_elevation.set_grid_name(GFIELD2_LAT_GRID, "Latitude");
  surface_elevation.set_grid(GFIELD2_LAT_GRID, Vector(1, 0));

  surface_elevation.set_grid_name(GFIELD2_LON_GRID, "Longitude");
  surface_elevation.set_grid(GFIELD2_LON_GRID, Vector(1, 0));

  surface_elevation.data.resize(1,1);
  surface_elevation.data(0,0) = elevation;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SurfaceElevationInterp(Numeric& elevation,
                            const GriddedField2& surface_elevation,
                            const Vector& pos,
                            const Verbosity&)
{
  surface_elevation.checksize_strict();
  chk_rte_pos("pos", pos);
  
  elevation = surface_z_at_pos(pos, surface_elevation);
}                  



void chk_rte_los(const String& name,
                 ConstVectorView los)
{
  ARTS_USER_ERROR_IF(los.nelem() != 2,
                      "The vector *", name, "* must have length 2.")
  ARTS_USER_ERROR_IF(los[0] < -90 || los[0] > 90,
      "The zenith angle in *", name, "* must be in the range [-90,90].")
  ARTS_USER_ERROR_IF(los[1] < -180 || los[1] > 180,
      "The azimuth angle in *", name, "* must be in the range [-180,180].")
}


void chk_rte_pos(const String& name,
                 ConstVectorView pos)
{
  ARTS_USER_ERROR_IF(pos.nelem() != 3,
                      "The vector *", name, "* must have length 3.")
  ARTS_USER_ERROR_IF(pos[1] < -90 || pos[1] > 90,
      "The latitude in *", name, "* must be in the range [-90,90].")
  ARTS_USER_ERROR_IF(pos[2] < -180 || pos[2] > 360,
      "The longitude in *", name, "* must be in the range [-180,360].")
}


void chk_refellipsoidZZZ(ConstVectorView refellipsoid)
{
  ARTS_USER_ERROR_IF(refellipsoid.nelem() != 2,
                     "*refellipsoid* must have two elements.");
  ARTS_USER_ERROR_IF(refellipsoid[0] <= 0 || refellipsoid[1] <= 0,
                     "All elements of *refellipsoid* must be > 0.");
  ARTS_USER_ERROR_IF(abs(refellipsoid[1]/refellipsoid[0]-1) > 0.5,
      "The ratio of the two radii in *refellipsoid* is outisde of [0.5,1.5].\n"
      "Do you really want to have such a flat reference ellipsoid?");
}


void chk_sensor_pos(const String& name,
                    ConstMatrixView sensor_pos)
{
  ARTS_USER_ERROR_IF(sensor_pos.ncols() != 3,
                     "*", name, "* must have three columns.");
  ARTS_USER_ERROR_IF(sensor_pos.nrows() == 0,
                     "*", name, "*must have at least one row.");
  for (Index i=0; i<sensor_pos.nrows(); i++) {
    ARTS_USER_ERROR_IF(sensor_pos(i,1) < -90 || sensor_pos(i,1) > 90,
                       "Unvalid latitude in *", name, "*.\n",
                       "Latitudes must be inside the range [-90,90],\n",
                       "but ", name, "(", i, ",1) is ", sensor_pos(i,1));
    ARTS_USER_ERROR_IF(sensor_pos(i,2) < -180 || sensor_pos(i,1) > 360,
                       "Unvalid longitude in *", name, "*.\n",
                       "Longitudes must be inside the range [-180,360],\n",
                       "but ", name, "(", i, ",2) is ", sensor_pos(i,2));
  }
}


void chk_sensor_los(const String& name,
                    ConstMatrixView sensor_los)
{
  ARTS_USER_ERROR_IF(sensor_los.ncols() != 2,
                     "*", name, "* must have two columns.");
  ARTS_USER_ERROR_IF(sensor_los.nrows() == 0,
                     "*", name, "* must have at least one row.");
  for (Index i=0; i<sensor_los.nrows(); i++) {
    ARTS_USER_ERROR_IF(sensor_los(i,0) < 0 || sensor_los(i,0) > 180,
                       "Unvalid zenith angle in *", name, "*.\n",
                       "Zenith angles must be inside the range [0,180],\n"
                       "but ", name, "(", i, ",0) is ", sensor_los(i,0));
    ARTS_USER_ERROR_IF(sensor_los(i,1) < -180 || sensor_los(i,1) > 180,
                       "Unvalid azimuth angle in *", name, "*.\n",
                       "Azimuth angles must be inside the range [-180,180],\n"
                       "but ", name, "(", i, ",1) is ", sensor_los(i,1));
  }
}

void chk_sensor_poslos(const String& name1,
                       ConstMatrixView sensor_pos,
                       const String& name2,
                       ConstMatrixView sensor_los) {
  chk_sensor_pos(name1, sensor_pos);
  chk_sensor_los(name2, sensor_los);
  ARTS_USER_ERROR_IF(sensor_los.nrows() != sensor_pos.nrows(),
                     "*", name1, "* and *", name2,
                     "* must have the same number of rows.");
}


void chk_surface_elevation(const GriddedField2& surface_elevation) {
  const Vector& lat_grid = surface_elevation.get_numeric_grid(GFIELD2_LAT_GRID);
  ARTS_USER_ERROR_IF(surface_elevation.data.nrows() != lat_grid.nelem(),
                     "Inconsistent latitude size in *surface_elevation*\n"
                     "Length of latitude grid: ", lat_grid.nelem(), "\n"
                     "Latitude size of data: ", surface_elevation.data.nrows());
  const Vector& lon_grid = surface_elevation.get_numeric_grid(GFIELD2_LON_GRID);
  ARTS_USER_ERROR_IF(surface_elevation.data.ncols() != lon_grid.nelem(),
                     "Inconsistent longitude size in *surface_elevation*\n"
                     "Length of longitude grid: ", lon_grid.nelem(), "\n"
                     "Longitude size of data: ", surface_elevation.data.ncols());
  ARTS_USER_ERROR_IF(surface_elevation.data.empty(),
                     "The data in *surface_elevation* are empty. Not allowed!");
}


Numeric find_crossing_with_surface_z(const Vector rte_pos,
                                     const Vector rte_los,
                                     const Vector ecef,
                                     const Vector decef,
                                     const Vector& refellipsoid,
                                     const GriddedField2& surface_elevation,
                                     const Numeric& surface_search_accuracy,
                                     const Index& surface_search_safe)
{
  // Find min and max surface altitude
  const Numeric z_min = min(surface_elevation.data);
  const Numeric z_max = max(surface_elevation.data);

  // Catch upward looking cases that can not have a surface intersection
  if (rte_pos[0] >= z_max && rte_los[0] <= 90) {
    return -1;
  }

  // Check that observation position is above ground
  if (rte_pos[0] < z_max) {
    Numeric z_surf = surface_z_at_pos(rte_pos, surface_elevation);
    if (rte_pos[0] < z_surf - surface_search_accuracy)
      ARTS_USER_ERROR(
          "The sensor is below the surface. Not allowed!\n"
          "The sensor altitude is at ", rte_pos[0], " m\n"
          "The surface altitude is ", z_surf, " m\n"
          "The position is (lat,lon): (", rte_pos[1], ",", rte_pos[2], ")");
  }

  // Constant surface altitude (in comparison to *surface_search_accuracy*)
  if (z_max - z_min < surface_search_accuracy / 100) {
    // Catch cases with position on the ground, as they can fail if
    // intersection_altitude is used
    if (rte_pos[0] <= z_max) {
      return 0.0;
    } else {
      return intersection_altitude(ecef, decef, refellipsoid, z_min);
    }

    // The general case
  } else {
    // Find a distance that is guaranteed above or at surface
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
    if (surface_search_safe) {
      Numeric l_test =
          l_min - surface_search_accuracy / 2;  // Remove l/2 to get exact result
      bool above_surface = true;   // if true l_test is 0
      while (above_surface && l_test < l_max) {
        l_test += surface_search_accuracy;
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_test);
        Numeric z_surf = surface_z_at_pos(pos, surface_elevation);
        if (pos[0] < z_surf) above_surface = false;
      }
      if (above_surface) {
        return -1;
      } else {
        return l_test - surface_search_accuracy / 2;
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
        Numeric z_surf = surface_z_at_pos(pos, surface_elevation);
        if (pos[0] > z_surf) return -1;
      }
      // Start bisection
      while (l_max - l_min > 2 * surface_search_accuracy) {
        const Numeric l_test = (l_min + l_max) / 2;
        Vector pos(3);
        pos_at_distance(pos, ecef, decef, refellipsoid, l_test);
        Numeric z_surf =
            surface_z_at_pos(pos, surface_elevation);
        if (pos[0] >= z_surf)
          l_min = l_test;
        else
          l_max = l_test;
      }
      return (l_min + l_max) / 2;
    }
  }
}


Numeric surface_z_at_pos(const Vector& pos,
                         const GriddedField2& surface_elevation)
{
  // Lat and lon grids
  const Vector& lat_grid = surface_elevation.get_numeric_grid(GFIELD2_LAT_GRID);
  const Vector& lon_grid = surface_elevation.get_numeric_grid(GFIELD2_LON_GRID);
  const Index& nlat = lat_grid.nelem();
  const Index& nlon = lon_grid.nelem();
  //
  ARTS_ASSERT(surface_elevation.data.nrows() == nlat);
  ARTS_ASSERT(surface_elevation.data.ncols() == nlon);

  // There should be special functions for this interpolation
  
  // Handle the case of data size (1, 1)
  if (nlat == 1 && nlon == 1) {
    return surface_elevation.data(0, 0);
  }
  
  // Extract latitude, handling extrapolation on the same time
  Numeric lat = pos[1];
  if (lat < lat_grid[0])
    lat = lat_grid[0];
  else if (lat > lat_grid[nlat - 1])
    lat = lat_grid[nlat - 1];

  // Data size is (nlat, 1)
  if (nlon == 1) {
    ArrayOfGridPos gp_lat(1);
    gridpos(gp_lat, lat_grid, lat);
    Vector itw(2);
    interpweights(itw, gp_lat[0]);
    return interp(itw, surface_elevation.data(joker, 0), gp_lat[0]);
  }
  
  // Extract longitude, handling extrapolation and [-180,180] vs [0,360]
  Numeric lon = move_lon_to_range(pos[2], lon_grid[0], lon_grid[nlon - 1]);
  if (lon < lon_grid[0])
    lon = lon_grid[0];
  else if (lon > lon_grid[nlon])
    lon = lon_grid[nlon - 1];

  // Data size is (1, nlon)
  if (nlat == 1) {
    ArrayOfGridPos gp_lon(1);
    gridpos(gp_lon, lon_grid, lon);
    Vector itw(2);
    interpweights(itw, gp_lon[0]);
    return interp(itw, surface_elevation.data(0, joker), gp_lon[0]);
    
  // Data size is (nlat, nlon)
  } else {
    ArrayOfGridPos gp_lat(1);
    gridpos(gp_lat, lat_grid, lat);
    ArrayOfGridPos gp_lon(1);
    gridpos(gp_lon, lon_grid, lon);
    Vector itw(4);
    interpweights(itw, gp_lat[0], gp_lon[0]);
    return interp(itw, surface_elevation.data, gp_lat[0], gp_lon[0]);
  }
}
