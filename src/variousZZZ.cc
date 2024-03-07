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

#include "variousZZZ.h"

#include <workspace.h>

#include "arts_conversions.h"
#include "check_input.h"
#include "interpolation.h"
#include "logic.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);


/* Workspace method: Doxygen documentation will be auto-generated */
void AltLatLonFieldSet(GriddedField3& gfield3,
                       const Vector& altitude_grid,
                       const Vector& latitude_grid,
                       const Vector& longitude_grid,
                       const Tensor3& data,
                       const String& name) {
  ARTS_USER_ERROR_IF(!is_increasing(altitude_grid),
                     "*altitude_grid* must be strictly increasing.");
  ARTS_USER_ERROR_IF(!is_increasing(latitude_grid),
                     "*latitude_grid* must be strictly increasing.");
  ARTS_USER_ERROR_IF(min(latitude_grid) < -90 || max(latitude_grid) > 90,
                     "All values in *latitude_grid* must be inside [-90,90].");
  ARTS_USER_ERROR_IF(!is_increasing(longitude_grid),
                     "*longitude_grid* must be strictly increasing.");
  ARTS_USER_ERROR_IF(min(longitude_grid) < -180 || max(longitude_grid) >= 180,
                     "All values in *latitude_grid* must be inside [-180,180[.");
  ARTS_USER_ERROR_IF(data.npages() != altitude_grid.size(),
                     "Inconsistent altitude size!\n"
                     "Length of altitude grid: ", altitude_grid.size(), "\n"
                     "Altitude size of data: ", data.npages());
  ARTS_USER_ERROR_IF(data.nrows() != latitude_grid.size(),
                     "Inconsistent latitude size!\n"
                     "Length of latitude grid: ", latitude_grid.size(), "\n"
                     "Latitude size of data: ", data.nrows());
  ARTS_USER_ERROR_IF(data.ncols() != longitude_grid.size(),
                     "Inconsistent longitude size!\n"
                     "Length of longitude grid: ", longitude_grid.size(), "\n"
                     "Longitude size of data: ", data.ncols());

  gfield3.data_name = name;
  gfield3.grid_names = {"Altitude", "Latitude", "Longitude"};
  gfield3.grids = {altitude_grid, latitude_grid, longitude_grid};

  gfield3.data = data;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AltLatLonFieldSetConstant(GriddedField3& gfield3,
                               const Numeric& value,
                               const String& name) {
  gfield3.data_name = (name);
  gfield3.grid_names = {"Altitude", "Latitude", "Longitude"};
  gfield3.grids = {Vector(1, 0), Vector(1, 0), Vector(1, 0)};

  gfield3.data.resize(1, 1, 1);
  gfield3.data(0, 0, 0) = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void LatLonFieldSet(GriddedField2& gfield2,
                    const Vector& latitude_grid,
                    const Vector& longitude_grid,
                    const Matrix& data,
                    const String& name) {
  ARTS_USER_ERROR_IF(!is_increasing(latitude_grid),
                     "*latitude_grid* must be strictly increasing.");
  ARTS_USER_ERROR_IF(min(latitude_grid) < -90 || max(latitude_grid) > 90,
                     "All values in *latitude_grid* must be inside [-90,90].");
  ARTS_USER_ERROR_IF(!is_increasing(longitude_grid),
                     "*longitude_grid* must be strictly increasing.");
  ARTS_USER_ERROR_IF(min(longitude_grid) < -180 || max(longitude_grid) >= 180,
                     "All values in *latitude_grid* must be inside [-180,180[.");
  ARTS_USER_ERROR_IF(data.nrows() != latitude_grid.size(),
                     "Inconsistent latitude size!\n"
                     "Length of latitude grid: ", latitude_grid.size(), "\n"
                     "Latitude size of data: ", data.nrows());
  ARTS_USER_ERROR_IF(data.ncols() != longitude_grid.size(),
                     "Inconsistent longitude size!\n"
                     "Length of longitude grid: ", longitude_grid.size(), "\n"
                     "Longitude size of data: ", data.ncols());

  gfield2.data_name = name;
  gfield2.grid_names = {"Latitude", "Longitude"};
  gfield2.grids = {latitude_grid, longitude_grid};

  gfield2.data = data;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void LatLonFieldSetConstant(GriddedField2& gfield2,
                            const Numeric& value,
                            const String& name) {
  gfield2.data_name = name;
  gfield2.grid_names = {"Latitude", "Longitude"};
  gfield2.grids = {Vector(1, 0), Vector(1, 0)};

  gfield2.data.resize(1, 1);
  gfield2.data(0, 0) = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericInterpAltLatLonField(Numeric& value,
                                 const GriddedField3& gfield3,
                                 const Vector& pos)
{
  ARTS_USER_ERROR_IF(gfield3.gridname<0>() != "Altitude",
                     "Name of first grid must be \"Altitude\".");
  ARTS_USER_ERROR_IF(gfield3.gridname<1>() != "Latitude",
                     "Name of second grid must be \"Latitude\".");
  ARTS_USER_ERROR_IF(gfield3.gridname<2>() != "Longitude",
                     "Name of third grid must be \"Longitude\".");
  ARTS_USER_ERROR_IF(not gfield3.ok(), "Inconsistent grid sizes for:\n", gfield3);
  chk_rte_pos("pos", pos);
  
  value = interp_gfield3(gfield3, pos);
}                  


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericInterpLatLonField(Numeric& value,
                              const GriddedField2& gfield2,
                              const Vector& pos)
{
  ARTS_USER_ERROR_IF(gfield2.gridname<0>() != "Latitude",
                     "Name of first grid must be \"Latitude\".");
  ARTS_USER_ERROR_IF(gfield2.gridname<1>() != "Longitude",
                     "Name of second grid must be \"Longitude\".");
  ARTS_USER_ERROR_IF(not gfield2.ok(), "Inconsistent grid sizes for:\n", gfield2);
  chk_rte_pos("pos", pos);
  
  value = interp_gfield2(gfield2, Vector{pos[Range(1, 2)]});
}                  


//
// So far a local solution 
//
void gridpos_local(ArrayOfGridPos& gp,
                 ConstVectorView grid, 
                 ConstVectorView points) {
  const Size n = points.size();
  ARTS_ASSERT(gp.size() == n);

  // To save time in case of grid length 1
  const Index l = grid.size();
  if (l == 1) {
    gp4length1grid(gp);
    return;
  }
    
  for (Size i=0; i<n; ++i) {
    Numeric x = points[i];
    
    // Handle nearest extrapolation
    if (x < grid[0])
      x = grid[0];
    else if (x > grid[l - 1])
      x = grid[l - 1];

    gridpos(gp[i], grid, x);
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericInterpVector(Numeric& value,
                         const Vector& x,
                         const Vector& y,
                         const Numeric& xv)
{
  ArrayOfGridPos gp(1);
  gridpos_local(gp, x, Vector(1, xv));
  Vector itw(2);
  interpweights(itw, gp[0]);
  value = interp(itw, y, gp[0]);
}


Numeric interp_gfield2(const GriddedField2& G,
                       const Vector& pos2D)
{
  ARTS_ASSERT(pos2D.size() == 2);
  ARTS_ASSERT(G.ok());

  // Sizes
  const Index ncols = G.data.ncols();
  const Index nrows = G.data.nrows();

  // To save time in case of data size (1, 1)
  if (nrows == 1 && ncols == 1) {
    return G.data(0, 0);
  }
  
  // Get grid positions
  ArrayOfGridPos gp_row(1), gp_col(1);
  gridpos_local(gp_row, G.grid<0>(), Vector(1, pos2D[0]));
  gridpos_local(gp_col, G.grid<1>(), Vector(1, pos2D[1]));

  // Interpolate
  Vector itw(4);
  interpweights(itw, gp_row[0], gp_col[0]);
  return interp(itw, G.data, gp_row[0], gp_col[0]);  
}


Numeric interp_gfield3(const GriddedField3& G,
                       const Vector& pos)
{
  ARTS_ASSERT(pos.size() == 3);
  ARTS_ASSERT(G.ok());

  // Sizes
  const Index npages = G.data.npages();
  const Index ncols = G.data.ncols();
  const Index nrows = G.data.nrows();

  // To save time in case of data size (1, 1)
  if (npages == 1 && nrows == 1 && ncols == 1) {
    return G.data(0, 0, 0);
  }

  // Get grid positions
  ArrayOfGridPos gp_pages(1), gp_row(1), gp_col(1);
  gridpos_local(gp_pages, G.grid<0>(), Vector(1, pos[0]));
  gridpos_local(gp_row, G.grid<1>(), Vector(1, pos[1]));
  gridpos_local(gp_col, G.grid<2>(), Vector(1, pos[2]));

  // Interpolate
  Vector itw(8);
  interpweights(itw, gp_pages[0], gp_row[0], gp_col[0]);
  return interp(itw, G.data, gp_pages[0], gp_row[0], gp_col[0]);  
}
