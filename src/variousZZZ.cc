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
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "lin_alg.h"
#include "logic.h"
#include "messages.h"
#include "variousZZZ.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);


/* Workspace method: Doxygen documentation will be auto-generated */
void AltLatLonFieldSet(GriddedField3& gfield3,
                       const Vector& altitude_grid,
                       const Vector& latitude_grid,
                       const Vector& longitude_grid,
                       const Tensor3& data,
                       const String& name,
                       const Verbosity&) {
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
  ARTS_USER_ERROR_IF(data.npages() != altitude_grid.nelem(),
                     "Inconsistent altitude size!\n"
                     "Length of altitude grid: ", altitude_grid.nelem(), "\n"
                     "Altitude size of data: ", data.npages());
  ARTS_USER_ERROR_IF(data.nrows() != latitude_grid.nelem(),
                     "Inconsistent latitude size!\n"
                     "Length of latitude grid: ", latitude_grid.nelem(), "\n"
                     "Latitude size of data: ", data.nrows());
  ARTS_USER_ERROR_IF(data.ncols() != longitude_grid.nelem(),
                     "Inconsistent longitude size!\n"
                     "Length of longitude grid: ", longitude_grid.nelem(), "\n"
                     "Longitude size of data: ", data.ncols());

  gfield3.set_name(name);

  gfield3.set_grid_name(0, "Altitude");
  gfield3.set_grid(0, altitude_grid);
  gfield3.set_grid_name(1, "Latitude");
  gfield3.set_grid(1, latitude_grid);
  gfield3.set_grid_name(2, "Longitude");
  gfield3.set_grid(2, longitude_grid);

  gfield3.data = data;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AltLatLonFieldSetConstant(GriddedField3& gfield3,
                               const Numeric& value,
                               const String& name,
                               const Verbosity&) {
  gfield3.set_name(name);

  gfield3.set_grid_name(0, "Altitude");
  gfield3.set_grid(0, Vector(0, 0));
  gfield3.set_grid_name(1, "Latitude");
  gfield3.set_grid(1, Vector(1, 0));
  gfield3.set_grid_name(2, "Longitude");
  gfield3.set_grid(2, Vector(1, 0));

  gfield3.data.resize(1, 1, 1);
  gfield3.data(0, 0, 0) = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void LatLonFieldSet(GriddedField2& gfield2,
                    const Vector& latitude_grid,
                    const Vector& longitude_grid,
                    const Matrix& data,
                    const String& name,
                    const Verbosity&) {
  ARTS_USER_ERROR_IF(!is_increasing(latitude_grid),
                     "*latitude_grid* must be strictly increasing.");
  ARTS_USER_ERROR_IF(min(latitude_grid) < -90 || max(latitude_grid) > 90,
                     "All values in *latitude_grid* must be inside [-90,90].");
  ARTS_USER_ERROR_IF(!is_increasing(longitude_grid),
                     "*longitude_grid* must be strictly increasing.");
  ARTS_USER_ERROR_IF(min(longitude_grid) < -180 || max(longitude_grid) >= 180,
                     "All values in *latitude_grid* must be inside [-180,180[.");
  ARTS_USER_ERROR_IF(data.nrows() != latitude_grid.nelem(),
                     "Inconsistent latitude size!\n"
                     "Length of latitude grid: ", latitude_grid.nelem(), "\n"
                     "Latitude size of data: ", data.nrows());
  ARTS_USER_ERROR_IF(data.ncols() != longitude_grid.nelem(),
                     "Inconsistent longitude size!\n"
                     "Length of longitude grid: ", longitude_grid.nelem(), "\n"
                     "Longitude size of data: ", data.ncols());

  gfield2.set_name(name);

  gfield2.set_grid_name(0, "Latitude");
  gfield2.set_grid(0, latitude_grid);

  gfield2.set_grid_name(1, "Longitude");
  gfield2.set_grid(1, longitude_grid);

  gfield2.data = data;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void LatLonFieldSetConstant(GriddedField2& gfield2,
                            const Numeric& value,
                            const String& name,
                            const Verbosity&) {
  gfield2.set_name(name);

  gfield2.set_grid_name(0, "Latitude");
  gfield2.set_grid(0, Vector(1, 0));
  gfield2.set_grid_name(1, "Longitude");
  gfield2.set_grid(1, Vector(1, 0));

  gfield2.data.resize(1, 1);
  gfield2.data(0, 0) = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericInterpAltLatLonField(Numeric& value,
                                 const GriddedField3& gfield3,
                                 const Vector& pos,
                                 const Verbosity&)
{
  ARTS_USER_ERROR_IF(gfield3.get_grid_name(0) != "Altitude",
                     "Name of first grid must be \"Altitude\".");
  ARTS_USER_ERROR_IF(gfield3.get_grid_name(1) != "Latitude",
                     "Name of second grid must be \"Latitude\".");
  ARTS_USER_ERROR_IF(gfield3.get_grid_name(2) != "Longitude",
                     "Name of third grid must be \"Longitude\".");
  gfield3.checksize_strict();
  chk_rte_pos("pos", pos);
  
  value = interp_gfield3(gfield3, pos);
}                  


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericInterpLatLonField(Numeric& value,
                              const GriddedField2& gfield2,
                              const Vector& pos,
                              const Verbosity&)
{
  ARTS_USER_ERROR_IF(gfield2.get_grid_name(0) != "Latitude",
                     "Name of first grid must be \"Latitude\".");
  ARTS_USER_ERROR_IF(gfield2.get_grid_name(1) != "Longitude",
                     "Name of second grid must be \"Longitude\".");
  gfield2.checksize_strict();
  chk_rte_pos("pos", pos);
  
  value = interp_gfield2(gfield2, pos[Range(1, 2)]);
}                  





//
// So far a local solution 
//
void gridpos_local(ArrayOfGridPos& gp,
                 ConstVectorView grid, 
                 ConstVectorView points) {
  const Index n = points.nelem();
  ARTS_ASSERT(gp.nelem() == n);

  // To save time in case of grid length 1
  const Index l = grid.nelem();
  if (l == 1) {
    gp4length1grid(gp);
    return;
  }
    
  for (Index i=0; i<n; ++i) {
    Numeric x = points[i];
    
    // Handle nearest extrapolation
    if (x < grid[0])
      x = grid[0];
    else if (x > grid[l - 1])
      x = grid[l - 1];

    gridpos(gp[i], grid, x);
  }
}


Numeric interp_gfield2(const GriddedField2& G,
                       const Vector& pos2D)
{
  ARTS_ASSERT(pos2D.nelem() == 2);
  ARTS_ASSERT(G.checksize());

  // Sizes
  const Index ncols = G.data.ncols();
  const Index nrows = G.data.nrows();

  // To save time in case of data size (1, 1)
  if (nrows == 1 && ncols == 1) {
    return G.data(0, 0);
  }
  
  // Get grid positions
  ArrayOfGridPos gp_row(1), gp_col(1);
  gridpos_local(gp_row, G.get_numeric_grid(0), Vector(1, pos2D[0]));
  gridpos_local(gp_col, G.get_numeric_grid(1), Vector(1, pos2D[1]));

  // Interpolate
  Vector itw(4);
  interpweights(itw, gp_row[0], gp_col[0]);
  return interp(itw, G.data, gp_row[0], gp_col[0]);  
}


Numeric interp_gfield3(const GriddedField3& G,
                       const Vector& pos)
{
  ARTS_ASSERT(pos.nelem() == 3);
  ARTS_ASSERT(G.checksize());

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
  gridpos_local(gp_pages, G.get_numeric_grid(0), Vector(1, pos[0]));
  gridpos_local(gp_row, G.get_numeric_grid(1), Vector(1, pos[1]));
  gridpos_local(gp_col, G.get_numeric_grid(2), Vector(1, pos[2]));

  // Interpolate
  Vector itw(8);
  interpweights(itw, gp_pages[0], gp_row[0], gp_col[0]);
  return interp(itw, G.data, gp_pages[0], gp_row[0], gp_col[0]);  
}


void refr_index_and_its_gradients(Numeric& refr_index_air,
                                  Numeric& refr_index_air_group,
                                  Numeric& dndz,
                                  Numeric& dndlat,
                                  Numeric& dndlon,
                                  Workspace& ws,
                                  const Agenda& refr_index_air_agenda,
                                  ConstVectorView pos,
                                  const bool& do_horizontal_gradients,
                                  const bool& do_twosided_perturb)
{
  // n at pos itself
  refr_index_air_ZZZ_agendaExecute(ws,
                                   refr_index_air,
                                   refr_index_air_group,
                                   pos,
                                   refr_index_air_agenda);

  // Altitude gradient
  Numeric n, dummy;
  {
    const Numeric dz = 1;
    Vector pos_shifted = pos;
    pos_shifted[0] += dz;
    refr_index_air_ZZZ_agendaExecute(ws,
                                     n,
                                     dummy,
                                     pos_shifted,
                                     refr_index_air_agenda);
    if (do_twosided_perturb) {
      pos_shifted[0] = pos[0] - dz;
      Numeric n2;
      refr_index_air_ZZZ_agendaExecute(ws,
                                       n2,
                                       dummy,
                                       pos_shifted,
                                       refr_index_air_agenda);
      dndz = (n - n2) / (2 * dz);
    } else {
      dndz = (n - refr_index_air) / dz;
    }
  }

  // Latitude and longitide gradients
  if (do_horizontal_gradients) {
    // Latitude
    {
      const Numeric dlat = 1e-4;
      Vector pos_shifted = pos;
      pos_shifted[1] += dlat;
      if (pos_shifted[1] > 90)  // We simply cut at 90 deg, will underestimate 
        pos_shifted[1] = 90;    // gradient but should not be of any practical
      refr_index_air_ZZZ_agendaExecute(ws,    // consequence
                                       n,
                                       dummy,
                                       pos_shifted,
                                       refr_index_air_agenda);
      if (do_twosided_perturb) {
        pos_shifted[1] = pos[1] - dlat;
        if (pos_shifted[1] < -90)  // As above, but -90 deg
          pos_shifted[1] = -90;    
        Numeric n2;
        refr_index_air_ZZZ_agendaExecute(ws,
                                         n2,
                                         dummy,
                                         pos_shifted,
                                         refr_index_air_agenda);
        dndlat = (n - n2) / (2 * dlat);
      } else {
        dndlat = (n - refr_index_air) / dlat;
      }
    }
    // Longitude
    { 
      const Numeric dlon = 1e-4;
      Vector pos_shifted = pos;
      pos_shifted[2] += dlon;
      if (pos_shifted[2] >= 180)  // We can't go above 180 deg
        pos_shifted[2] -= 180;
      refr_index_air_ZZZ_agendaExecute(ws,
                                       n,
                                       dummy,
                                       pos_shifted,
                                       refr_index_air_agenda);
      if (do_twosided_perturb) {
        pos_shifted[2] = pos[2] - dlon;
        if (pos_shifted[2] < -180)  // We can't go below -180 deg
          pos_shifted[2] += 360;
        Numeric n2;
        refr_index_air_ZZZ_agendaExecute(ws,
                                         n2,
                                         dummy,
                                         pos_shifted,
                                         refr_index_air_agenda);
        dndlon = (n - n2) / (2 * dlon);
      } else {
        dndlon = (n - refr_index_air) / dlon;
      }
    }
  } else {
    dndlat = 0.0;
    dndlon = 0.0;
  }
}
