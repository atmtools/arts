/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
  \file   m_atmosphere.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-16

  \brief  Workspace functions to set variables defining the atmosphere 
          (excluding the surface).

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cfloat>
#include <cmath>
#include <vector>
#include "arts_constants.h"
#include "arts_conversions.h"
#include "debug.h"
#include "species_tags.h"
#include "absorption.h"
#include "agenda_class.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "geodetic_OLD.h"
#include "global_data.h"
#include "gridded_fields.h"
#include "igrf13.h"
#include "interpolation.h"
#include "interp.h"
#include "linescaling.h"
#include "matpack_data.h"
#include "messages.h"
#include "rte.h"
#include "special_interp.h"
#include "xml_io.h"

using GriddedFieldGrids::GFIELD3_P_GRID;
using GriddedFieldGrids::GFIELD3_LAT_GRID;
using GriddedFieldGrids::GFIELD3_LON_GRID;
using GriddedFieldGrids::GFIELD4_FIELD_NAMES;
using GriddedFieldGrids::GFIELD4_P_GRID;
using GriddedFieldGrids::GFIELD4_LAT_GRID;
using GriddedFieldGrids::GFIELD4_LON_GRID;


inline constexpr Numeric GAS_CONSTANT=Constant::ideal_gas_constant;

//! Data value accuracy requirement for values at 0 and 360 deg if longitudes are cyclic
/*!
 */
inline constexpr Numeric EPSILON_LON_CYCLIC = 2 * DBL_EPSILON;

/*===========================================================================
 *=== Helper functions
 *===========================================================================*/

//! atm_fields_compactExpand
/*!
   Add a species to an *atm_fields_compact*. Does not add any content, but only
   resizes the data and adds a field to the *ArrayOfString* respresenting the
   species for this *GriddedField4*. This helper function is used by e.g.
   *atm_fields_compactAddSpecies* and *atm_fields_compactAddConstant*.
  

   \retval  af         The new atm_fields_compact 
   \retval  nf         The new number of fields
   \param   name       Name of new field
   \param   prepend    0 = append, otherwise prepend

   \author Gerrit Holl
   \date   2011-05-04
*/

void atm_fields_compactExpand(GriddedField4& af,
                              Index& nf,
                              const String& name,
                              const Index& prepend,
                              const Verbosity&) {
  // Number of fields already present:
  nf = af.get_string_grid(GFIELD4_FIELD_NAMES).nelem();

  if (prepend) {
    // Add name of new field to field name list:
    ArrayOfString& as = af.get_string_grid(GFIELD4_FIELD_NAMES);
    as.insert(as.begin(), name);
  } else {
    // Add name of new field to field name list:
    af.get_string_grid(GFIELD4_FIELD_NAMES).push_back(name);
  }

  // Save original fields:
  const Tensor4 dummy = af.data;

  // Adjust size:
  af.resize(nf + 1, dummy.npages(), dummy.nrows(), dummy.ncols());
  nf++;  // set to new number of fields

  // Copy back original field:
  af.data(Range((prepend && nf > 1) ? 1 : 0, nf - 1), joker, joker, joker) =
      dummy;
}

//! Calculate grid positions and interpolations weights for AtmFieldPRegrid
/*
 This helper function is used by all AtmFieldPRegrid WSMs to do the common
 calculation of the grid positions and interpolation weights for pressures.
 It is an adaptation of GriddedFieldPRegridHelper for Tensors instead of
 GriddedFields.

 \param[out]    ing_min       Index in the new grid with first value covered
                              by the old grid.
 \param[out]    ing_max       Index in the new grid with last value covered
                              by the old grid.
 \param[out]    gp_p          Pressure grid positions
 \param[out]    itw           Interpolation weights
 \param[in]     p_grid_out    New pressure grid
 \param[in]     p_grid_in     Old pressure grid
 \param[in]     interp_order  Interpolation order
 \param[in]     verbosity     Verbosity levels
 */
void AtmFieldPRegridHelper(Index& ing_min,
                           Index& ing_max,
                           ArrayOfLagrangeLogInterpolation& lag_p,
                           Matrix& itw,
                           ConstVectorView p_grid_out,
                           ConstVectorView p_grid_in,
                           const Index& interp_order,
                           const Verbosity& verbosity) {
  CREATE_OUT2;

  out2 << "  Interpolation order: " << interp_order << "\n";

  ing_min = 0;
  ing_max = p_grid_out.nelem() - 1;
  chk_interpolation_pgrids(
      "Atmospheric field to p_grid_out", p_grid_in, p_grid_out, interp_order);

  Index nelem_in_range = ing_max - ing_min + 1;

  // Calculate grid positions:
  if (nelem_in_range > 0) {
    lag_p = my_interp::lagrange_interpolation_list<LagrangeLogInterpolation>(
      p_grid_out, p_grid_in, interp_order, 0.5);
    itw = interpweights(lag_p);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldPRegrid(  // WS Generic Output:
    Tensor3& atmtensor_out,
    // WS Generic Input:
    const Tensor3& atmtensor_in_orig,
    const Vector& p_grid_new,
    const Vector& p_grid_old,
    const Index& interp_order,
    const Verbosity& verbosity) {
  // check that p_grid_old is the p_grid associated with atmtensor_in_orig. we can
  // only check for consistent size with p_grid_old.
  ARTS_USER_ERROR_IF (atmtensor_in_orig.npages() != p_grid_old.nelem(),
      "p_grid_old is supposed to be the p_grid associated with the "
      "atmospheric field.\n"
      "However, it is not as their sizes are inconsistent.\n")

  const Tensor3* atmtensor_in_pnt;
  Tensor3 atmtensor_in_copy;

  if (&atmtensor_in_orig == &atmtensor_out) {
    atmtensor_in_copy = atmtensor_in_orig;
    atmtensor_in_pnt = &atmtensor_in_copy;
  } else
    atmtensor_in_pnt = &atmtensor_in_orig;

  const Tensor3& atmtensor_in = *atmtensor_in_pnt;

  // Resize output tensor
  atmtensor_out.resize(
      p_grid_new.nelem(), atmtensor_in.nrows(), atmtensor_in.ncols());

  ArrayOfLagrangeLogInterpolation lag_p;
  Matrix itw;  // nb. it is invalid to use this as it stands here...

  Index ing_min, ing_max;

  AtmFieldPRegridHelper(ing_min,
                        ing_max,
                        lag_p,
                        itw,
                        p_grid_new,
                        p_grid_old,
                        interp_order,
                        verbosity);

  // Interpolate:
  ARTS_USER_ERROR_IF ((ing_max - ing_min < 0) ||
      (ing_max - ing_min + 1 != p_grid_new.nelem()),
      "New grid seems not to be sufficiently covered by old grid.\n")

  for (Index i = 0; i < atmtensor_in.nrows(); i++)
    for (Index j = 0; j < atmtensor_in.ncols(); j++)
      reinterp(atmtensor_out(joker, i, j), atmtensor_in(joker, i, j), itw, lag_p);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldPRegrid(  // WS Generic Output:
    Tensor4& atmtensor_out,
    // WS Generic Input:
    const Tensor4& atmtensor_in_orig,
    const Vector& p_grid_new,
    const Vector& p_grid_old,
    const Index& interp_order,
    const Verbosity& verbosity) {
  const Tensor4* atmtensor_in_pnt;
  Tensor4 atmtensor_in_copy;

  if (&atmtensor_in_orig == &atmtensor_out) {
    atmtensor_in_copy = atmtensor_in_orig;
    atmtensor_in_pnt = &atmtensor_in_copy;
  } else
    atmtensor_in_pnt = &atmtensor_in_orig;

  const Tensor4& atmtensor_in = *atmtensor_in_pnt;

  // Resize output tensor
  atmtensor_out.resize(atmtensor_in.nbooks(),
                       p_grid_new.nelem(),
                       atmtensor_in.nrows(),
                       atmtensor_in.ncols());
  
  ArrayOfLagrangeLogInterpolation lag_p;
  Matrix itw;  // nb. it is invalid to use this as it stands here...

  Index ing_min, ing_max;

  AtmFieldPRegridHelper(ing_min,
                        ing_max,
                        lag_p,
                        itw,
                        p_grid_new,
                        p_grid_old,
                        interp_order,
                        verbosity);

  // Interpolate:
  ARTS_USER_ERROR_IF ((ing_max - ing_min < 0) ||
                      (ing_max - ing_min + 1 != p_grid_new.nelem()),
    "New grid seems not to be sufficiently covered by old grid.\n")

  for (Index b = 0; b < atmtensor_in.nbooks(); b++)
    for (Index i = 0; i < atmtensor_in.nrows(); i++)
      for (Index j = 0; j < atmtensor_in.ncols(); j++)
        reinterp(atmtensor_out(b, joker, i, j),
                 atmtensor_in(b, joker, i, j),
                 itw,
                 lag_p);
}

//! Check for correct grid dimensions
/*
 Helper function for FieldFromGriddedField functions to ensure correct
 dimension and values of latitude/longitude grids

 \param[in]  lat_grid   Latitude grid
 \param[in]  lon_grid   Longitude grid
 \param[in]  ilat       Latitude grid index in gfield
 \param[in]  ilon       Longitude grid index in gfield
 \param[in]  gfield     GriddedField
 */
void FieldFromGriddedFieldCheckLatLonHelper(const Vector& lat_grid,
                                            const Vector& lon_grid,
                                            const Index ilat,
                                            const Index ilon,
                                            const GriddedField& gfield) {
  chk_griddedfield_gridname(gfield, ilat, "Latitude");
  chk_griddedfield_gridname(gfield, ilon, "Longitude");

  if (lon_grid.empty()) {
    chk_size("gfield.lon_grid", gfield.get_numeric_grid(ilon), 1);

    if (lat_grid.empty())
      chk_size("gfield.lat_grid", gfield.get_numeric_grid(ilat), 1);
    else
      chk_if_equal("lat_grid",
                   "gfield.lat_grid",
                   lat_grid,
                   gfield.get_numeric_grid(ilat));
  } else {
    chk_if_equal(
        "lat_grid", "gfield.lat_grid", lat_grid, gfield.get_numeric_grid(ilat));
    chk_if_equal(
        "lon_grid", "gfield.lon_grid", lon_grid, gfield.get_numeric_grid(ilon));
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(  // WS Generic Output:
    Matrix& field_out,
    // WS Input:
    const Vector& p_grid _U_,
    const Vector& lat_grid,
    const Vector& lon_grid,
    // WS Generic Input:
    const GriddedField2& gfraw_in,
    const Verbosity&) {
  FieldFromGriddedFieldCheckLatLonHelper(lat_grid, lon_grid, 0, 1, gfraw_in);

  field_out = gfraw_in.data;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(  // WS Generic Output:
    Tensor3& field_out,
    // WS Input:
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    // WS Generic Input:
    const GriddedField3& gfraw_in,
    const Verbosity&) {
  chk_griddedfield_gridname(gfraw_in, 0, "Pressure");
  chk_if_equal("p_grid", "gfield.p_grid", p_grid, gfraw_in.get_numeric_grid(0));

  FieldFromGriddedFieldCheckLatLonHelper(lat_grid, lon_grid, 1, 2, gfraw_in);

  field_out = gfraw_in.data;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(  // WS Generic Output:
    Tensor4& field_out,
    // WS Input:
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    // WS Generic Input:
    const GriddedField4& gfraw_in,
    const Verbosity&) {
  chk_griddedfield_gridname(gfraw_in, 1, "Pressure");
  chk_if_equal("p_grid", "gfield.p_grid", p_grid, gfraw_in.get_numeric_grid(0));

  FieldFromGriddedFieldCheckLatLonHelper(lat_grid, lon_grid, 2, 3, gfraw_in);

  field_out = gfraw_in.data;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(  // WS Generic Output:
    Tensor4& field_out,
    // WS Input:
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    // WS Generic Input:
    const ArrayOfGriddedField3& gfraw_in,
    const Verbosity& verbosity) {
  if (!gfraw_in.nelem()) {
    CREATE_OUT1;
    out1 << "   Warning: gfraw_in is empty, proceeding anyway\n";
    field_out.resize(0, 0, 0, 0);
  } else {
    field_out.resize(gfraw_in.nelem(),
                     p_grid.nelem(),
                     gfraw_in[0].data.nrows(),
                     gfraw_in[0].data.ncols());
  }

  for (Index i = 0; i < gfraw_in.nelem(); i++) {
    chk_griddedfield_gridname(gfraw_in[i], 0, "Pressure");
    chk_if_equal(
        "p_grid", "gfield.p_grid", p_grid, gfraw_in[i].get_numeric_grid(0));

    FieldFromGriddedFieldCheckLatLonHelper(
        lat_grid, lon_grid, 1, 2, gfraw_in[i]);

    field_out(i, joker, joker, joker) = gfraw_in[i].data;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(  // WS Generic Output:
    GriddedField2& gfraw_out,
    // WS Generic Input:
    const GriddedField2& gfraw_in_orig,
    const Verbosity&) {
  const GriddedField2* gfraw_in_pnt;
  GriddedField2 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField2& gfraw_in = *gfraw_in_pnt;

  chk_griddedfield_gridname(gfraw_in, 0, "Latitude");
  chk_griddedfield_gridname(gfraw_in, 1, "Longitude");

  ARTS_USER_ERROR_IF (gfraw_in.data.ncols() != 1 && gfraw_in.data.nrows() != 1,
        "Can't expand data because number of Latitudes and Longitudes is greater than 1");

  gfraw_out.set_grid_name(0, "Latitude");
  gfraw_out.set_grid_name(1, "Longitude");

  Vector v(2);
  if (gfraw_in.data.nrows() == 1 && gfraw_in.data.ncols() != 1) {
    v[0] = -90;
    v[1] = 90;
    gfraw_out.set_grid(0, v);
    gfraw_out.resize(2, gfraw_in.data.ncols());

    for (Index j = 0; j < gfraw_in.data.ncols(); j++)
      gfraw_out.data(joker, j) = gfraw_in.data(0, j);
  } else if (gfraw_in.data.nrows() != 1 && gfraw_in.data.ncols() == 1) {
    v[0] = 0;
    v[1] = 360;
    gfraw_out.set_grid(1, v);
    gfraw_out.resize(gfraw_in.data.nrows(), 2);

    for (Index j = 0; j < gfraw_in.data.nrows(); j++)
      gfraw_out.data(j, joker) = gfraw_in.data(j, 0);
  } else {
    v[0] = -90;
    v[1] = 90;
    gfraw_out.set_grid(0, v);
    v[0] = 0;
    v[1] = 360;
    gfraw_out.set_grid(1, v);
    gfraw_out.resize(2, 2);

    gfraw_out.data = gfraw_in.data(0, 0);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(  // WS Generic Output:
    GriddedField3& gfraw_out,
    // WS Generic Input:
    const GriddedField3& gfraw_in_orig,
    const Verbosity&) {
  const GriddedField3* gfraw_in_pnt;
  GriddedField3 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField3& gfraw_in = *gfraw_in_pnt;

  chk_griddedfield_gridname(gfraw_in, 1, "Latitude");
  chk_griddedfield_gridname(gfraw_in, 2, "Longitude");

  ARTS_USER_ERROR_IF (gfraw_in.data.ncols() != 1 && gfraw_in.data.nrows() != 1,
        "Can't expand data because number of Latitudes and Longitudes is greater than 1");

  gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
  gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));
  gfraw_out.set_grid_name(1, "Latitude");
  gfraw_out.set_grid_name(2, "Longitude");

  Vector v(2);
  if (gfraw_in.data.nrows() == 1 && gfraw_in.data.ncols() != 1) {
    v[0] = -90;
    v[1] = 90;
    gfraw_out.set_grid(1, v);
    gfraw_out.resize(gfraw_in.data.npages(), 2, gfraw_in.data.ncols());

    for (Index i = 0; i < gfraw_in.data.npages(); i++)
      for (Index j = 0; j < gfraw_in.data.ncols(); j++)
        gfraw_out.data(i, joker, j) = gfraw_in.data(i, 0, j);
  } else if (gfraw_in.data.nrows() != 1 && gfraw_in.data.ncols() == 1) {
    v[0] = 0;
    v[1] = 360;
    gfraw_out.set_grid(2, v);
    gfraw_out.resize(gfraw_in.data.npages(), gfraw_in.data.nrows(), 2);

    for (Index i = 0; i < gfraw_in.data.npages(); i++)
      for (Index j = 0; j < gfraw_in.data.nrows(); j++)
        gfraw_out.data(i, j, joker) = gfraw_in.data(i, j, 0);
  } else {
    v[0] = -90;
    v[1] = 90;
    gfraw_out.set_grid(1, v);
    v[0] = 0;
    v[1] = 360;
    gfraw_out.set_grid(2, v);
    gfraw_out.resize(gfraw_in.data.npages(), 2, 2);

    for (Index i = 0; i < gfraw_in.data.npages(); i++)
      gfraw_out.data(i, joker, joker) = gfraw_in.data(i, 0, 0);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(  // WS Generic Output:
    GriddedField4& gfraw_out,
    // WS Generic Input:
    const GriddedField4& gfraw_in_orig,
    const Verbosity&) {
  const GriddedField4* gfraw_in_pnt;
  GriddedField4 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField4& gfraw_in = *gfraw_in_pnt;

  chk_griddedfield_gridname(gfraw_in, 2, "Latitude");
  chk_griddedfield_gridname(gfraw_in, 3, "Longitude");

  ARTS_USER_ERROR_IF (gfraw_in.data.ncols() != 1 && gfraw_in.data.nrows() != 1,
        "Can't expand data because number of Latitudes and Longitudes is greater than 1");

  gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
  gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));

  gfraw_out.set_grid(1, gfraw_in.get_numeric_grid(1));
  gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));

  gfraw_out.set_grid_name(2, "Latitude");
  gfraw_out.set_grid_name(3, "Longitude");

  Vector v(2);
  if (gfraw_in.data.nrows() == 1 && gfraw_in.data.ncols() != 1) {
    v[0] = -90;
    v[1] = 90;
    gfraw_out.set_grid(2, v);
    gfraw_out.resize(gfraw_in.data.nbooks(),
                     gfraw_in.data.npages(),
                     2,
                     gfraw_in.data.ncols());

    for (Index k = 0; k < gfraw_in.data.nbooks(); k++)
      for (Index i = 0; i < gfraw_in.data.npages(); i++)
        for (Index j = 0; j < gfraw_in.data.ncols(); j++)
          gfraw_out.data(k, i, joker, j) = gfraw_in.data(k, i, 0, j);
  } else if (gfraw_in.data.nrows() != 1 && gfraw_in.data.ncols() == 1) {
    v[0] = 0;
    v[1] = 360;
    gfraw_out.set_grid(3, v);
    gfraw_out.resize(gfraw_in.data.nbooks(),
                     gfraw_in.data.npages(),
                     gfraw_in.data.nrows(),
                     2);

    for (Index k = 0; k < gfraw_in.data.nbooks(); k++)
      for (Index i = 0; i < gfraw_in.data.npages(); i++)
        for (Index j = 0; j < gfraw_in.data.nrows(); j++)
          gfraw_out.data(k, i, j, joker) = gfraw_in.data(k, i, j, 0);
  } else {
    v[0] = -90;
    v[1] = 90;
    gfraw_out.set_grid(2, v);
    v[0] = 0;
    v[1] = 360;
    gfraw_out.set_grid(3, v);
    gfraw_out.resize(gfraw_in.data.nbooks(), gfraw_in.data.npages(), 2, 2);

    for (Index k = 0; k < gfraw_in.data.nbooks(); k++)
      for (Index i = 0; i < gfraw_in.data.npages(); i++)
        gfraw_out.data(k, i, joker, joker) = gfraw_in.data(k, i, 0, 0);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(  // WS Generic Output:
    ArrayOfGriddedField3& gfraw_out,
    // WS Generic Input:
    const ArrayOfGriddedField3& gfraw_in,
    const Verbosity& verbosity) {
  gfraw_out.resize(gfraw_in.nelem());

  for (Index i = 0; i < gfraw_in.nelem(); i++)
    GriddedFieldLatLonExpand(gfraw_out[i], gfraw_in[i], verbosity);
}

//! Calculate grid positions and interpolations weights for GriddedFieldPRegrid
/*
 This helper function is used by all GriddedFieldPRegrid WSMs to do the common
 calculation of the grid positions and interpolation weights for pressures.

 \param[out]    ing_min       Index in the new grid with first value covered
                              by the old grid.
 \param[out]    ing_max       Index in the new grid with last value covered
                              by the old grid.
 \param[out]    gp_p          Pressure grid positions
 \param[out]    itw           Interpolation weights
 \param[in,out] gfraw_out     Output GriddedField
 \param[in]     gfraw_in      Input GriddedField
 \param[in]     p_grid_index  Index of pressure grid
 \param[in]     p_grid        New pressure grid
 \param[in]     interp_order  Interpolation order
 \param[in]     zeropadding   Allow zero padding
 \param[in]     verbosity     Verbosity levels
 */
void GriddedFieldPRegridHelper(Index& ing_min,
                               Index& ing_max,
                               ArrayOfLagrangeLogInterpolation& lag_p,
                               Matrix& itw,
                               GriddedField& gfraw_out,
                               const GriddedField& gfraw_in,
                               const Index p_grid_index,
                               ConstVectorView p_grid,
                               const Index& interp_order,
                               const Index& zeropadding,
                               const Verbosity& verbosity) {
  CREATE_OUT2;

  chk_griddedfield_gridname(gfraw_in, p_grid_index, "Pressure");

  out2 << "  Interpolation order: " << interp_order << "\n";

  const Vector& in_p_grid = gfraw_in.get_numeric_grid(p_grid_index);

  // Initialize output field. Set grids and copy grid names
  gfraw_out.set_grid(p_grid_index, Vector{p_grid});
  gfraw_out.set_grid_name(p_grid_index, gfraw_in.get_grid_name(p_grid_index));

  if (zeropadding) {
    if (in_p_grid[0] < p_grid[p_grid.nelem() - 1] ||
        in_p_grid[in_p_grid.nelem() - 1] > p_grid[0]) {
      ing_min = 0;
      ing_max = ing_min - 1;
    } else
      chk_interpolation_pgrids_loose_no_data_check(ing_min,
                                                   ing_max,
                                                   "Raw field to p_grid",
                                                   in_p_grid,
                                                   p_grid,
                                                   interp_order);
  } else {
    ing_min = 0;
    ing_max = p_grid.nelem() - 1;
    chk_interpolation_pgrids(
        "Raw field to p_grid", in_p_grid, p_grid, interp_order);
  }
  Index nelem_in_range = ing_max - ing_min + 1;

  // Calculate grid positions:
  if (nelem_in_range > 0) {
    lag_p = my_interp::lagrange_interpolation_list<LagrangeLogInterpolation>(p_grid[Range(ing_min, nelem_in_range)], in_p_grid, interp_order, 0.5);
    itw = interpweights(lag_p);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldPRegrid(  // WS Generic Output:
    GriddedField3& gfraw_out,
    // WS Input:
    const Vector& p_grid,
    // WS Generic Input:
    const GriddedField3& gfraw_in_orig,
    const Index& interp_order,
    const Index& zeropadding,
    const Verbosity& verbosity) {
  const GriddedField3* gfraw_in_pnt;
  GriddedField3 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField3& gfraw_in = *gfraw_in_pnt;

  const Index p_grid_index = 0;

  // Resize output GriddedField and copy all non-latitude/longitude grids
  gfraw_out.resize(
      p_grid.nelem(), gfraw_in.data.nrows(), gfraw_in.data.ncols());
  gfraw_out.set_grid(1, gfraw_in.get_numeric_grid(1));
  gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));
  gfraw_out.set_grid(2, gfraw_in.get_numeric_grid(2));
  gfraw_out.set_grid_name(2, gfraw_in.get_grid_name(2));
  
  ArrayOfLagrangeLogInterpolation lag_p;
  Matrix itw;  // nb. it is invalid to use this as it stands here...

  Index ing_min, ing_max;

  GriddedFieldPRegridHelper(ing_min,
                            ing_max,
                            lag_p,
                            itw,
                            gfraw_out,
                            gfraw_in,
                            p_grid_index,
                            p_grid,
                            interp_order,
                            zeropadding,
                            verbosity);

  // Interpolate:
  if (ing_max - ing_min < 0)
    gfraw_out.data = 0.;
  else if (ing_max - ing_min + 1 != p_grid.nelem()) {
    gfraw_out.data = 0.;
    for (Index i = 0; i < gfraw_in.data.nrows(); i++)
      for (Index j = 0; j < gfraw_in.data.ncols(); j++) {
        //                chk_interpolation_grids_loose_check_data(ing_min, ing_max,
        //                                                         "Raw field to p_grid",
        //                                                         gfraw_in.get_numeric_grid(p_grid_index),
        //                                                         p_grid, gfraw_in.data(joker, i, j));
        reinterp(gfraw_out.data(Range(ing_min, ing_max - ing_min + 1), i, j),
                 gfraw_in.data(joker, i, j),
                 itw,
                 lag_p);
      }
  } else
    for (Index i = 0; i < gfraw_in.data.nrows(); i++)
      for (Index j = 0; j < gfraw_in.data.ncols(); j++)
        reinterp(
          gfraw_out.data(joker, i, j), gfraw_in.data(joker, i, j), itw, lag_p);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldPRegrid(  // WS Generic Output:
    GriddedField4& gfraw_out,
    // WS Input:
    const Vector& p_grid,
    // WS Generic Input:
    const GriddedField4& gfraw_in_orig,
    const Index& interp_order,
    const Index& zeropadding,
    const Verbosity& verbosity) {
  const GriddedField4* gfraw_in_pnt;
  GriddedField4 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField4& gfraw_in = *gfraw_in_pnt;

  const Index p_grid_index = 1;

  // Resize output GriddedField and copy all non-latitude/longitude grids
  gfraw_out.resize(gfraw_in.data.nbooks(),
                   p_grid.nelem(),
                   gfraw_in.data.nrows(),
                   gfraw_in.data.ncols());
  gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
  gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));
  gfraw_out.set_grid(2, gfraw_in.get_numeric_grid(2));
  gfraw_out.set_grid_name(2, gfraw_in.get_grid_name(2));
  gfraw_out.set_grid(3, gfraw_in.get_numeric_grid(3));
  gfraw_out.set_grid_name(3, gfraw_in.get_grid_name(3));
  
  ArrayOfLagrangeLogInterpolation lag_p;
  Matrix itw;  // nb. it is invalid to use this as it stands here...

  Index ing_min, ing_max;

  GriddedFieldPRegridHelper(ing_min,
                            ing_max,
                            lag_p,
                            itw,
                            gfraw_out,
                            gfraw_in,
                            p_grid_index,
                            p_grid,
                            interp_order,
                            zeropadding,
                            verbosity);

  // Interpolate:
  if (ing_max - ing_min < 0)
    gfraw_out.data = 0.;
  else if (ing_max - ing_min + 1 != p_grid.nelem()) {
    gfraw_out.data = 0.;
    for (Index b = 0; b < gfraw_in.data.nbooks(); b++)
      for (Index i = 0; i < gfraw_in.data.nrows(); i++)
        for (Index j = 0; j < gfraw_in.data.ncols(); j++) {
          //                    chk_interpolation_grids_loose_check_data(ing_min, ing_max,
          //                                                             "Raw field to p_grid",
          //                                                             gfraw_in.get_numeric_grid(p_grid_index),
          //                                                             p_grid, gfraw_in.data(b, joker, i, j));
          reinterp(gfraw_out.data(b, Range(ing_min, ing_max - ing_min), i, j),
                   gfraw_in.data(b, joker, i, j),
                   itw,
                   lag_p);
        }
  } else
    for (Index b = 0; b < gfraw_in.data.nbooks(); b++)
      for (Index i = 0; i < gfraw_in.data.nrows(); i++)
        for (Index j = 0; j < gfraw_in.data.ncols(); j++)
          reinterp(gfraw_out.data(b, joker, i, j),
                   gfraw_in.data(b, joker, i, j),
                   itw,
                   lag_p);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldPRegrid(  // WS Generic Output:
    ArrayOfGriddedField3& agfraw_out,
    // WS Input:
    const Vector& p_grid,
    // WS Generic Input:
    const ArrayOfGriddedField3& agfraw_in,
    const Index& interp_order,
    const Index& zeropadding,
    const Verbosity& verbosity) {
  agfraw_out.resize(agfraw_in.nelem());

  for (Index i = 0; i < agfraw_in.nelem(); i++) {
    GriddedFieldPRegrid(agfraw_out[i],
                        p_grid,
                        agfraw_in[i],
                        interp_order,
                        zeropadding,
                        verbosity);
  }
}

//! Calculate grid positions and interpolations weights for GriddedFieldLatLonRegrid
/*
 This helper function is used by all GriddedFieldLatLonRegrid WSMs to do the common
 calculation of the grid positions and interpolation weights for latitudes and longitudes.

 \param[out]    gp_lat          Latitude grid positions
 \param[out]    gp_lon          Longitude grid positions
 \param[out]    itw             Interpolation weights
 \param[in,out] gfraw_out       Output GriddedField
 \param[in]     gfraw_in        Input GriddedField
 \param[in]     lat_grid_index  Index of latitude grid
 \param[in]     lon_grid_index  Index of longitude grid
 \param[in]     lat_true        New latitude grid
 \param[in]     lon_true        New longitude grid
 \param[in]     interp_order    Interpolation order
 \param[in]     verbosity       Verbosity levels
 */
void GriddedFieldLatLonRegridHelper(ArrayOfLagrangeInterpolation& lag_lat,
                                    ArrayOfLagrangeCyclic0to360Interpolation& lag_lon,
                                    Tensor4& itw,
                                    GriddedField& gfraw_out,
                                    const GriddedField& gfraw_in,
                                    const Index lat_grid_index,
                                    const Index lon_grid_index,
                                    ConstVectorView lat_true,
                                    ConstVectorView lon_true,
                                    const Index& interp_order,
                                    const Verbosity& verbosity) {
  CREATE_OUT2;

  ARTS_USER_ERROR_IF (!lat_true.nelem(),
                      "The new latitude grid is not allowed to be empty.");
  ARTS_USER_ERROR_IF (!lon_true.nelem(),
                      "The new longitude grid is not allowed to be empty.");

  chk_griddedfield_gridname(gfraw_in, lat_grid_index, "Latitude");
  chk_griddedfield_gridname(gfraw_in, lon_grid_index, "Longitude");
  ARTS_USER_ERROR_IF (gfraw_in.get_grid_size(lat_grid_index) == 1 ||
      gfraw_in.get_grid_size(lon_grid_index) == 1,
                      "Raw data has to be true 3D data (nlat>1 and nlon>1).");

  out2 << "  Interpolation order: " << interp_order << "\n";

  const Vector& in_lat_grid = gfraw_in.get_numeric_grid(lat_grid_index);
  const Vector& in_lon_grid = gfraw_in.get_numeric_grid(lon_grid_index);

  // Initialize output field. Set grids and copy grid names
  gfraw_out.set_grid(lat_grid_index, Vector{lat_true});
  gfraw_out.set_grid_name(lat_grid_index,
                          gfraw_in.get_grid_name(lat_grid_index));
  gfraw_out.set_grid(lon_grid_index, Vector{lon_true});
  gfraw_out.set_grid_name(lon_grid_index,
                          gfraw_in.get_grid_name(lon_grid_index));

  chk_interpolation_grids(
      "Raw field to lat_grid, 3D case", in_lat_grid, lat_true, interp_order);

  // Calculate grid positions:
  lag_lat = my_interp::lagrange_interpolation_list<LagrangeInterpolation>(lat_true, in_lat_grid, interp_order);
  lag_lon = my_interp::lagrange_interpolation_list<LagrangeCyclic0to360Interpolation>(lon_true, in_lon_grid, interp_order, 0.5);
  itw = interpweights(lag_lat, lag_lon);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(  // WS Generic Output:
    GriddedField2& gfraw_out,
    // WS Input:
    const Vector& lat_true,
    const Vector& lon_true,
    // WS Generic Input:
    const GriddedField2& gfraw_in_orig,
    const Index& interp_order,
    const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF (!lat_true.nelem(),
                      "The new latitude grid is not allowed to be empty.");
  ARTS_USER_ERROR_IF (!lon_true.nelem(),
                      "The new longitude grid is not allowed to be empty.");

  const GriddedField2* gfraw_in_pnt;
  GriddedField2 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField2& gfraw_in = *gfraw_in_pnt;

  const Index lat_grid_index = 0;
  const Index lon_grid_index = 1;

  ARTS_USER_ERROR_IF (gfraw_in.get_grid_size(lat_grid_index) < 2 ||
      gfraw_in.get_grid_size(lon_grid_index) < 2,
      "Raw data has to be true 3D data (nlat>1 and nlon>1).\n"
      "Use GriddedFieldLatLonExpand to convert 1D or 2D data to 3D!\n")

  // Resize output GriddedField
  gfraw_out.resize(lat_true.nelem(), lon_true.nelem());

  ArrayOfLagrangeInterpolation lag_lat;
  ArrayOfLagrangeCyclic0to360Interpolation lag_lon;
  Tensor4 itw;

  // If lon grid is cyclic, the data values at 0 and 360 must match
  const Vector& in_lat_grid =
      gfraw_in.get_numeric_grid(lat_grid_index);
  const Vector& in_lon_grid =
      gfraw_in.get_numeric_grid(lon_grid_index);

  if (is_lon_cyclic(in_lon_grid)) {
    for (Index lat = 0; lat < in_lat_grid.nelem(); lat++) {
      ARTS_USER_ERROR_IF (!is_same_within_epsilon(gfraw_in.data(lat, 0),
                                  gfraw_in.data(lat, in_lon_grid.nelem() - 1),
                                  EPSILON_LON_CYCLIC),
           "Data values at 0 and 360 degrees for a cyclic longitude grid must match: \n"
          , "Mismatch at latitude index    : ", lat, " ("
          , in_lat_grid[lat], " degrees)\n"
          , "Value at 0 degrees longitude  : ", gfraw_in.data(lat, 0)
          , "\n"
          , "Value at 360 degrees longitude: "
          , gfraw_in.data(lat, in_lon_grid.nelem() - 1), "\n"
          , "Difference                    : "
          , gfraw_in.data(lat, in_lon_grid.nelem() - 1) -
                  gfraw_in.data(lat, 0)
          , "\n"
          , "Allowed difference            : ", EPSILON_LON_CYCLIC)
    }
  }

  GriddedFieldLatLonRegridHelper(lag_lat,
                                 lag_lon,
                                 itw,
                                 gfraw_out,
                                 gfraw_in,
                                 lat_grid_index,
                                 lon_grid_index,
                                 lat_true,
                                 lon_true,
                                 interp_order,
                                 verbosity);

  // Interpolate:
  reinterp(gfraw_out.data, gfraw_in.data, itw, lag_lat, lag_lon);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(  // WS Generic Output:
    GriddedField3& gfraw_out,
    // WS Input:
    const Vector& lat_true,
    const Vector& lon_true,
    // WS Generic Input:
    const GriddedField3& gfraw_in_orig,
    const Index& interp_order,
    const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF (!lat_true.nelem(),
    "The new latitude grid is not allowed to be empty.");
  ARTS_USER_ERROR_IF (!lon_true.nelem(),
    "The new longitude grid is not allowed to be empty.");

  const GriddedField3* gfraw_in_pnt;
  GriddedField3 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField3& gfraw_in = *gfraw_in_pnt;

  const Index lat_grid_index = 1;
  const Index lon_grid_index = 2;

  ARTS_USER_ERROR_IF (gfraw_in.get_grid_size(lat_grid_index) < 2 ||
      gfraw_in.get_grid_size(lon_grid_index) < 2,
      "Raw data has to be true 3D data (nlat>1 and nlon>1).\n"
      "Use GriddedFieldLatLonExpand to convert 1D or 2D data to 3D!\n")

  // Resize output GriddedField and copy all non-latitude/longitude grids
  gfraw_out.resize(gfraw_in.data.npages(), lat_true.nelem(), lon_true.nelem());
  gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
  gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));
  
  ArrayOfLagrangeInterpolation lag_lat;
  ArrayOfLagrangeCyclic0to360Interpolation lag_lon;
  Tensor4 itw;

  // If lon grid is cyclic, the data values at 0 and 360 must match
  const Vector& in_grid0 = gfraw_in.get_numeric_grid(0);
  const Vector& in_lat_grid =
      gfraw_in.get_numeric_grid(lat_grid_index);
  const Vector& in_lon_grid =
      gfraw_in.get_numeric_grid(lon_grid_index);

  if (is_lon_cyclic(in_lon_grid)) {
    for (Index g0 = 0; g0 < in_grid0.nelem(); g0++)
      for (Index lat = 0; lat < in_lat_grid.nelem(); lat++) {
        ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                gfraw_in.data(g0, lat, 0),
                gfraw_in.data(g0, lat, in_lon_grid.nelem() - 1),
                EPSILON_LON_CYCLIC),
             "Data values at 0 and 360 degrees for a cyclic longitude grid must match: \n"
             , "Mismatch at 1st grid index    : " , g0 , " (" , in_grid0[g0]
             , ")\n"
             , "         at latitude index    : " , lat , " ("
             , in_lat_grid[lat] , " degrees)\n"
             , "Value at 0 degrees longitude  : " , gfraw_in.data(g0, lat, 0)
             , "\n"
             , "Value at 360 degrees longitude: "
             , gfraw_in.data(g0, lat, in_lon_grid.nelem() - 1) , "\n"
             , "Difference                    : "
             , gfraw_in.data(g0, lat, in_lon_grid.nelem() - 1) -
                    gfraw_in.data(g0, lat, 0)
             , "\n"
             , "Allowed difference            : " , EPSILON_LON_CYCLIC)
      }
  }

  GriddedFieldLatLonRegridHelper(lag_lat,
                                 lag_lon,
                                 itw,
                                 gfraw_out,
                                 gfraw_in,
                                 lat_grid_index,
                                 lon_grid_index,
                                 lat_true,
                                 lon_true,
                                 interp_order,
                                 verbosity);

  // Interpolate:
  for (Index i = 0; i < gfraw_in.data.npages(); i++)
    reinterp(gfraw_out.data(i, joker, joker),
             gfraw_in.data(i, joker, joker),
             itw,
             lag_lat,
             lag_lon);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(  // WS Generic Output:
    GriddedField4& gfraw_out,
    // WS Input:
    const Vector& lat_true,
    const Vector& lon_true,
    // WS Generic Input:
    const GriddedField4& gfraw_in_orig,
    const Index& interp_order,
    const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF (!lat_true.nelem(),
                      "The new latitude grid is not allowed to be empty.");
  ARTS_USER_ERROR_IF (!lon_true.nelem(),
                      "The new longitude grid is not allowed to be empty.");

  const GriddedField4* gfraw_in_pnt;
  GriddedField4 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  const GriddedField4& gfraw_in = *gfraw_in_pnt;

  const Index lat_grid_index = 2;
  const Index lon_grid_index = 3;

  ARTS_USER_ERROR_IF (gfraw_in.get_grid_size(lat_grid_index) < 2 ||
                      gfraw_in.get_grid_size(lon_grid_index) < 2,
      "Raw data has to be true 3D data (nlat>1 and nlon>1).\n"
      "Use GriddedFieldLatLonExpand to convert 1D or 2D data to 3D!\n")

  // Resize output GriddedField and copy all non-latitude/longitude grids
  gfraw_out.resize(gfraw_in.data.nbooks(),
                   gfraw_in.data.npages(),
                   lat_true.nelem(),
                   lon_true.nelem());
  gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
  gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));
  gfraw_out.set_grid(1, gfraw_in.get_numeric_grid(1));
  gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));
  
  ArrayOfLagrangeInterpolation lag_lat;
  ArrayOfLagrangeCyclic0to360Interpolation lag_lon;
  Tensor4 itw;

  GriddedFieldLatLonRegridHelper(lag_lat,
                                 lag_lon,
                                 itw,
                                 gfraw_out,
                                 gfraw_in,
                                 lat_grid_index,
                                 lon_grid_index,
                                 lat_true,
                                 lon_true,
                                 interp_order,
                                 verbosity);

  // If lon grid is cyclic, the data values at 0 and 360 must match
  const Vector& in_grid0 = gfraw_in.get_numeric_grid(0);
  const Vector& in_grid1 = gfraw_in.get_numeric_grid(1);
  const Vector& in_lat_grid =
      gfraw_in.get_numeric_grid(lat_grid_index);
  const Vector& in_lon_grid =
      gfraw_in.get_numeric_grid(lon_grid_index);

  if (is_lon_cyclic(in_lon_grid)) {
    for (Index g0 = 0; g0 < in_grid0.nelem(); g0++)
      for (Index g1 = 0; g1 < in_grid1.nelem(); g1++)
        for (Index lat = 0; lat < in_lat_grid.nelem(); lat++) {
          ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                  gfraw_in.data(g0, g1, lat, 0),
                  gfraw_in.data(g0, g1, lat, in_lon_grid.nelem() - 1),
                  EPSILON_LON_CYCLIC),
              "Data values at 0 and 360 degrees for a cyclic longitude grid must match: \n"
              "Mismatch at 1st grid index    : ", g0, " ("
              , in_grid0[g0], ")\n"
              , "         at 2nd grid index    : ", g1, " ("
              , in_grid1[g1], ")\n"
              , "         at latitude index    : ", lat, " ("
              , in_lat_grid[lat], " degrees)\n"
              , "Value at 0 degrees longitude  : "
              , gfraw_in.data(g0, g1, lat, 0), "\n"
              , "Value at 360 degrees longitude: "
              , gfraw_in.data(g0, g1, lat, in_lon_grid.nelem() - 1), "\n"
              , "Difference                    : "
              , gfraw_in.data(g0, g1, lat, in_lon_grid.nelem() - 1) -
                      gfraw_in.data(g0, g1, lat, 0)
              , "\n"
              , "Allowed difference            : ", EPSILON_LON_CYCLIC)
        }
  }

  // Interpolate:
  for (Index i = 0; i < gfraw_in.data.nbooks(); i++)
    for (Index j = 0; j < gfraw_in.data.npages(); j++)
      reinterp(gfraw_out.data(i, j, joker, joker),
               gfraw_in.data(i, j, joker, joker),
               itw,
               lag_lat,
               lag_lon);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(  // WS Generic Output:
    ArrayOfGriddedField3& agfraw_out,
    // WS Input:
    const Vector& lat_true,
    const Vector& lon_true,
    // WS Generic Input:
    const ArrayOfGriddedField3& agfraw_in,
    const Index& interp_order,
    const Verbosity& verbosity) {
  agfraw_out.resize(agfraw_in.nelem());

  for (Index i = 0; i < agfraw_in.nelem(); i++) {
    GriddedFieldLatLonRegrid(agfraw_out[i],
                             lat_true,
                             lon_true,
                             agfraw_in[i],
                             interp_order,
                             verbosity);
  }
}

//! Calculate grid positions and interpolations weights for GriddedFieldZToPRegrid
/*
 This helper function is used by GriddedFieldZToPRegrid WSM to do the common
 calculation of the grid positions and interpolation weights for latitudes and longitudes.

 \param[out]    ing_min       Index in the new grid with first value covered
 by the old grid.
 \param[out]    ing_max       Index in the new grid with last value covered
 by the old grid.
 \param[out]    gp_p          Altitude grid positions
 \param[out]    itw           Interpolation weights
 \param[in]     gfraw_in      Input GriddedField
 \param[in]     z_grid_index  Index of altitude grid
 \param[in]     z_grid        New z_grid grid
 \param[in]     interp_order  Interpolation order
 \param[in]     zeropadding   Allow zero padding
 \param[in]     verbosity     Verbosity levels
 */
void GriddedFieldZToPRegridHelper(Index& ing_min,
                                  Index& ing_max,
                                  ArrayOfLagrangeInterpolation& lag_p,
                                  Matrix& itw,
                                  const GriddedField& gfraw_in,
                                  const Index z_grid_index,
                                  ConstVectorView z_grid,
                                  const Index& interp_order,
                                  const Index& zeropadding,
                                  const Verbosity& verbosity) {
  CREATE_OUT2;

  chk_griddedfield_gridname(gfraw_in, z_grid_index, "Altitude");

  out2 << "  Interpolation order: " << interp_order << "\n";

  const Vector& in_z_grid = gfraw_in.get_numeric_grid(z_grid_index);

  if (zeropadding) {
    if (in_z_grid[0] > z_grid[z_grid.nelem() - 1] ||
        in_z_grid[in_z_grid.nelem() - 1] < z_grid[0]) {
      ing_min = 0;
      ing_max = ing_min - 1;
    } else
      chk_interpolation_grids_loose_no_data_check(ing_min,
                                                  ing_max,
                                                  "Raw field to z_grid",
                                                  in_z_grid,
                                                  z_grid,
                                                  interp_order);
  } else {
    ing_min = 0;
    ing_max = z_grid.nelem() - 1;
    chk_interpolation_pgrids(
        "Raw field to p_grid", in_z_grid, z_grid, interp_order);
  }

  Index nelem_in_range = ing_max - ing_min + 1;

  // Calculate grid positions:
  if (nelem_in_range > 0) {
    lag_p = my_interp::lagrange_interpolation_list<LagrangeInterpolation>(z_grid[Range(ing_min, nelem_in_range)], in_z_grid, interp_order);
    itw = interpweights(lag_p);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldZToPRegrid(   // WS Generic Output:
    GriddedField3& gfraw_out,  //grid in P
    // WS Input:
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& z_field,
    // WS Generic Input:
    const GriddedField3& gfraw_in_orig,  //grid in Z
    const Index& interp_order,           // Only linear interpolation allowed
    const Index& zeropadding,
    const Verbosity& verbosity) {
  // z_field must be of the same size as its grids
  ARTS_USER_ERROR_IF (!((z_field.npages() == p_grid.nelem() &&
         z_field.nrows() == lat_grid.nelem()) &&
        z_field.ncols() == lon_grid.nelem()),
        "*z_field* must be of the same size as *p_grid*, *lat_grid*, and *lon_grid* in *GriddedFieldZToPRegrid*.");

  // Must name the dimension "Altitude" to ensure user is aware of what they are doing.
  chk_griddedfield_gridname(gfraw_in_orig, 0, "Altitude");

  // Each lat and lon grid must be identical between z_field and in field
  const Vector& lat_in = gfraw_in_orig.get_numeric_grid(1);
  const Vector& lon_in = gfraw_in_orig.get_numeric_grid(2);

  ARTS_USER_ERROR_IF (lat_grid.nelem() != lat_in.nelem() || lon_grid.nelem() != lon_in.nelem(),
        "Gridding of field to regrid is bad.\n*GriddedFieldZToPRegrid* requires latitude and longitude to be on the same grid as *z_field*.");
  for (Index ii = 0; ii < lat_grid.nelem(); ii++)
    ARTS_USER_ERROR_IF (lat_grid[ii] != lat_in[ii],
          "Gridding of field to regrid is bad.\n*GriddedFieldZToPRegrid* requires latitude and longitude of the gridded field to be the same as for *z_field*.");
  for (Index ii = 0; ii < lon_grid.nelem(); ii++)
    ARTS_USER_ERROR_IF (lon_grid[ii] != lon_in[ii],
          "Gridding of field to regrid is bad.\n*GriddedFieldZToPRegrid* requires latitude and longitude of the gridded field to be the same as for *z_field*.");

  // Pointer in case output is input variable (same memory allocated)
  const GriddedField3* gfraw_in_pnt;
  GriddedField3 gfraw_in_copy;

  if (&gfraw_in_orig == &gfraw_out) {
    gfraw_in_copy = gfraw_in_orig;
    gfraw_in_pnt = &gfraw_in_copy;
  } else
    gfraw_in_pnt = &gfraw_in_orig;

  // Now output and input are separate variables (not allocating the same memory)
  const GriddedField3& gfraw_in = *gfraw_in_pnt;

  // Right size and order
  gfraw_out.resize(p_grid.nelem(), lat_grid.nelem(), lon_grid.nelem());
  gfraw_out.set_grid(0, p_grid);
  gfraw_out.set_grid_name(0, "Pressure");
  gfraw_out.set_grid(1, lat_grid);
  gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));
  gfraw_out.set_grid(2, lon_grid);
  gfraw_out.set_grid_name(2, gfraw_in.get_grid_name(2));
  gfraw_out.data = 0.;

  ArrayOfLagrangeInterpolation lag_p;
  Matrix itw;  // nb. it is invalid to use this as it stands here...

  Index ing_min, ing_max;

  for (Index lat_index = 0; lat_index < lat_grid.nelem(); lat_index++) {
    for (Index lon_index = 0; lon_index < lon_grid.nelem(); lon_index++) {
      const Vector z_out{z_field(joker, lat_index, lon_index)};

      GriddedFieldZToPRegridHelper(ing_min,
                                   ing_max,
                                   lag_p,
                                   itw,
                                   gfraw_in,
                                   0,
                                   z_out,
                                   interp_order,
                                   zeropadding,
                                   verbosity);

      if (ing_max - ing_min >= 0) {
        Range r = joker;

        if (ing_max - ing_min + 1 != z_out.nelem()) {
          r = Range(ing_min, ing_max - ing_min + 1);
        }

        reinterp(gfraw_out.data(r, lat_index, lon_index),
                 gfraw_in.data(joker, lat_index, lon_index),
                 itw,
                 lag_p);
      }
    }
  }
}

// Workspace method, doxygen header will be auto-generated.
// 2007-07-25 Stefan Buehler
void atm_fields_compactFromMatrix(  // WS Output:
    GriddedField4& af,              // atm_fields_compact
    // WS Input:
    // WS Generic Input:
    const Matrix& im,
    // Control Parameters:
    const ArrayOfString& field_names,
    const Verbosity&) {
  ARTS_USER_ERROR_IF (1 != 3,
    "Atmospheric dimension must be 1.")

  const Index np = im.nrows();      // Number of pressure levels.
  const Index nf = im.ncols() - 1;  // Total number of fields.
  ArrayOfIndex f_1;                 // indices of non-ignored fields
      // All fields called "ignore" will be ignored.
  String fn_upper;  // Temporary variable to hold upper case field_names.

  ARTS_USER_ERROR_IF (field_names.nelem() != nf,
    "Cannot extract fields from Matrix.\n"
    "*field_names* must have one element less than there are\n"
    "matrix columns.")

  // Remove additional fields from the field_names. All fields that
  // are flagged by 'ignore' in the field names, small or large letters,
  // are removed.
  for (Index f = 0; f < field_names.nelem(); f++) {
    fn_upper = field_names[f];
    fn_upper.toupper();
    //cout << "fieldname[" << f << "]: " << fn_upper;
    if (fn_upper != "IGNORE") {
      f_1.push_back(f);
      //cout << " put as element " << f_1.size()-1 << " into selection\n";
    }
  }

  // Copy required field_names to a new variable called field_names_1
  Index nf_1 = f_1.size();  // Number of required fields.
  ArrayOfString field_names_1(nf_1);
  for (Index f = 0; f < nf_1; f++) {
    field_names_1[f] = field_names[f_1[f]];
    //cout << "new fieldname[" << f << "] (formerly [" << f_1[f] << "]): "
    //     << field_names_1[f] << "\n";
  }

  //  out3 << "Copying *" << im_name << "* to *atm_fields_compact*.\n";

  af.set_grid(GFIELD4_FIELD_NAMES, field_names_1);

  af.set_grid(GFIELD4_P_GRID, Vector{im(Range(joker), 0)});

  af.set_grid(GFIELD4_LAT_GRID, Vector());
  af.set_grid(GFIELD4_LON_GRID, Vector());

  af.resize(nf_1, np, 1, 1);  // Resize it according to the required fields
  for (Index f = 0; f < nf_1; f++)
    af.data(f, Range(joker), 0, 0) = im(Range(joker), f_1[f] + 1);
}

// Workspace method, doxygen header is auto-generated.
// 2007-07-31 Stefan Buehler
// 2011-05-04 Adapted by Gerrit Holl
void atm_fields_compactAddConstant(  // WS Output:
    GriddedField4& af,
    // Control Parameters:
    const String& name,
    const Numeric& value,
    const Index& prepend,
    const ArrayOfString& condensibles,
    const Verbosity& verbosity) {
  Index nf;  // Will hold new size

  // Add book
  atm_fields_compactExpand(af, nf, name, prepend, verbosity);

  if (condensibles.nelem()) {
    const Tensor4& vmrs = af.data;
    const ArrayOfString& species = af.get_string_grid(GFIELD4_FIELD_NAMES);
    Tensor3 condensible_sum(vmrs.npages(), vmrs.nrows(), vmrs.ncols(), 1.);
    for (Index c = 0; c < condensibles.nelem(); c++) {
      bool species_found = false;
      for (Index i = 0; !species_found && i < species.nelem(); i++) {
        if (species[i] == condensibles[c]) {
          condensible_sum -= vmrs(i, joker, joker, joker);
          species_found = true;
        }
      }
      ARTS_USER_ERROR_IF (!species_found,
        "Condensible species \"", condensibles[c], "\" not found "
        "in input data.")
    }
    condensible_sum *= value;
    // Add the constant value:
    if (prepend)
      af.data(0, joker, joker, joker) = condensible_sum;
    else
      af.data(nf - 1, joker, joker, joker) = condensible_sum;
  } else {
    // Add the constant value:
    if (prepend)
      af.data(0, joker, joker, joker) = value;
    else
      af.data(nf - 1, joker, joker, joker) = value;
  }
}

// Workspace method, doxygen header is auto-generated
// 2011-05-02 Gerrit Holl
void atm_fields_compactAddSpecies(  // WS Output:
    GriddedField4& atm_fields_compact,
    // WS Generic Input:
    const String& name,
    const GriddedField3& species,
    const Index& prepend,
    const Verbosity& verbosity) {
  ARTS_ASSERT(atm_fields_compact.checksize());
  ARTS_ASSERT(species.checksize());

  ConstVectorView af_p_grid =
      atm_fields_compact.get_numeric_grid(GFIELD4_P_GRID);
  ConstVectorView af_lat_grid =
      atm_fields_compact.get_numeric_grid(GFIELD4_LAT_GRID);
  ConstVectorView af_lon_grid =
      atm_fields_compact.get_numeric_grid(GFIELD4_LON_GRID);
  ConstVectorView sp_p_grid = species.get_numeric_grid(GFIELD3_P_GRID);
  ConstVectorView sp_lat_grid = species.get_numeric_grid(GFIELD3_LAT_GRID);
  ConstVectorView sp_lon_grid = species.get_numeric_grid(GFIELD3_LON_GRID);

  Index new_n_fields;  // To be set in next line
  atm_fields_compactExpand(
      atm_fields_compact, new_n_fields, name, prepend, verbosity);

  const Index insert_pos = (prepend) ? 0 : new_n_fields - 1;

  // Interpolate species to atm_fields_compact

  // Common for all dim
  chk_interpolation_grids(
      "species p_grid to atm_fields_compact p_grid", sp_p_grid, af_p_grid);
  ArrayOfGridPos p_gridpos(af_p_grid.nelem());
  // gridpos(p_gridpos, sp_p_grid, af_p_grid);
  p2gridpos(p_gridpos, sp_p_grid, af_p_grid);

  if (sp_lat_grid.nelem() > 1) {
    // Common for all dim>=2
    chk_interpolation_grids("species lat_grid to atm_fields_compact lat_grid",
                            sp_lat_grid,
                            af_lat_grid);
    ArrayOfGridPos lat_gridpos(af_lat_grid.nelem());
    gridpos(lat_gridpos, sp_lat_grid, af_lat_grid);

    if (sp_lon_grid.nelem() > 1) {  // 3D-case
      chk_interpolation_grids("species lon_grid to atm_fields_compact lon_grid",
                              sp_lon_grid,
                              af_lon_grid);
      ArrayOfGridPos lon_gridpos(af_lon_grid.nelem());
      gridpos(lon_gridpos, sp_lon_grid, af_lon_grid);

      Tensor4 itw(
          p_gridpos.nelem(), lat_gridpos.nelem(), lon_gridpos.nelem(), 8);
      interpweights(itw, p_gridpos, lat_gridpos, lon_gridpos);

      Tensor3 newfield(
          af_p_grid.nelem(), af_lat_grid.nelem(), af_lon_grid.nelem());
      interp(newfield, itw, species.data, p_gridpos, lat_gridpos, lon_gridpos);

      atm_fields_compact.data(insert_pos, joker, joker, joker) = newfield;
    } else {  // 2D-case

      Tensor3 itw(p_gridpos.nelem(), lat_gridpos.nelem(), 4);
      interpweights(itw, p_gridpos, lat_gridpos);

      Matrix newfield(af_p_grid.nelem(), af_lat_grid.nelem());
      interp(
          newfield, itw, species.data(joker, joker, 0), p_gridpos, lat_gridpos);

      atm_fields_compact.data(insert_pos, joker, joker, 0) = newfield;
    }
  } else {  // 1D-case
    Matrix itw(p_gridpos.nelem(), 2);
    interpweights(itw, p_gridpos);

    Vector newfield(af_p_grid.nelem());
    interp(newfield, itw, species.data(joker, 0, 0), p_gridpos);

    atm_fields_compact.data(insert_pos, joker, 0, 0) = newfield;
  }
}

// Workspace method, doxygen header is auto-generated
// 2015-06-30 Jana Mendrok
void atm_fields_compactCleanup(  //WS Output:
    GriddedField4& atm_fields_compact,
    //WS Input:
    const Numeric& threshold,
    const Verbosity&) {
  ARTS_ASSERT(atm_fields_compact.checksize());
  Tensor4View afd = atm_fields_compact.data;

  // Check that the data tensor does not contain realistically low (e.g.
  // negative) values. Values smaller than threshold will be set to 0).
  // Ignore T and z, though.
  for (Index i = 0; i < afd.nbooks(); i++)
    if (atm_fields_compact.get_string_grid(GFIELD4_FIELD_NAMES)[i] != "T" &&
        atm_fields_compact.get_string_grid(GFIELD4_FIELD_NAMES)[i] != "z")
      for (Index j = 0; j < afd.npages(); j++)
        for (Index k = 0; k < afd.nrows(); k++)
          for (Index l = 0; l < afd.ncols(); l++)
            if (afd(i, j, k, l) < threshold) afd(i, j, k, l) = 0.0;
}

// Workspace method, doxygen header is auto-generated
void atm_fields_compactCreateFromField(  // WS Output:
    GriddedField4& atm_fields_compact,
    // WS Generic Input:
    const String& name,
    const GriddedField3& field,
    const Verbosity&) {
  ARTS_ASSERT(field.checksize());

  ConstVectorView sp_p_grid = field.get_numeric_grid(GFIELD3_P_GRID);
  ConstVectorView sp_lat_grid = field.get_numeric_grid(GFIELD3_LAT_GRID);
  ConstVectorView sp_lon_grid = field.get_numeric_grid(GFIELD3_LON_GRID);
  ArrayOfString sp_name_grid(1);
  sp_name_grid[0] = name;

  atm_fields_compact.set_grid(0, sp_name_grid);
  atm_fields_compact.set_grid(1, Vector{sp_p_grid});
  atm_fields_compact.set_grid(2, Vector{sp_lat_grid});
  atm_fields_compact.set_grid(3, Vector{sp_lon_grid});

  atm_fields_compact.data.resize(
      1, sp_p_grid.nelem(), sp_lat_grid.nelem(), sp_lon_grid.nelem());

  atm_fields_compact.data(0, joker, joker, joker) = field.data;
}

// Workspace method, doxygen header is auto-generated
// 2011-05-11 Gerrit Holl
void batch_atm_fields_compactAddConstant(  // WS Output:
    ArrayOfGriddedField4& batch_atm_fields_compact,
    // WS Generic Input:
    const String& name,
    const Numeric& value,
    const Index& prepend,
    const ArrayOfString& condensibles,
    const Verbosity& verbosity) {
  for (Index i = 0; i < batch_atm_fields_compact.nelem(); i++) {
    atm_fields_compactAddConstant(batch_atm_fields_compact[i],
                                  name,
                                  value,
                                  prepend,
                                  condensibles,
                                  verbosity);
  }
}

// Workspace method, doxygen header is auto-generated
// 2011-05-09 Gerrit Holl
void batch_atm_fields_compactAddSpecies(  // WS Output:
    ArrayOfGriddedField4& batch_atm_fields_compact,
    // WS Generic Input:
    const String& name,
    const GriddedField3& species,
    const Index& prepend,
    const Verbosity& verbosity) {
  const Index nelem = batch_atm_fields_compact.nelem();

  String fail_msg;
  bool failed = false;

  // Parallelise this for-loop (some interpolation is being done, so it may
  // be beneficial)
#pragma omp parallel for if (!arts_omp_in_parallel() && \
                             nelem >= arts_omp_get_max_threads())
  for (Index i = 0; i < nelem; i++) {
    try {
      atm_fields_compactAddSpecies(
          batch_atm_fields_compact[i], name, species, prepend, verbosity);
    } catch (const std::exception& e) {
#pragma omp critical(batch_atm_fields_compactAddSpecies_fail)
      {
        fail_msg = e.what();
        failed = true;
      }
    }
  }

  ARTS_USER_ERROR_IF (failed, fail_msg);
}

// Workspace method, doxygen header is auto-generated
// 2015-06-30 Jana Mendrok
void batch_atm_fields_compactCleanup(  //WS Output:
    ArrayOfGriddedField4& batch_atm_fields_compact,
    //WS Input:
    const Numeric& threshold,
    const Verbosity& verbosity) {
  for (Index i = 0; i < batch_atm_fields_compact.nelem(); i++) {
    atm_fields_compactCleanup(
        batch_atm_fields_compact[i], threshold, verbosity);
  }
}

// Workspace method, doxygen header is auto-generated.
void batch_atm_fields_compactFromArrayOfMatrix(  // WS Output:
    ArrayOfGriddedField4& batch_atm_fields_compact,
    // WS Input:
    // WS Generic Input:
    const ArrayOfMatrix& am,
    // Control Parameters:
    const ArrayOfString& field_names,
    const Verbosity& verbosity) {
  const Index amnelem = am.nelem();

  ARTS_USER_ERROR_IF (amnelem == 0,
      "No elements in atmospheric scenario batch.\n"
      "Check, whether any batch atmosphere file has been read!")

  // We use the existing WSMs atm_fields_compactFromMatrix and
  // atm_fields_compactAddConstant to do most of the work.

  // Make output variable the proper size:
  batch_atm_fields_compact.resize(amnelem);

  String fail_msg;
  bool failed = false;

  // Loop the batch cases:
#pragma omp parallel for if (!arts_omp_in_parallel() && \
                             amnelem >= arts_omp_get_max_threads())
  for (Index i = 0; i < amnelem; ++i) {
    // Skip remaining iterations if an error occurred
    if (failed) continue;

    // All the input variables are visible here, despite the
    // "default(none)". The reason is that they are return by
    // reference arguments of this function, which are shared
    // automatically.

    // The try block here is necessary to correctly handle
    // exceptions inside the parallel region.
    try {
      atm_fields_compactFromMatrix(batch_atm_fields_compact[i],
                                   am[i],
                                   field_names,
                                   verbosity);
    } catch (const std::exception& e) {
#pragma omp critical(batch_atm_fields_compactFromArrayOfMatrix_fail)
      {
        fail_msg = e.what();
        failed = true;
      }
    }
  }

  ARTS_USER_ERROR_IF (failed, fail_msg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsAndParticleBulkPropFieldFromCompact(  // WS Output:
    Vector& p_grid,
    Vector& lat_grid,
    Vector& lon_grid,
    Tensor3& t_field,
    Tensor3& z_field,
    Tensor4& vmr_field,
    Tensor4& particle_bulkprop_field,
    ArrayOfString& particle_bulkprop_names,

    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const GriddedField4& atm_fields_compact,
    const String& delim,
    const Numeric& p_min,
    // Control parameters:
    const Index& check_gridnames,
    const Verbosity&) {
  // Make a handle on atm_fields_compact to save typing:
  const GriddedField4& c = atm_fields_compact;

  // Check if the grids in our data match 3
  // (throws an error if the dimensionality is not correct):
  chk_atm_grids(3,
                c.get_numeric_grid(GFIELD4_P_GRID),
                c.get_numeric_grid(GFIELD4_LAT_GRID),
                c.get_numeric_grid(GFIELD4_LON_GRID));

  // Optional check for gridnames.
  if (check_gridnames == 1) {
    chk_griddedfield_gridname(c, 1, "Pressure");
    chk_griddedfield_gridname(c, 2, "Latitude");
    chk_griddedfield_gridname(c, 3, "Longitude");
  }

  const Index nf = c.get_grid_size(GFIELD4_FIELD_NAMES);
  const Index np = c.get_grid_size(GFIELD4_P_GRID);
  //  const Index nlat = c.get_grid_size(GFIELD4_LAT_GRID);
  //  const Index nlon = c.get_grid_size(GFIELD4_LON_GRID);
  Index nlat = c.get_grid_size(GFIELD4_LAT_GRID);
  Index nlon = c.get_grid_size(GFIELD4_LON_GRID);
  if (nlat == 0) nlat++;
  if (nlon == 0) nlon++;

  // Grids:
  p_grid = c.get_numeric_grid(GFIELD4_P_GRID);
  lat_grid = c.get_numeric_grid(GFIELD4_LAT_GRID);
  lon_grid = c.get_numeric_grid(GFIELD4_LON_GRID);

  // Reduce p_grid to region below p_min
  Index l = np - 1;
  bool search_toa = true;
  while (search_toa && l > 0) {
    if (p_grid[l - 1] < p_min)
      l--;
    else
      search_toa = false;
  }
  ARTS_USER_ERROR_IF (search_toa,
      "At least one atmospheric level with pressure larger p_min (=",
      p_min, ")\n"
      "is needed, but none is found.")
  const Index npn = l + 1;
  p_grid = p_grid[Range(0, npn)];

  const Index nsa = abs_species.nelem();

  // Check that there is at least one VMR species:
  ARTS_USER_ERROR_IF (nsa < 1,
    "There must be at least one absorption species.")

  // In contrast to older versions, we allow the data entries to be in arbitrary
  // order. that is, we have to look for all fields individually. we require the
  // abs_species related fields as well as the scat_species related fields to
  // have leading identifiers, namely 'abs_species' and 'scat_species'. it is
  // not mandatory that all fields available in atm_field_compact are applied
  // through abs_species and scat_species; left over fields are silently
  // ignored.
  // For temperature and altitude, occurence of exactly one field entry is
  // ensured. For other fields, the first match is used.
  // FIXME: or should a match be removed such that a second species occurence
  // would match a later-in-line field?

  bool found;
  const String as_type = "abs_species";
  const String ss_type = "scat_species";

  // Find temperature field:
  found = false;
  t_field.resize(npn, nlat, nlon);
  for (Index i = 0; i < nf; ++i) {
    if (c.get_string_grid(GFIELD4_FIELD_NAMES)[i] == "T") {
      ARTS_USER_ERROR_IF (found,
          "Only one temperature ('T') field allowed,\n"
          "but found at least 2.")
      found = true;
      t_field = c.data(i, Range(0, npn), Range(joker), Range(joker));
    }
  }
  ARTS_USER_ERROR_IF (!found,
    "One temperature ('T') field required, but none found")

  // Find Altitude field:
  found = false;
  z_field.resize(npn, nlat, nlon);
  for (Index i = 0; i < nf; ++i) {
    if (c.get_string_grid(GFIELD4_FIELD_NAMES)[i] == "z") {
      ARTS_USER_ERROR_IF (found,
          "Only one altitude ('z') field allowed,\n"
          "but found at least 2.")
      found = true;
      z_field = c.data(i, Range(0, npn), Range(joker), Range(joker));
    }
  }
  ARTS_USER_ERROR_IF (!found,
    "One altitude ('z') field required, but none found")

  // Extracting the required abs_species fields:

  vmr_field.resize(nsa, npn, nlat, nlon);
  for (Index j = 0; j < nsa; ++j) {
    const String as_name = Species::toShortName(abs_species[j][0].Spec());
    found = false;
    Index i = 0;
    String species_type;
    String species_name;
    while (!found && i < nf) {
      parse_atmcompact_speciestype(
          species_type, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);
      // do we have an abs_species type field?
      if (species_type == as_type) {
        parse_atmcompact_speciesname(
            species_name, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);
        if (species_name == as_name) {
          found = true;
          vmr_field(j, Range(joker), Range(joker), Range(joker)) =
              c.data(i, Range(0, npn), Range(joker), Range(joker));
        }
      }
      i++;
    }
    ARTS_USER_ERROR_IF (!found,
      "No field for absorption species '", as_name, "' found.")
  }

  //get the number of scat_species entries
  std::vector<Index> Idx;

  String species_type;
  for (Index i = 0; i < nf; ++i) {
    parse_atmcompact_speciestype(
        species_type, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);

    if (species_type == ss_type) {
      Idx.push_back(i);
    }
  }

  const Index nsp = Idx.size();

  // Extracting the required scattering species fields:
  particle_bulkprop_field.resize(nsp, npn, nlat, nlon);
  particle_bulkprop_field = NAN;
  particle_bulkprop_names.resize(nsp);

  // put scat_species entries in particle_bulkprop_field

  for (Index j = 0; j < nsp; ++j) {
    String species_name;
    String scat_type;

    parse_atmcompact_scattype(
        scat_type, c.get_string_grid(GFIELD4_FIELD_NAMES)[Idx[j]], delim);

    parse_atmcompact_speciesname(
        species_name, c.get_string_grid(GFIELD4_FIELD_NAMES)[Idx[j]], delim);

    particle_bulkprop_field(j, Range(joker), Range(joker), Range(joker)) =
        c.data(Idx[j], Range(0, npn), Range(joker), Range(joker));

    particle_bulkprop_names[j] = species_name + delim + scat_type;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void z_surfaceFromFileAndGrid(Matrix& z_surface,
                              const Vector& lat_grid,
                              const Vector& lon_grid,
                              const String& filename,
                              const Index& interp_order,
                              const Index& set_lowest_altitude_to_zero,
                              const Verbosity& verbosity) {
  CREATE_OUT3;

  out3 << "Reading GriddedField2 surface altitude from " << filename << "\n";
  GriddedField2 z_surface_field;
  xml_read_from_file(filename, z_surface_field, verbosity);

  out3 << "Surface altitude field interpolated back to lat_grid and lon_grid\n";
  GriddedFieldLatLonRegrid(z_surface_field,
                           lat_grid,
                           lon_grid,
                           z_surface_field,
                           interp_order,
                           verbosity);
  z_surface = z_surface_field.data;
  if (set_lowest_altitude_to_zero) {
    z_surface -= min(z_surface);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void z_surfaceConstantAltitude(Matrix& z_surface,
                               const Vector& lat_grid,
                               const Vector& lon_grid,
                               const Numeric& altitude,
                               const Verbosity& verbosity) {
  CREATE_OUT3;
  out3 << "Setting surface to constant altitude of " << altitude << " m\n";
  z_surface = Matrix(lat_grid.nelem() ? lat_grid.nelem() : 1, lon_grid.nelem() ? lon_grid.nelem() : 1, altitude);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void InterpAtmFieldToPosition(Numeric& outvalue,
                              const Vector& p_grid,
                              const Vector& lat_grid,
                              const Vector& lon_grid,
                              const Tensor3& z_field,
                              const Vector& rtp_pos,
                              const Tensor3& field,
                              const Verbosity& verbosity) {
  // Input checks
  chk_atm_grids(3, p_grid, lat_grid, lon_grid);
  chk_atm_field("input argument *field*",
                field,
                3,
                p_grid,
                lat_grid,
                lon_grid);
  chk_rte_pos(3, rtp_pos);

  // Determine grid positions
  GridPos gp_p, gp_lat, gp_lon;
  rte_pos2gridpos(gp_p,
                  gp_lat,
                  gp_lon,
                  3,
                  p_grid,
                  lat_grid,
                  lon_grid,
                  z_field,
                  rtp_pos);

  // Interpolate
  outvalue = interp_atmfield_by_gp(3, field, gp_p, gp_lat, gp_lon);

  CREATE_OUT3;
  out3 << "    Result = " << outvalue << "\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridDensify(  // WS Output:
    Vector& p_grid,
    Index& atmfields_checked,
    Index& atmgeom_checked,
    Index& cloudbox_checked,
    // WS Input:
    const Vector& p_grid_old,
    // Control Parameters:
    const Index& nfill,
    const Verbosity& verbosity) {
  // Check that p_grid and p_grid_old are not the same variable (pointing to the
  // same memory space). this as p_grid will be overwritten, but we will need
  // both data later on for data regridding.
  ARTS_USER_ERROR_IF (&p_grid == &p_grid_old,
      "The old and new grids (p_grid and p_grid_old) are not allowed\n"
      "to be identical (pointing to same memory space).\n"
      "But they are doing in your case.")

  // as we manipoulate the overall vertical grid (but not simultaneously the
  // atmospheric fields), we reset all atmfields related checked WSV to
  // unchecked, forcing the user to do the checks again.
  atmfields_checked = 0;
  atmgeom_checked = 0;
  cloudbox_checked = 0;

  // Check the keyword argument:
  ARTS_USER_ERROR_IF (nfill < 0, "Argument *nfill* must be >= 0.");

  // Nothing to do if nfill=0
  if (nfill > 0) {
    // Allocate new size for p_grid
    const Index n0 = p_grid_old.nelem();
    p_grid.resize((n0 - 1) * (1 + nfill) + 1);

    Index iout = 0;
    p_grid[0] = p_grid_old[0];

    for (Index i = 1; i < n0; i++) {
      Vector pnew;
      VectorNLogSpace(
          pnew, 2 + nfill, p_grid_old[i - 1], p_grid_old[i], verbosity);
      for (Index j = 1; j < nfill + 2; j++) {
        iout += 1;
        p_grid[iout] = pnew[j];
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridRefine(  // WS Output:
    Vector& p_grid,
    Index& atmfields_checked,
    Index& atmgeom_checked,
    Index& cloudbox_checked,
    // WS Input:
    const Vector& p_grid_old,
    // Control Parameters:
    const Numeric& p_step10,
    const Verbosity&) {
  // Check that p_grid and p_grid_old are not the same variable (pointing to the
  // same memory space). this as p_grid will be overwritten, but we will need
  // both data later on for data regridding.
  ARTS_USER_ERROR_IF (&p_grid == &p_grid_old,
      "The old and new grids (p_grid and p_grid_old) are not allowed\n"
      "to be identical (pointing to same memory space).\n"
      "But they are doing in your case.")

  // as we manipoulate the overall vertical grid (but not simultaneously the
  // atmospheric fields), we reset all atmfields related checked WSV to
  // unchecked, forcing the user to do the checks again.
  atmfields_checked = 0;
  atmgeom_checked = 0;
  cloudbox_checked = 0;

  // Check the keyword argument:
  ARTS_USER_ERROR_IF (p_step10 <= 0,
    "The keyword argument p_step must be >0.")

  // For consistency with other code around arts (e.g., correlation
  // lengths in atmlab), p_step is given as log10(p[Pa]). However, we
  // convert it here to the natural log:
  const Numeric p_step = log(pow(10.0, p_step10));

  // Now starting modification of p_grid

  // We will need the log of the pressure grid:
  Vector log_p_old(p_grid_old.nelem());
  transform(log_p_old, log, p_grid_old);

  //  const Numeric epsilon = 0.01 * p_step; // This is the epsilon that
  //                                         // we use for comparing p grid spacings.

  // Construct p_grid (new)
  // ----------------------

  ArrayOfNumeric log_p_new;  // We take log_p_new as an array of Numeric, so
                             // that we can easily build it up by appending new
                             // elements to the end.

  // Check whether there are pressure levels that are further apart
  // (in log(p)) than p_step, and insert additional levels if
  // necessary:

  log_p_new.push_back(log_p_old[0]);

  for (Index i = 1; i < log_p_old.nelem(); ++i) {
    const Numeric dp =
        log_p_old[i - 1] - log_p_old[i];  // The grid is descending.

    const Numeric dp_by_p_step = dp / p_step;
    //          cout << "dp_by_p_step: " << dp_by_p_step << "\n";

    // How many times does p_step fit into dp?
    const Index n = (Index)ceil(dp_by_p_step);
    // n is the number of intervals that we want to have in the
    // new grid. The number of additional points to insert is
    // n-1. But we have to insert the original point as well.
    //          cout << n << "\n";

    const Numeric ddp = dp / (Numeric)n;
    //          cout << "ddp: " << ddp << "\n";

    for (Index j = 1; j <= n; ++j)
      log_p_new.push_back(log_p_old[i - 1] - (Numeric)j * ddp);
  }

  // Copy ArrayOfNumeric to proper vector.
  Vector log_p(log_p_new.nelem());
  for (Index i = 0; i < log_p_new.nelem(); ++i) log_p[i] = log_p_new[i];

  // Copy the new grid to abs_p, removing the log:
  p_grid.resize(log_p.nelem());
  transform(p_grid, exp, log_p);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridFromZRaw(  //WS Output
    Vector& p_grid,
    //WS Input
    const GriddedField3& z_field_raw,
    const Index& no_negZ,
    const Verbosity&) {
  // original version excludes negative z. not clear, why this is. maybe is
  // currently a convention somehwere in ARTS (DOIT?). negative z seem, however,
  // to work fine for clear-sky cases. so we make the negative z exclude an
  // option (for consistency until unclarities solved, default: do exclude)
  Vector p_grid_raw = z_field_raw.get_numeric_grid(GFIELD3_P_GRID);

  Index i;
  if (is_increasing(z_field_raw.data(joker, 0, 0))) {
    i = 0;
    if (no_negZ) {
      while (z_field_raw.data(i, 0, 0) < 0.0) i++;
    }
    p_grid = p_grid_raw[Range(i, joker)];
  } else if (is_decreasing(z_field_raw.data(joker, 0, 0))) {
    i = z_field_raw.data.npages() - 1;
    if (no_negZ) {
      while (z_field_raw.data(i, 0, 0) < 0.0) i--;
    }
    p_grid = reverse(p_grid_raw);
  } else {
    ARTS_USER_ERROR ("z_field_raw needs to be monotonous, but this is not the case.\n")
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void lat_gridFromZRaw(  //WS Output
    Vector& lat_grid,
    //WS Input
    const GriddedField3& z_field_raw,
    const Verbosity&) {
  lat_grid = z_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void lon_gridFromZRaw(  //WS Output
    Vector& lon_grid,
    //WS Input
    const GriddedField3& z_field_raw,
    const Verbosity&) {
  lon_grid = z_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atm_gridsFromZRaw(  //WS Output
    Vector& p_grid,
    Vector& lat_grid,
    Vector& lon_grid,
    //WS Input
    const GriddedField3& z_field_raw,
    const Index& no_negZ,
    const Verbosity& v) {
  p_gridFromZRaw(p_grid, z_field_raw, no_negZ, v);
  lat_gridFromZRaw(lat_grid, z_field_raw, v);
  lon_gridFromZRaw(lon_grid, z_field_raw, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void lat_gridFromRawField(  //WS Output
    Vector& lat_grid,
    //WS Input
    const GriddedField3& field_raw,
    const Verbosity&) {
  lat_grid = field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void lon_gridFromRawField(  //WS Output
    Vector& lon_grid,
    //WS Input
    const GriddedField3& field_raw,
    const Verbosity&) {
  lon_grid = field_raw.get_numeric_grid(GFIELD3_LON_GRID);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void wind_u_fieldIncludePlanetRotation(Tensor3& wind_u_field,
                                       const Vector& p_grid,
                                       const Vector& lat_grid,
                                       const Vector& lon_grid,
                                       const Vector& refellipsoid,
                                       const Tensor3& z_field,
                                       const Numeric& planet_rotation_period,
                                       const Verbosity&) {
  const Index np = p_grid.nelem();
  const Index na = lat_grid.nelem();
  const Index no = lon_grid.nelem();

  chk_atm_field("z_field", z_field, 3, p_grid, lat_grid, lon_grid);
  if (wind_u_field.npages() > 0) {
    chk_atm_field("wind_u_field",
                  wind_u_field,
                  3,
                  p_grid,
                  lat_grid,
                  lon_grid);
  } else {
    wind_u_field.resize(np, na, no);
    wind_u_field = 0.;
  }

  const Numeric k1 = 2 * Constant::pi / planet_rotation_period;

  for (Index a = 0; a < na; a++) {
    const Numeric k2 = k1 * Conversion::cosd(lat_grid[a]);
    const Numeric re = refell2r(refellipsoid, lat_grid[a]);

    for (Index o = 0; o < no; o++) {
      for (Index p = 0; p < np; p++) {
        wind_u_field(p, a, o) += k2 * (re + z_field(p, a, o));
      }
    }
  }
}

// A small help function
void z2g(Numeric& g, const Numeric& r, const Numeric& g0, const Numeric& z) {
  const Numeric x = r / (r + z);
  g = g0 * x * x;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void z_fieldFromHSE(Workspace& ws,
                    Tensor3& z_field,
                    const Vector& p_grid,
                    const Vector& lat_grid,
                    const Vector& lon_grid,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const ArrayOfArrayOfSpeciesTag& abs_species,
                    const Tensor3& t_field,
                    const Tensor4& vmr_field,
                    const Vector& refellipsoid,
                    const Matrix& z_surface,
                    const Index& atmfields_checked,
                    const Agenda& g0_agenda,
                    const Numeric& molarmass_dry_air,
                    const Numeric& p_hse,
                    const Numeric& z_hse_accuracy,
                    const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");

  // Some general variables
  const Index np = p_grid.nelem();
  const Index nlat = t_field.nrows();
  const Index nlon = t_field.ncols();
  //
  const Index firstH2O = find_first_species(
      abs_species, Species::fromShortName("H2O"));

  if (firstH2O < 0) {
    CREATE_OUT1;
    out1 << "No water vapour tag group in *abs_species*.\n"
         << "Be aware that this leads to significant variations in atmospheres\n"
         << "that contain considerable amounts of water vapour (e.g. Earth)!\n";
  }
  //
  ARTS_USER_ERROR_IF (p_hse > p_grid[0] || p_hse < p_grid[np - 1],
      "The value of *p_hse* must be inside the range of *p_grid*:"
      "  p_hse  = ", p_hse, " Pa\n"
      "  p_grid = ", p_grid[np-1],
      " - ", p_grid[0], " Pa\n")
  //
  ARTS_USER_ERROR_IF (z_hse_accuracy <= 0,
                      "The value of *z_hse_accuracy* must be > 0.");
  //
  chk_latlon_true(3, lat_grid, lat_true, lon_true);

  // Determine interpolation weights for p_hse
  //
  ArrayOfGridPos gp(1);
  Matrix itw(1, 2);
  p2gridpos(gp, p_grid, Vector(1, p_hse));
  interpweights(itw, gp);

  // // Molecular weight of water vapour
  const Numeric mw = 18.016;

  // mw/molarmass_dry_air matches eps in Eq. 3.14 in Wallace&Hobbs:
  const Numeric k = 1 - mw / molarmass_dry_air;

  // Gas constant for 1kg dry air:
  const Numeric rd = 1e3 * GAS_CONSTANT / molarmass_dry_air;

  // The calculations
  //
  for (Index ilat = 0; ilat < nlat; ilat++) {
    // The reference ellipsoid is already adjusted to internal 1D or 2D
    // views, and lat_grid is the relevant grid for *refellipsoid*, also
    // for 2D. On the other hand, extraction of g0 requires that the true
    // position is determined.

    // Radius of reference ellipsoid
    Numeric re;
    if (3 == 1) {
      re = refellipsoid[0];
    } else {
      re = refell2r(refellipsoid, lat_grid[ilat]);
    }

    for (Index ilon = 0; ilon < nlon; ilon++) {
      // Determine true latitude and longitude
      Numeric lat, lon;
      Vector pos(3);  // pos[0] can be a dummy value
      if (3 >= 2) {
        pos[1] = lat_grid[ilat];
        if (3 == 3) {
          pos[2] = lon_grid[ilon];
        }
      }
      pos2true_latlon(
          lat, lon, 3, lat_grid, lat_true, lon_true, pos);

      // Get g0
      Numeric g0;
      g0_agendaExecute(ws, g0, lat, lon, g0_agenda);

      // Determine altitude for p_hse
      Vector z_hse(1);
      interp(z_hse, itw, z_field(joker, ilat, ilon), gp);

      Numeric z_acc = 2 * z_hse_accuracy;

      while (z_acc > z_hse_accuracy) {
        z_acc = 0;
        Numeric g2 = g0;

        for (Index ip = 0; ip < np - 1; ip++) {
          // Calculate average g
          if (ip == 0) {
            z2g(g2, re, g0, z_field(ip, ilat, ilon));
          }
          const Numeric g1 = g2;
          z2g(g2, re, g0, z_field(ip + 1, ilat, ilon));
          //
          const Numeric g = (g1 + g2) / 2.0;

          //Average water vapour VMR
          Numeric hm;
          if (firstH2O < 0) {
            hm = 0.0;
          } else {
            hm = 0.5 * (vmr_field(firstH2O, ip, ilat, ilon) +
                        vmr_field(firstH2O, ip + 1, ilat, ilon));
          }

          // Average virtual temperature (no liquid water)
          // (cf. e.g. Eq. 3.16 in Wallace&Hobbs)
          const Numeric tv =
              (1 / (2 * (1 - hm * k))) *
              (t_field(ip, ilat, ilon) + t_field(ip + 1, ilat, ilon));

          // Change in vertical altitude from ip to ip+1
          //  (cf. e.g. Eq. 3.24 in Wallace&Hobbs)
          const Numeric dz = rd * (tv / g) * log(p_grid[ip] / p_grid[ip + 1]);

          // New altitude
          Numeric znew = z_field(ip, ilat, ilon) + dz;
          z_acc = max(z_acc, fabs(znew - z_field(ip + 1, ilat, ilon)));
          z_field(ip + 1, ilat, ilon) = znew;
        }

        // Adjust to z_hse
        Vector z_tmp(1);
        interp(z_tmp, itw, z_field(joker, ilat, ilon), gp);
        z_field(joker, ilat, ilon) -= z_tmp[0] - z_hse[0];
      }
    }
  }

  // Check that there is no gap between the surface and lowest pressure
  // level
  // (This code is copied from *basics_checkedCalc*. Make this to an internal
  // function if used in more places.)
  for (Index row = 0; row < z_surface.nrows(); row++) {
    for (Index col = 0; col < z_surface.ncols(); col++) {
      ARTS_USER_ERROR_IF (z_surface(row, col) < z_field(0, row, col) ||
                          z_surface(row, col) >= z_field(z_field.npages() - 1, row, col),
          "The surface altitude (*z_surface*) cannot be outside "
          "of the altitudes in *z_field*.")
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldSetConstant(Tensor4& vmr_field,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const String& species,
                          const Numeric& vmr_value,
                          const Verbosity&) {
  // Check input
  chk_if_in_range("vmr_value", vmr_value, 0, 1);
  //
  ARTS_USER_ERROR_IF (abs_species.nelem() != vmr_field.nbooks(),
        "Size of *vmr_field* and length of *abs_species* do not agree.");

  // Find position for this species.
  const ArrayOfSpeciesTag tag(species);
  Index si = chk_contains("species", abs_species, tag);

  // Apply value
  vmr_field(si, joker, joker, joker) = vmr_value;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldSetAllConstant(Tensor4& vmr_field,
                             const ArrayOfArrayOfSpeciesTag& abs_species,
                             const Vector& vmr_values,
                             const Verbosity& verbosity) {
  CREATE_OUT3;

  const Index nspecies = abs_species.nelem();

  ARTS_USER_ERROR_IF (vmr_values.nelem() not_eq nspecies,
                      "Not same number of vmr_values as abs_species.");

  out3 << "Setting all " << nspecies << " species to constant VMR\n";

  for (Index i = 0; i < nspecies; i++) {
    const ArrayOfSpeciesTag& a_abs_species = abs_species[i];
    const String species_tag_name = a_abs_species.Name();
    vmr_fieldSetConstant(
        vmr_field, abs_species, species_tag_name, vmr_values[i], verbosity);
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldSetRh(Workspace& ws,
                    Tensor4& vmr_field,
                    const ArrayOfArrayOfSpeciesTag& abs_species,
                    const Tensor3& t_field,
                    const Vector& p_grid,
                    const Agenda& water_p_eq_agenda,
                    const Numeric& rh,
                    const Numeric& vmr_threshold,
                    const Verbosity&)
{
  // Locate H2O
  const Index ih2o = find_first_species(abs_species, Species::fromShortName("H2O"));
  ARTS_USER_ERROR_IF (ih2o < 0, "There is no H2O species in *abs_species*.")

  // Calculate partial pressure matching selected RH
  Tensor3 p_rh;
  water_p_eq_agendaExecute(ws, p_rh, t_field, water_p_eq_agenda);
  p_rh *= rh;

  // Put in VMR
  for (Index p=0; p<vmr_field.npages(); ++p) {
    for (Index a=0; a<vmr_field.nrows(); ++a) {
      for (Index o=0; o<vmr_field.ncols(); ++o) {
        if (vmr_field(ih2o, p, a, o) > vmr_threshold) {
          vmr_field(ih2o, p, a, o) = p_rh(p, a, o) / p_grid[p];
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atm_fieldLteExternalPartitionFunction(
    Index& nlte_do,
    AtmField& atm_field,
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& nlte_quantum_identifiers,
    const Verbosity& verbosity) {
  using Constant::h;

  ARTS_USER_ERROR_IF(not atm_field.regularized, "Only for regular atmospheric field")
  ARTS_USER_ERROR_IF(not atm_field.has(Atm::Key::t), "Atmospheric field must have temperature field")
  const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();

  CREATE_OUT2;
  const Index nn = nlte_quantum_identifiers.nelem(), np = t_field.npages(),
              nlat = t_field.nrows(), nlon = t_field.ncols();
  if (nn == 0) return;

  Tensor4 nlte_tensor4(nn, np, nlat, nlon);
  nlte_do = 1;
  ArrayOfIndex checked(nn, 0);

  for (Index in = 0; in < nn; in++) {
    const QuantumIdentifier& qi = nlte_quantum_identifiers[in];
    Tensor3View lte = nlte_tensor4(in, joker, joker, joker);

    for (auto& abs_lines : abs_lines_per_species) {
      for (auto& band : abs_lines) {
        for (Index k=0; k<band.NumLines(); k++) {
          const Quantum::Number::StateMatch lt(qi, band.lines[k].localquanta, band.quantumidentity);
          if (lt == Quantum::Number::StateMatchType::Level and lt.low) {
            band.population = Absorption::PopulationType::NLTE;
            
            if (not checked[in]) {
              checked[in] = 1;
              
              for (Index ip = 0; ip < np; ip++) {
                for (Index ilat = 0; ilat < nlat; ilat++) {
                  for (Index ilon = 0; ilon < nlon; ilon++) {
                    lte(ip, ilat, ilon) =
                    boltzman_factor(t_field(ip, ilat, ilon), band.lines[k].E0) *
                    band.lines[k].glow / single_partition_function(
                      t_field(ip, ilat, ilon), band.Isotopologue());
                  }
                }
              }
            }
          }
          
          if (lt == Quantum::Number::StateMatchType::Level and lt.upp) {
            band.population = Absorption::PopulationType::NLTE;
            
            if (not checked[in]) {
              checked[in] = 1;
              for (Index ip = 0; ip < np; ip++) {
                for (Index ilat = 0; ilat < nlat; ilat++) {
                  for (Index ilon = 0; ilon < nlon; ilon++) {
                    lte(ip, ilat, ilon) =
                    boltzman_factor(t_field(ip, ilat, ilon), band.lines[k].E0 + h*band.lines[k].F0) *
                    band.lines[k].gupp / single_partition_function(
                      t_field(ip, ilat, ilon), band.Isotopologue());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  for (Index in = 0; in < nn; in++) {
    if (not checked[in]) {
      out2 << "Did not find match among lines for: "
           << nlte_quantum_identifiers[in] << "\n";
    }
    atm_field.nlte[nlte_quantum_identifiers[in]] = Tensor3{nlte_tensor4[in]};
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atm_fieldLteInternalPartitionFunction(
    Index& nlte_do,
    AtmField& atm_field,
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& nlte_quantum_identifiers,
    const Verbosity& verbosity) {
  using Constant::h;

  ARTS_USER_ERROR_IF(not atm_field.regularized, "Only for regular atmospheric field")
  ARTS_USER_ERROR_IF(not atm_field.has(Atm::Key::t), "Atmospheric field must have temperature field")
  const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();

  CREATE_OUT2;
  const Index nn = nlte_quantum_identifiers.nelem(), np = t_field.npages(),
              nlat = t_field.nrows(), nlon = t_field.ncols();
  if (nn == 0) return;

  // Find where they are positioned and how many different molecules there are for the NLTE fields
  ArrayOfIndex part_fun_pos(nn, 0);
  Index x = 1;
  for (Index in = 1; in < nn; in++) {
    bool found = false;
    for (Index ix = 0; ix < in; ix++) {
      if (nlte_quantum_identifiers[in].Species() ==
              nlte_quantum_identifiers[ix].Species() and
          nlte_quantum_identifiers[in].Isotopologue() ==
              nlte_quantum_identifiers[ix].Isotopologue()) {
        part_fun_pos[in] = part_fun_pos[ix];
        found = true;
        break;
      }
    }
    if (not found) {
      part_fun_pos[in] = x;
      x++;
    }
  }
  
  Tensor4 part_fun(x, np, nlat, nlon, 0.0);
  Tensor4 nlte_tensor4(nn, np, nlat, nlon, 0);
  nlte_do = 1;
  ArrayOfIndex checked(nn, 0);

  for (Index in = 0; in < nn; in++) {
    const QuantumIdentifier& qi = nlte_quantum_identifiers[in];
    Tensor3View lte = nlte_tensor4(in, joker, joker, joker);

    for (auto& abs_lines : abs_lines_per_species) {
      for (auto& band : abs_lines) {
        for (Index k=0; k<band.NumLines(); k++) {
          const Quantum::Number::StateMatch lt(qi, band.lines[k].localquanta, band.quantumidentity);
          if (lt == Quantum::Number::StateMatchType::Level and lt.low) {
            band.population = Absorption::PopulationType::NLTE;
            
            if (not checked[in]) {
              checked[in] = 1;
              
              for (Index ip = 0; ip < np; ip++) {
                for (Index ilat = 0; ilat < nlat; ilat++) {
                  for (Index ilon = 0; ilon < nlon; ilon++) {
                    lte(ip, ilat, ilon) =
                      boltzman_factor(t_field(ip, ilat, ilon), band.lines[k].E0) *
                                      band.lines[k].glow;
                    part_fun(part_fun_pos[in], ip, ilat, ilon) +=
                      lte(ip, ilat, ilon);
                  }
                }
              }
            }
          }
          
          if (lt == Quantum::Number::StateMatchType::Level and lt.upp) {
            band.population = Absorption::PopulationType::NLTE;
            
            if (not checked[in]) {
              checked[in] = 1;
              
              for (Index ip = 0; ip < np; ip++) {
                for (Index ilat = 0; ilat < nlat; ilat++) {
                  for (Index ilon = 0; ilon < nlon; ilon++) {
                    lte(ip, ilat, ilon) =
                      boltzman_factor(t_field(ip, ilat, ilon),
                                      band.lines[k].E0 + h * band.lines[k].F0) *
                                      band.lines[k].gupp;
                    part_fun(part_fun_pos[in], ip, ilat, ilon) +=
                      lte(ip, ilat, ilon);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  for (Index in = 0; in < nn; in++) {
    if (not checked[in]) {
      out2 << "Did not find match among lines for: "
           << nlte_quantum_identifiers[in] << "\n";
    }
  }
  
  for (Index in = 0; in < nn; in++) {
    if (checked[in]) {
      nlte_tensor4(in, joker, joker, joker) /=
        part_fun(part_fun_pos[in], joker, joker, joker);
    }
    atm_field.nlte[nlte_quantum_identifiers[in]] = Tensor3{nlte_tensor4[in]};
  }
}
