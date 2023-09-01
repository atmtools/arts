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
#include "atm.h"
#include "debug.h"
#include "species_tags.h"
#include "absorption.h"
#include <workspace.h>
#include "check_input.h"
#include "cloudbox.h"
#include "geodetic.h"
#include "gridded_fields.h"
#include "igrf13.h"
#include "interpolation.h"
#include "interp.h"
#include "linescaling.h"
#include "matpack_data.h"
#include "rte.h"
#include "special_interp.h"
#include "surf.h"
#include "xml_io.h"
#include "arts_omp.h"

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
                              const Index& prepend) {
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
 \param levels
 */
void AtmFieldPRegridHelper(Index& ing_min,
                           Index& ing_max,
                           ArrayOfLagrangeLogInterpolation& lag_p,
                           Matrix& itw,
                           ConstVectorView p_grid_out,
                           ConstVectorView p_grid_in,
                           const Index& interp_order) {
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
    const Index& interp_order) {
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
                        interp_order);

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
    const Index& interp_order) {
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
                        interp_order);

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
    const GriddedField2& gfraw_in) {
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
    const GriddedField3& gfraw_in) {
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
    const GriddedField4& gfraw_in) {
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
    const ArrayOfGriddedField3& gfraw_in) {
  if (!gfraw_in.nelem()) {
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
    const GriddedField2& gfraw_in_orig) {
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
    const GriddedField3& gfraw_in_orig) {
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
    const GriddedField4& gfraw_in_orig) {
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
    const ArrayOfGriddedField3& gfraw_in) {
  gfraw_out.resize(gfraw_in.nelem());

  for (Index i = 0; i < gfraw_in.nelem(); i++)
    GriddedFieldLatLonExpand(gfraw_out[i], gfraw_in[i]);
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
 \param levels
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
                               const Index& zeropadding) {
  chk_griddedfield_gridname(gfraw_in, p_grid_index, "Pressure");

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
    const Index& zeropadding) {
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
                            zeropadding);

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
    const Index& zeropadding) {
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
                            zeropadding);

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
    const Index& zeropadding) {
  agfraw_out.resize(agfraw_in.nelem());

  for (Index i = 0; i < agfraw_in.nelem(); i++) {
    GriddedFieldPRegrid(agfraw_out[i],
                        p_grid,
                        agfraw_in[i],
                        interp_order,
                        zeropadding);
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
 \param levels
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
                                    const Index& interp_order) {
  ARTS_USER_ERROR_IF (!lat_true.nelem(),
                      "The new latitude grid is not allowed to be empty.");
  ARTS_USER_ERROR_IF (!lon_true.nelem(),
                      "The new longitude grid is not allowed to be empty.");

  chk_griddedfield_gridname(gfraw_in, lat_grid_index, "Latitude");
  chk_griddedfield_gridname(gfraw_in, lon_grid_index, "Longitude");
  ARTS_USER_ERROR_IF (gfraw_in.get_grid_size(lat_grid_index) == 1 ||
      gfraw_in.get_grid_size(lon_grid_index) == 1,
                      "Raw data has to be true 3D data (nlat>1 and nlon>1).");

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
    const Index& interp_order) {
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
                                 interp_order);

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
    const Index& interp_order) {
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
                                 interp_order);

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
    const Index& interp_order) {
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
                                 interp_order);

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
    const Index& interp_order) {
  agfraw_out.resize(agfraw_in.nelem());

  for (Index i = 0; i < agfraw_in.nelem(); i++) {
    GriddedFieldLatLonRegrid(agfraw_out[i],
                             lat_true,
                             lon_true,
                             agfraw_in[i],
                             interp_order);
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
    const ArrayOfString& field_names) {
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
    const ArrayOfString& condensibles) {
  Index nf;  // Will hold new size

  // Add book
  atm_fields_compactExpand(af, nf, name, prepend);

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
    const Index& prepend) {
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
      atm_fields_compact, new_n_fields, name, prepend);

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
    const Numeric& threshold) {
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
    const GriddedField3& field) {
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
    const ArrayOfString& condensibles) {
  for (Index i = 0; i < batch_atm_fields_compact.nelem(); i++) {
    atm_fields_compactAddConstant(batch_atm_fields_compact[i],
                                  name,
                                  value,
                                  prepend,
                                  condensibles);
  }
}

// Workspace method, doxygen header is auto-generated
// 2011-05-09 Gerrit Holl
void batch_atm_fields_compactAddSpecies(  // WS Output:
    ArrayOfGriddedField4& batch_atm_fields_compact,
    // WS Generic Input:
    const String& name,
    const GriddedField3& species,
    const Index& prepend) {
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
          batch_atm_fields_compact[i], name, species, prepend);
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
    const Numeric& threshold) {
  for (Index i = 0; i < batch_atm_fields_compact.nelem(); i++) {
    atm_fields_compactCleanup(
        batch_atm_fields_compact[i], threshold);
  }
}

// Workspace method, doxygen header is auto-generated.
void batch_atm_fields_compactFromArrayOfMatrix(  // WS Output:
    ArrayOfGriddedField4& batch_atm_fields_compact,
    // WS Input:
    // WS Generic Input:
    const ArrayOfMatrix& am,
    // Control Parameters:
    const ArrayOfString& field_names) {
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
                                   field_names);
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
    AtmField& atm_field,
    ArrayOfString& particle_bulkprop_names,

    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const GriddedField4& atm_fields_compact,
    const String& delim,
    // Control parameters:
    const Index& check_gridnames) {
  // Make a handle on atm_fields_compact to save typing:
  const GriddedField4& c = atm_fields_compact;

    Vector z_grid;
    Vector lat_grid;
    Vector lon_grid;

  // Check if the grids in our data match 3
  // (throws an error if the dimensionality is not correct):
  chk_atm_grids(3,
                c.get_numeric_grid(GFIELD4_P_GRID),
                c.get_numeric_grid(GFIELD4_LAT_GRID),
                c.get_numeric_grid(GFIELD4_LON_GRID));

  // Optional check for gridnames.
  if (check_gridnames == 1) {
    chk_griddedfield_gridname(c, 1, "Altitude");
    chk_griddedfield_gridname(c, 2, "Latitude");
    chk_griddedfield_gridname(c, 3, "Longitude");
  }

  const Index nf = c.get_grid_size(GFIELD4_FIELD_NAMES);

  // Grids:
  z_grid = c.get_numeric_grid(GFIELD4_P_GRID);
  lat_grid = c.get_numeric_grid(GFIELD4_LAT_GRID);
  lon_grid = c.get_numeric_grid(GFIELD4_LON_GRID);

  const Index nsa = abs_species.nelem();

  // Check that there is at least one VMR species:
  ARTS_USER_ERROR_IF (nsa < 1,
    "There must be at least one absorption species.")
  
  // Set TOA
  atm_fieldInit(atm_field, max(z_grid));
  GriddedField3 field_data;
  field_data.set_grid_name(0, "Altitude");
  field_data.set_grid_name(1, "Latitude");
  field_data.set_grid_name(2, "Longitude");
  field_data.set_grid(0, z_grid);
  field_data.set_grid(1, lat_grid);
  field_data.set_grid(2, lon_grid);

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
  for (Index i = 0; i < nf; ++i) {
    if (c.get_string_grid(GFIELD4_FIELD_NAMES)[i] == "T") {
      ARTS_USER_ERROR_IF (found,
          "Only one temperature ('T') field allowed,\n"
          "but found at least 2.")
      found = true;
      field_data.data = Tensor3{c.data[i]};
      atm_field[Atm::Key::t] = field_data;
    }
  }
  ARTS_USER_ERROR_IF (!found,
    "One temperature ('T') field required, but none found")

  // Find Pressure field:
  found = false;
  for (Index i = 0; i < nf; ++i) {
    if (c.get_string_grid(GFIELD4_FIELD_NAMES)[i] == "p") {
      ARTS_USER_ERROR_IF (found,
          "Only one pressure ('p') field allowed,\n"
          "but found at least 2.")
      found = true;
      field_data.data = Tensor3{c.data[i]};
      atm_field[Atm::Key::p] = field_data;
    }
  }
  ARTS_USER_ERROR_IF (!found,
    "One pressure ('p') field required, but none found")

  // Extracting the required abs_species fields:
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
          field_data.data = Tensor3{c.data[i]};
          atm_field[abs_species[j]] = field_data;
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
  particle_bulkprop_names.resize(nsp);

  // put scat_species entries in particle_bulkprop_field

  for (Index j = 0; j < nsp; ++j) {
    String species_name;
    String scat_type;

    parse_atmcompact_scattype(
        scat_type, c.get_string_grid(GFIELD4_FIELD_NAMES)[Idx[j]], delim);

    parse_atmcompact_speciesname(
        species_name, c.get_string_grid(GFIELD4_FIELD_NAMES)[Idx[j]], delim);

    particle_bulkprop_names[j] = species_name + delim + scat_type;
    field_data.data = Tensor3{c.data[Idx[j]]};
    atm_field[ParticulatePropertyTag{particle_bulkprop_names[j]}] = field_data;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void InterpAtmFieldToPosition(AtmPoint& atm_point,
                              const AtmField& atm_field,
                              const Vector& rtp_pos) {
  atm_point = atm_field.at({rtp_pos[0]}, {rtp_pos[1]}, {rtp_pos[2]})[0];
}

// A small help function
void z2g(Numeric& g, const Numeric& r, const Numeric& g0, const Numeric& z) {
  const Numeric x = r / (r + z);
  g = g0 * x * x;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atm_fieldLteExternalPartitionFunction(
    Index& nlte_do,
    AtmField& atm_field,
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& nlte_quantum_identifiers) {
  using Constant::h;

  // FIXME: REQUIRES REGULAR GRIDS
  Vector z_grid, lat_grid, lon_grid;
  Tensor3 t_field;
  ARTS_USER_ERROR_IF(not atm_field.has(Atm::Key::t), "Atmospheric field must have temperature field")
  //const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();

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
    }
    //atm_field[nlte_quantum_identifiers[in]] = Tensor3{nlte_tensor4[in]};
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atm_fieldLteInternalPartitionFunction(
    Index& nlte_do,
    AtmField& atm_field,
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& nlte_quantum_identifiers) {
  using Constant::h;

  // FIXME: REQUIRES REGULAR GRIDS
  Vector z_grid, lat_grid, lon_grid;
  Tensor3 t_field;
  ARTS_USER_ERROR_IF(not atm_field.has(Atm::Key::t), "Atmospheric field must have temperature field")
  //const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();

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
    }
  }
  
  for (Index in = 0; in < nn; in++) {
    if (checked[in]) {
      nlte_tensor4(in, joker, joker, joker) /=
        part_fun(part_fun_pos[in], joker, joker, joker);
    }
    ARTS_USER_ERROR("ERROR: atm_fieldLteInternalPartitionFunction not implemented yet\n")
    //atm_field[nlte_quantum_identifiers[in]] = Tensor3{nlte_tensor4[in]};
  }
}
