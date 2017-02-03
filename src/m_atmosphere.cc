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
#include "absorption.h"
#include "abs_species_tags.h"
#include "agenda_class.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "geodetic.h"
#include "global_data.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "interpolation_poly.h"
#include "matpackIII.h"
#include "messages.h"
#include "rte.h"
#include "special_interp.h"
#include "xml_io.h"

extern const Index GFIELD3_P_GRID;
extern const Index GFIELD3_LAT_GRID;
extern const Index GFIELD3_LON_GRID;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_P_GRID;
extern const Index GFIELD4_LAT_GRID;
extern const Index GFIELD4_LON_GRID;

extern const Numeric DEG2RAD;
extern const Numeric PI;
extern const Numeric GAS_CONSTANT;

//! Data value accuracy requirement for values at 0 and 360 deg if longitudes are cyclic
/*!
 */
extern const Numeric EPSILON_LON_CYCLIC = 2*DBL_EPSILON;


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
                              const Verbosity&)
{
    // Number of fields already present:
    nf = af.get_string_grid(GFIELD4_FIELD_NAMES).nelem();

    if (prepend)
    {
        // Add name of new field to field name list:
        ArrayOfString& as = af.get_string_grid(GFIELD4_FIELD_NAMES);
        as.insert(as.begin(), name);
    }
    else
    {
        // Add name of new field to field name list:
        af.get_string_grid(GFIELD4_FIELD_NAMES).push_back(name);
    }

    // Save original fields:
    const Tensor4 dummy = af.data;

    // Adjust size:
    af.resize( nf+1, dummy.npages(), dummy.nrows(), dummy.ncols() );
    nf++; // set to new number of fields

    // Copy back original field:
    af.data( Range((prepend && nf > 1)?1:0, nf-1), joker, joker, joker ) = dummy;
}



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

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
                           ArrayOfGridPosPoly& gp_p,
                           Matrix& itw,
                           ConstVectorView p_grid_out,
                           ConstVectorView p_grid_in,
                           const Index& interp_order,
                           const Verbosity& verbosity)
{
    CREATE_OUT2;

    out2 << "  Interpolation order: " << interp_order << "\n";

    ing_min = 0;
    ing_max = p_grid_out.nelem()-1;
    chk_interpolation_pgrids("Atmospheric field to p_grid_out",
                             p_grid_in,  p_grid_out, interp_order);

    Index nelem_in_range = ing_max - ing_min + 1;

    // Calculate grid positions:
    if (nelem_in_range>0)
      {
        gp_p.resize(nelem_in_range);
        p2gridpos_poly(gp_p, p_grid_in, p_grid_out[Range(ing_min, nelem_in_range)],
                       interp_order);

        // Interpolation weights:
        itw.resize(nelem_in_range, interp_order+1);
        interpweights(itw, gp_p);
      }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldPRegrid(// WS Generic Output:
                     Tensor3& atmtensor_out,
                     // WS Generic Input:
                     const Tensor3& atmtensor_in_orig,
                     const Vector& p_grid_new,
                     const Vector& p_grid_old,
                     const Index& interp_order,
                     const Verbosity& verbosity)
{
    // check that p_grid_old is the p_grid associated with atmtensor_in_orig. we can
    // only check for consistent size with p_grid_old.
    if (atmtensor_in_orig.npages() != p_grid_old.nelem())
      {
        ostringstream os;
        os << "p_grid_old is supposed to be the p_grid associated with the "
           << "atmospheric field.\n"
           << "However, it is not as their sizes are inconsistent.\n"; 
        throw runtime_error(os.str());
      }

    const Tensor3* atmtensor_in_pnt;
    Tensor3 atmtensor_in_copy;

    if (&atmtensor_in_orig == &atmtensor_out)
    {
        atmtensor_in_copy = atmtensor_in_orig;
        atmtensor_in_pnt = &atmtensor_in_copy;
    }
    else
        atmtensor_in_pnt = &atmtensor_in_orig;

    const Tensor3 &atmtensor_in = *atmtensor_in_pnt;

    // Resize output tensor
    atmtensor_out.resize(p_grid_new.nelem(), atmtensor_in.nrows(), atmtensor_in.ncols());

    ArrayOfGridPosPoly gp_p;
    Matrix itw;

    Index ing_min, ing_max;

    AtmFieldPRegridHelper(ing_min, ing_max, gp_p, itw, 
                          p_grid_new, p_grid_old,
                          interp_order, verbosity);

    // Interpolate:
    if ( (ing_max - ing_min < 0) || (ing_max - ing_min + 1 != p_grid_new.nelem()) )
      {
        ostringstream os;
        os << "New grid seems not to be sufficiently covered by old grid.\n"; 
        throw runtime_error(os.str());
      }

    for (Index i = 0; i < atmtensor_in.nrows(); i++)
        for (Index j = 0; j < atmtensor_in.ncols(); j++)
            interp(atmtensor_out(joker, i, j),
                   itw, atmtensor_in(joker, i, j), gp_p);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldPRegrid(// WS Generic Output:
                     Tensor4& atmtensor_out,
                     // WS Generic Input:
                     const Tensor4& atmtensor_in_orig,
                     const Vector& p_grid_new,
                     const Vector& p_grid_old,
                     const Index& interp_order,
                     const Verbosity& verbosity)
{
    const Tensor4* atmtensor_in_pnt;
    Tensor4 atmtensor_in_copy;

    if (&atmtensor_in_orig == &atmtensor_out)
    {
        atmtensor_in_copy = atmtensor_in_orig;
        atmtensor_in_pnt = &atmtensor_in_copy;
    }
    else
        atmtensor_in_pnt = &atmtensor_in_orig;

    const Tensor4 &atmtensor_in = *atmtensor_in_pnt;

    // Resize output tensor
    atmtensor_out.resize(atmtensor_in.nbooks(), p_grid_new.nelem(),
                         atmtensor_in.nrows(), atmtensor_in.ncols());

    ArrayOfGridPosPoly gp_p;
    Matrix itw;

    Index ing_min, ing_max;

    AtmFieldPRegridHelper(ing_min, ing_max, gp_p, itw, 
                          p_grid_new, p_grid_old,
                          interp_order, verbosity);

    // Interpolate:
    if ( (ing_max - ing_min < 0) || (ing_max - ing_min + 1 != p_grid_new.nelem()) )
      {
        ostringstream os;
        os << "New grid seems not to be sufficiently covered by old grid.\n"; 
        throw runtime_error(os.str());
      }

    for (Index b = 0; b <  atmtensor_in.nbooks(); b++)
        for (Index i = 0; i <  atmtensor_in.nrows(); i++)
            for (Index j = 0; j <  atmtensor_in.ncols(); j++)
                interp(atmtensor_out(b, joker, i, j),
                       itw,  atmtensor_in(b, joker, i, j), gp_p);
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
                                            const GriddedField& gfield)
{
    chk_griddedfield_gridname(gfield, ilat, "Latitude");
    chk_griddedfield_gridname(gfield, ilon, "Longitude");

    if (lon_grid.empty())
    {
        chk_size("gfield.lon_grid", gfield.get_numeric_grid(ilon), 1);

        if (lat_grid.empty())
            chk_size("gfield.lat_grid", gfield.get_numeric_grid(ilat), 1);
        else
            chk_if_equal("lat_grid", "gfield.lat_grid", lat_grid, gfield.get_numeric_grid(ilat));
    }
    else
    {
        chk_if_equal("lat_grid", "gfield.lat_grid", lat_grid,
                     gfield.get_numeric_grid(ilat));
        chk_if_equal("lon_grid", "gfield.lon_grid", lon_grid,
                     gfield.get_numeric_grid(ilon));
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(// WS Generic Output:
                           Matrix& field_out,
                           // WS Input:
                           const Vector& p_grid _U_,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           // WS Generic Input:
                           const GriddedField2& gfraw_in,
                           const Verbosity&)
{
    FieldFromGriddedFieldCheckLatLonHelper(lat_grid, lon_grid, 0, 1, gfraw_in);

    field_out = gfraw_in.data;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(// WS Generic Output:
                           Tensor3& field_out,
                           // WS Input:
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           // WS Generic Input:
                           const GriddedField3& gfraw_in,
                           const Verbosity&)
{
    chk_griddedfield_gridname(gfraw_in, 0, "Pressure");
    chk_if_equal("p_grid", "gfield.p_grid", p_grid, gfraw_in.get_numeric_grid(0));

    FieldFromGriddedFieldCheckLatLonHelper(lat_grid, lon_grid, 1, 2, gfraw_in);

    field_out = gfraw_in.data;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(// WS Generic Output:
                           Tensor4& field_out,
                           // WS Input:
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           // WS Generic Input:
                           const GriddedField4& gfraw_in,
                           const Verbosity&)
{
    chk_griddedfield_gridname(gfraw_in, 1, "Pressure");
    chk_if_equal("p_grid", "gfield.p_grid", p_grid, gfraw_in.get_numeric_grid(0));

    FieldFromGriddedFieldCheckLatLonHelper(lat_grid, lon_grid, 2, 3, gfraw_in);

    field_out = gfraw_in.data;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void FieldFromGriddedField(// WS Generic Output:
                           Tensor4& field_out,
                           // WS Input:
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           // WS Generic Input:
                           const ArrayOfGriddedField3& gfraw_in,
                           const Verbosity& verbosity)
{
    if (!gfraw_in.nelem())
    {
        CREATE_OUT1;
        out1 << "   Warning: gfraw_in is empty, proceeding anyway\n";
        field_out.resize(0, 0, 0, 0);
    }
    else
    {
        field_out.resize(gfraw_in.nelem(), p_grid.nelem(),
                         gfraw_in[0].data.nrows(), gfraw_in[0].data.ncols());
    }

    for (Index i = 0; i < gfraw_in.nelem(); i++)
    {
        chk_griddedfield_gridname(gfraw_in[i], 0, "Pressure");
        chk_if_equal("p_grid", "gfield.p_grid", p_grid, gfraw_in[i].get_numeric_grid(0));

        FieldFromGriddedFieldCheckLatLonHelper(lat_grid, lon_grid, 1, 2, gfraw_in[i]);

        field_out(i, joker, joker, joker) = gfraw_in[i].data;
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(// WS Generic Output:
                              GriddedField2& gfraw_out,
                              // WS Generic Input:
                              const GriddedField2& gfraw_in_orig,
                              const Verbosity&)
{
    const GriddedField2* gfraw_in_pnt;
    GriddedField2 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField2 &gfraw_in = *gfraw_in_pnt;

    chk_griddedfield_gridname(gfraw_in, 0, "Latitude");
    chk_griddedfield_gridname(gfraw_in, 1, "Longitude");

    if (gfraw_in.data.ncols() != 1 && gfraw_in.data.nrows() != 1)
        throw runtime_error("Can't expand data because number of Latitudes and Longitudes is greater than 1");

    gfraw_out.set_grid_name(0, "Latitude");
    gfraw_out.set_grid_name(1, "Longitude");

    Vector v(2);
    if (gfraw_in.data.nrows() == 1 && gfraw_in.data.ncols() != 1)
    {
        v[0] = -90; v[1] = 90;
        gfraw_out.set_grid(0, v);
        gfraw_out.resize(2, gfraw_in.data.ncols());

        for (Index j = 0; j < gfraw_in.data.ncols(); j++)
            gfraw_out.data(joker, j) = gfraw_in.data(0, j);
    }
    else if (gfraw_in.data.nrows() != 1 && gfraw_in.data.ncols() == 1)
    {
        v[0] = 0; v[1] = 360;
        gfraw_out.set_grid(1, v);
        gfraw_out.resize(gfraw_in.data.nrows(), 2);

        for (Index j = 0; j < gfraw_in.data.nrows(); j++)
            gfraw_out.data(j, joker) = gfraw_in.data(j, 0);
    }
    else
    {
        v[0] = -90; v[1] = 90;
        gfraw_out.set_grid(0, v);
        v[0] = 0; v[1] = 360;
        gfraw_out.set_grid(1, v);
        gfraw_out.resize(2, 2);

        gfraw_out.data = gfraw_in.data(0, 0);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(// WS Generic Output:
                              GriddedField3& gfraw_out,
                              // WS Generic Input:
                              const GriddedField3& gfraw_in_orig,
                              const Verbosity&)
{
    const GriddedField3* gfraw_in_pnt;
    GriddedField3 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField3 &gfraw_in = *gfraw_in_pnt;

    chk_griddedfield_gridname(gfraw_in, 1, "Latitude");
    chk_griddedfield_gridname(gfraw_in, 2, "Longitude");

    if (gfraw_in.data.ncols() != 1 && gfraw_in.data.nrows() != 1)
      throw runtime_error("Can't expand data because number of Latitudes and Longitudes is greater than 1");

    gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
    gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));
    gfraw_out.set_grid_name(1, "Latitude");
    gfraw_out.set_grid_name(2, "Longitude");

    Vector v(2);
    if (gfraw_in.data.nrows() == 1 && gfraw_in.data.ncols() != 1)
    {
        v[0] = -90; v[1] = 90;
        gfraw_out.set_grid(1, v);
        gfraw_out.resize(gfraw_in.data.npages(), 2, gfraw_in.data.ncols());
        
        for (Index i = 0; i < gfraw_in.data.npages(); i++)
            for (Index j = 0; j < gfraw_in.data.ncols(); j++)
                gfraw_out.data(i, joker, j) = gfraw_in.data(i, 0, j);
    }
    else if (gfraw_in.data.nrows() != 1 && gfraw_in.data.ncols() == 1)
    {
        v[0] = 0; v[1] = 360;
        gfraw_out.set_grid(2, v);
        gfraw_out.resize(gfraw_in.data.npages(), gfraw_in.data.nrows(), 2);

        for (Index i = 0; i < gfraw_in.data.npages(); i++)
            for (Index j = 0; j < gfraw_in.data.nrows(); j++)
                gfraw_out.data(i, j, joker) = gfraw_in.data(i, j, 0);
    }
    else
    {
        v[0] = -90; v[1] = 90;
        gfraw_out.set_grid(1, v);
        v[0] = 0; v[1] = 360;
        gfraw_out.set_grid(2, v);
        gfraw_out.resize(gfraw_in.data.npages(), 2, 2);

        for (Index i = 0; i < gfraw_in.data.npages(); i++)
            gfraw_out.data(i, joker, joker) = gfraw_in.data(i, 0, 0);
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(// WS Generic Output:
                              GriddedField4& gfraw_out,
                              // WS Generic Input:
                              const GriddedField4& gfraw_in_orig,
                              const Verbosity&)
{
    const GriddedField4* gfraw_in_pnt;
    GriddedField4 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField4 &gfraw_in = *gfraw_in_pnt;

    chk_griddedfield_gridname(gfraw_in, 2, "Latitude");
    chk_griddedfield_gridname(gfraw_in, 3, "Longitude");

    if (gfraw_in.data.ncols() != 1 && gfraw_in.data.nrows() != 1)
        throw runtime_error("Can't expand data because number of Latitudes and Longitudes is greater than 1");
    
    gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
    gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));

    gfraw_out.set_grid(1, gfraw_in.get_numeric_grid(1));
    gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));

    gfraw_out.set_grid_name(2, "Latitude");
    gfraw_out.set_grid_name(3, "Longitude");

    Vector v(2);
    if (gfraw_in.data.nrows() == 1 && gfraw_in.data.ncols() != 1)
    {
        v[0] = -90; v[1] = 90;
        gfraw_out.set_grid(2, v);
        gfraw_out.resize(gfraw_in.data.nbooks(), gfraw_in.data.npages(), 2, gfraw_in.data.ncols());

        for (Index k = 0; k < gfraw_in.data.nbooks(); k++)
            for (Index i = 0; i < gfraw_in.data.npages(); i++)
                for (Index j = 0; j < gfraw_in.data.ncols(); j++)
                    gfraw_out.data(k, i, joker, j) = gfraw_in.data(k, i, 0, j);
    }
    else if (gfraw_in.data.nrows() != 1 && gfraw_in.data.ncols() == 1)
    {
        v[0] = 0; v[1] = 360;
        gfraw_out.set_grid(3, v);
        gfraw_out.resize(gfraw_in.data.nbooks(), gfraw_in.data.npages(), gfraw_in.data.nrows(), 2);

        for (Index k = 0; k < gfraw_in.data.nbooks(); k++)
            for (Index i = 0; i < gfraw_in.data.npages(); i++)
                for (Index j = 0; j < gfraw_in.data.nrows(); j++)
                    gfraw_out.data(k, i, j, joker) = gfraw_in.data(k, i, j, 0);
    }
    else
    {
        v[0] = -90; v[1] = 90;
        gfraw_out.set_grid(2, v);
        v[0] = 0; v[1] = 360;
        gfraw_out.set_grid(3, v);
        gfraw_out.resize(gfraw_in.data.nbooks(), gfraw_in.data.npages(), 2, 2);

        for (Index k = 0; k < gfraw_in.data.nbooks(); k++)
            for (Index i = 0; i < gfraw_in.data.npages(); i++)
                gfraw_out.data(k, i, joker, joker) = gfraw_in.data(k, i, 0, 0);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonExpand(// WS Generic Output:
                              ArrayOfGriddedField3& gfraw_out,
                              // WS Generic Input:
                              const ArrayOfGriddedField3& gfraw_in,
                              const Verbosity& verbosity)
{
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
                               ArrayOfGridPosPoly& gp_p,
                               Matrix& itw,
                               GriddedField& gfraw_out,
                               const GriddedField& gfraw_in,
                               const Index p_grid_index,
                               ConstVectorView p_grid,
                               const Index& interp_order,
                               const Index& zeropadding,
                               const Verbosity& verbosity)
{
    CREATE_OUT2;

    chk_griddedfield_gridname(gfraw_in, p_grid_index, "Pressure");

    out2 << "  Interpolation order: " << interp_order << "\n";

    const ConstVectorView in_p_grid = gfraw_in.get_numeric_grid(p_grid_index);

    // Initialize output field. Set grids and copy grid names
    gfraw_out.set_grid(p_grid_index, p_grid);
    gfraw_out.set_grid_name(p_grid_index, gfraw_in.get_grid_name(p_grid_index));

    if (zeropadding)
      {
        if (in_p_grid[0]<p_grid[p_grid.nelem()-1] ||
            in_p_grid[in_p_grid.nelem()-1]>p_grid[0])
          {
            ing_min = 0;
            ing_max = ing_min-1;
          }
        else
          chk_interpolation_pgrids_loose_no_data_check(ing_min, ing_max,
                                                       "Raw field to p_grid",
                                                       in_p_grid,
                                                       p_grid,
                                                       interp_order);
      }
    else
    {
        ing_min = 0;
        ing_max = p_grid.nelem()-1;
        chk_interpolation_pgrids("Raw field to p_grid",
                                 in_p_grid,
                                 p_grid,
                                 interp_order);
    }
    Index nelem_in_range = ing_max - ing_min + 1;

    // Calculate grid positions:
    if (nelem_in_range>0)
      {
        gp_p.resize(nelem_in_range);
        p2gridpos_poly(gp_p, in_p_grid, p_grid[Range(ing_min, nelem_in_range)],
                       interp_order);

        // Interpolation weights:
        itw.resize(nelem_in_range, interp_order+1);
        interpweights(itw, gp_p);
      }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldPRegrid(// WS Generic Output:
                              GriddedField3& gfraw_out,
                              // WS Input:
                              const Vector& p_grid,
                              // WS Generic Input:
                              const GriddedField3& gfraw_in_orig,
                              const Index& interp_order,
                              const Index& zeropadding,
                              const Verbosity& verbosity)
{
    const GriddedField3* gfraw_in_pnt;
    GriddedField3 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField3 &gfraw_in = *gfraw_in_pnt;

    const Index p_grid_index = 0;

    // Resize output GriddedField and copy all non-latitude/longitude grids
    gfraw_out.resize(p_grid.nelem(), gfraw_in.data.nrows(), gfraw_in.data.ncols());
    gfraw_out.set_grid(1, gfraw_in.get_numeric_grid(1));
    gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));
    gfraw_out.set_grid(2, gfraw_in.get_numeric_grid(2));
    gfraw_out.set_grid_name(2, gfraw_in.get_grid_name(2));

    ArrayOfGridPosPoly gp_p;
    Matrix itw;

    Index ing_min, ing_max;

    GriddedFieldPRegridHelper(ing_min, ing_max,
                              gp_p, itw, gfraw_out,
                              gfraw_in, p_grid_index,
                              p_grid, interp_order, zeropadding,
                              verbosity);

    // Interpolate:
    if (ing_max - ing_min < 0)
        gfraw_out.data = 0.;
    else
      if (ing_max - ing_min + 1 != p_grid.nelem())
      {
        gfraw_out.data = 0.;
        for (Index i = 0; i < gfraw_in.data.nrows(); i++)
            for (Index j = 0; j < gfraw_in.data.ncols(); j++)
            {
//                chk_interpolation_grids_loose_check_data(ing_min, ing_max,
//                                                         "Raw field to p_grid",
//                                                         gfraw_in.get_numeric_grid(p_grid_index),
//                                                         p_grid, gfraw_in.data(joker, i, j));
                interp(gfraw_out.data(Range(ing_min, ing_max-ing_min+1), i, j),
                       itw, gfraw_in.data(joker, i, j), gp_p);
            }
      }
      else
        for (Index i = 0; i < gfraw_in.data.nrows(); i++)
            for (Index j = 0; j < gfraw_in.data.ncols(); j++)
                interp(gfraw_out.data(joker, i, j),
                       itw, gfraw_in.data(joker, i, j), gp_p);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldPRegrid(// WS Generic Output:
                              GriddedField4& gfraw_out,
                              // WS Input:
                              const Vector& p_grid,
                              // WS Generic Input:
                              const GriddedField4& gfraw_in_orig,
                              const Index& interp_order,
                              const Index& zeropadding,
                              const Verbosity& verbosity)
{
    const GriddedField4* gfraw_in_pnt;
    GriddedField4 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField4 &gfraw_in = *gfraw_in_pnt;

    const Index p_grid_index = 1;

    // Resize output GriddedField and copy all non-latitude/longitude grids
    gfraw_out.resize(gfraw_in.data.nbooks(), p_grid.nelem(),
                     gfraw_in.data.nrows(), gfraw_in.data.ncols());
    gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
    gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));
    gfraw_out.set_grid(2, gfraw_in.get_numeric_grid(2));
    gfraw_out.set_grid_name(2, gfraw_in.get_grid_name(2));
    gfraw_out.set_grid(3, gfraw_in.get_numeric_grid(3));
    gfraw_out.set_grid_name(3, gfraw_in.get_grid_name(3));

    ArrayOfGridPosPoly gp_p;
    Matrix itw;

    Index ing_min, ing_max;

    GriddedFieldPRegridHelper(ing_min, ing_max,
                              gp_p, itw, gfraw_out,
                              gfraw_in, p_grid_index,
                              p_grid, interp_order, zeropadding,
                              verbosity);

    // Interpolate:
    if (ing_max - ing_min < 0)
        gfraw_out.data = 0.;
    else
      if (ing_max - ing_min + 1 != p_grid.nelem())
      {
        gfraw_out.data = 0.;
        for (Index b = 0; b < gfraw_in.data.nbooks(); b++)
            for (Index i = 0; i < gfraw_in.data.nrows(); i++)
                for (Index j = 0; j < gfraw_in.data.ncols(); j++)
                {
//                    chk_interpolation_grids_loose_check_data(ing_min, ing_max,
//                                                             "Raw field to p_grid",
//                                                             gfraw_in.get_numeric_grid(p_grid_index),
//                                                             p_grid, gfraw_in.data(b, joker, i, j));
                    interp(gfraw_out.data(b, Range(ing_min, ing_max-ing_min), i, j),
                           itw, gfraw_in.data(b, joker, i, j), gp_p);
                }
      }
      else
        for (Index b = 0; b < gfraw_in.data.nbooks(); b++)
            for (Index i = 0; i < gfraw_in.data.nrows(); i++)
                for (Index j = 0; j < gfraw_in.data.ncols(); j++)
                    interp(gfraw_out.data(b, joker, i, j),
                           itw, gfraw_in.data(b, joker, i, j), gp_p);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldPRegrid(// WS Generic Output:
                              ArrayOfGriddedField3& agfraw_out,
                              // WS Input:
                              const Vector& p_grid,
                              // WS Generic Input:
                              const ArrayOfGriddedField3& agfraw_in,
                              const Index& interp_order,
                              const Index& zeropadding,
                              const Verbosity& verbosity)
{
    agfraw_out.resize(agfraw_in.nelem());

    for (Index i = 0; i < agfraw_in.nelem(); i++)
    {
        GriddedFieldPRegrid(agfraw_out[i], p_grid, agfraw_in[i],
                            interp_order, zeropadding,
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
void GriddedFieldLatLonRegridHelper(ArrayOfGridPosPoly& gp_lat,
                                    ArrayOfGridPosPoly& gp_lon,
                                    Tensor3& itw,
                                    GriddedField& gfraw_out,
                                    const GriddedField& gfraw_in,
                                    const Index lat_grid_index,
                                    const Index lon_grid_index,
                                    ConstVectorView lat_true,
                                    ConstVectorView lon_true,
                                    const Index& interp_order,
                                    const Verbosity& verbosity)
{
    CREATE_OUT2;

    if (!lat_true.nelem()) throw runtime_error("The new latitude grid is not allowed to be empty.");
    if (!lon_true.nelem()) throw runtime_error("The new longitude grid is not allowed to be empty.");

    chk_griddedfield_gridname(gfraw_in, lat_grid_index, "Latitude");
    chk_griddedfield_gridname(gfraw_in, lon_grid_index, "Longitude");
    if (gfraw_in.get_grid_size(lat_grid_index) == 1 || gfraw_in.get_grid_size(lon_grid_index) == 1)
        throw runtime_error("Raw data has to be true 3D data (nlat>1 and nlon>1).");

    out2 << "  Interpolation order: " << interp_order << "\n";

    const ConstVectorView in_lat_grid = gfraw_in.get_numeric_grid(lat_grid_index);
    const ConstVectorView in_lon_grid = gfraw_in.get_numeric_grid(lon_grid_index);

    // Initialize output field. Set grids and copy grid names
    gfraw_out.set_grid(lat_grid_index, lat_true);
    gfraw_out.set_grid_name(lat_grid_index, gfraw_in.get_grid_name(lat_grid_index));
    gfraw_out.set_grid(lon_grid_index, lon_true);
    gfraw_out.set_grid_name(lon_grid_index, gfraw_in.get_grid_name(lon_grid_index));

    chk_interpolation_grids("Raw field to lat_grid, 3D case",
                            in_lat_grid,
                            lat_true,
                            interp_order);

    // Calculate grid positions:
    gp_lat.resize(lat_true.nelem());
    gp_lon.resize(lon_true.nelem());
    gridpos_poly(gp_lat, in_lat_grid, lat_true, interp_order);

    gridpos_poly_longitudinal("Raw field to lon_grid, 3D case",
                              gp_lon, in_lon_grid, lon_true, interp_order);

    // Interpolation weights:
    itw.resize(lat_true.nelem(), lon_true.nelem(), (interp_order+1)*(interp_order+1));
    interpweights(itw, gp_lat, gp_lon);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(// WS Generic Output:
                              GriddedField2& gfraw_out,
                              // WS Input:
                              const Vector& lat_true,
                              const Vector& lon_true,
                              // WS Generic Input:
                              const GriddedField2& gfraw_in_orig,
                              const Index& interp_order,
                              const Verbosity& verbosity)
{
    if (!lat_true.nelem()) throw runtime_error("The new latitude grid is not allowed to be empty.");
    if (!lon_true.nelem()) throw runtime_error("The new longitude grid is not allowed to be empty.");

    const GriddedField2* gfraw_in_pnt;
    GriddedField2 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField2 &gfraw_in = *gfraw_in_pnt;

    const Index lat_grid_index = 0;
    const Index lon_grid_index = 1;

    if (gfraw_in.get_grid_size(lat_grid_index) < 2 || 
        gfraw_in.get_grid_size(lon_grid_index) < 2)
      {
        ostringstream os;
        os << "Raw data has to be true 3D data (nlat>1 and nlon>1).\n"
           << "Use GriddedFieldLatLonExpand to convert 1D or 2D data to 3D!\n"; 
        throw runtime_error(os.str());
      }

    // Resize output GriddedField
    gfraw_out.resize(lat_true.nelem(), lon_true.nelem());

    ArrayOfGridPosPoly gp_lat;
    ArrayOfGridPosPoly gp_lon;
    Tensor3 itw;

    // If lon grid is cyclic, the data values at 0 and 360 must match
    const ConstVectorView& in_lat_grid = gfraw_in.get_numeric_grid(lat_grid_index);
    const ConstVectorView& in_lon_grid = gfraw_in.get_numeric_grid(lon_grid_index);

    if (is_lon_cyclic(in_lon_grid))
    {
        for (Index lat = 0; lat < in_lat_grid.nelem(); lat++)
        {
            if (!is_same_within_epsilon(gfraw_in.data(lat, 0),
                                        gfraw_in.data(lat, in_lon_grid.nelem()-1),
                                        EPSILON_LON_CYCLIC))
            {
                ostringstream os;
                os << "Data values at 0 and 360 degrees for a cyclic longitude grid must match: \n"
                << "Mismatch at latitude index    : "
                << lat << " (" << in_lat_grid[lat] << " degrees)\n"
                << "Value at 0 degrees longitude  : "
                << gfraw_in.data(lat, 0) << "\n"
                << "Value at 360 degrees longitude: "
                << gfraw_in.data(lat, in_lon_grid.nelem()-1) << "\n"
                << "Difference                    : "
                << gfraw_in.data(lat, in_lon_grid.nelem()-1) - gfraw_in.data(lat, 0) << "\n"
                << "Allowed difference            : " << EPSILON_LON_CYCLIC;
                throw runtime_error(os.str());
            }
        }
    }

    GriddedFieldLatLonRegridHelper(gp_lat, gp_lon, itw, gfraw_out,
                                   gfraw_in, lat_grid_index, lon_grid_index,
                                   lat_true, lon_true, interp_order,
                                   verbosity);

    // Interpolate:
    interp(gfraw_out.data, itw, gfraw_in.data, gp_lat, gp_lon);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(// WS Generic Output:
                              GriddedField3& gfraw_out,
                              // WS Input:
                              const Vector& lat_true,
                              const Vector& lon_true,
                              // WS Generic Input:
                              const GriddedField3& gfraw_in_orig,
                              const Index& interp_order,
                              const Verbosity& verbosity)
{
    if (!lat_true.nelem()) throw runtime_error("The new latitude grid is not allowed to be empty.");
    if (!lon_true.nelem()) throw runtime_error("The new longitude grid is not allowed to be empty.");

    const GriddedField3* gfraw_in_pnt;
    GriddedField3 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField3 &gfraw_in = *gfraw_in_pnt;

    const Index lat_grid_index = 1;
    const Index lon_grid_index = 2;

    if (gfraw_in.get_grid_size(lat_grid_index) < 2 || 
        gfraw_in.get_grid_size(lon_grid_index) < 2)
      {
        ostringstream os;
        os << "Raw data has to be true 3D data (nlat>1 and nlon>1).\n"
           << "Use GriddedFieldLatLonExpand to convert 1D or 2D data to 3D!\n"; 
        throw runtime_error(os.str());
      }

    // Resize output GriddedField and copy all non-latitude/longitude grids
    gfraw_out.resize(gfraw_in.data.npages(), lat_true.nelem(), lon_true.nelem());
    gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
    gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));

    ArrayOfGridPosPoly gp_lat;
    ArrayOfGridPosPoly gp_lon;
    Tensor3 itw;

    // If lon grid is cyclic, the data values at 0 and 360 must match
    const ConstVectorView& in_grid0 = gfraw_in.get_numeric_grid(0);
    const ConstVectorView& in_lat_grid = gfraw_in.get_numeric_grid(lat_grid_index);
    const ConstVectorView& in_lon_grid = gfraw_in.get_numeric_grid(lon_grid_index);

    if (is_lon_cyclic(in_lon_grid))
    {
        for (Index g0 = 0; g0 < in_grid0.nelem(); g0++)
            for (Index lat = 0; lat < in_lat_grid.nelem(); lat++)
            {
                if (!is_same_within_epsilon(gfraw_in.data(g0, lat, 0),
                                            gfraw_in.data(g0, lat, in_lon_grid.nelem()-1),
                                            EPSILON_LON_CYCLIC))
                {
                    ostringstream os;
                    os << "Data values at 0 and 360 degrees for a cyclic longitude grid must match: \n"
                    << "Mismatch at 1st grid index    : "
                    << g0 << " (" << in_grid0[g0] << ")\n"
                    << "         at latitude index    : "
                    << lat << " (" << in_lat_grid[lat] << " degrees)\n"
                    << "Value at 0 degrees longitude  : "
                    << gfraw_in.data(g0, lat, 0) << "\n"
                    << "Value at 360 degrees longitude: "
                    << gfraw_in.data(g0, lat, in_lon_grid.nelem()-1) << "\n"
                    << "Difference                    : "
                    << gfraw_in.data(g0, lat, in_lon_grid.nelem()-1) - gfraw_in.data(g0, lat, 0) << "\n"
                    << "Allowed difference            : " << EPSILON_LON_CYCLIC;
                    throw runtime_error(os.str());
                }
            }
    }

    GriddedFieldLatLonRegridHelper(gp_lat, gp_lon, itw, gfraw_out,
                                   gfraw_in, lat_grid_index, lon_grid_index,
                                   lat_true, lon_true, interp_order,
                                   verbosity);

    // Interpolate:
    for (Index i = 0; i < gfraw_in.data.npages(); i++)
        interp(gfraw_out.data(i, joker, joker),
               itw, gfraw_in.data(i, joker, joker), gp_lat, gp_lon);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(// WS Generic Output:
                              GriddedField4& gfraw_out,
                              // WS Input:
                              const Vector& lat_true,
                              const Vector& lon_true,
                              // WS Generic Input:
                              const GriddedField4& gfraw_in_orig,
                              const Index& interp_order,
                              const Verbosity& verbosity)
{
    if (!lat_true.nelem()) throw runtime_error("The new latitude grid is not allowed to be empty.");
    if (!lon_true.nelem()) throw runtime_error("The new longitude grid is not allowed to be empty.");

    const GriddedField4* gfraw_in_pnt;
    GriddedField4 gfraw_in_copy;

    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;

    const GriddedField4 &gfraw_in = *gfraw_in_pnt;

    const Index lat_grid_index = 2;
    const Index lon_grid_index = 3;

    if (gfraw_in.get_grid_size(lat_grid_index) < 2 || 
        gfraw_in.get_grid_size(lon_grid_index) < 2)
      {
        ostringstream os;
        os << "Raw data has to be true 3D data (nlat>1 and nlon>1).\n"
           << "Use GriddedFieldLatLonExpand to convert 1D or 2D data to 3D!\n"; 
        throw runtime_error(os.str());
      }

    // Resize output GriddedField and copy all non-latitude/longitude grids
    gfraw_out.resize(gfraw_in.data.nbooks(), gfraw_in.data.npages(),
                     lat_true.nelem(), lon_true.nelem());
    gfraw_out.set_grid(0, gfraw_in.get_numeric_grid(0));
    gfraw_out.set_grid_name(0, gfraw_in.get_grid_name(0));
    gfraw_out.set_grid(1, gfraw_in.get_numeric_grid(1));
    gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));

    ArrayOfGridPosPoly gp_lat;
    ArrayOfGridPosPoly gp_lon;
    Tensor3 itw;

    GriddedFieldLatLonRegridHelper(gp_lat, gp_lon, itw, gfraw_out,
                                   gfraw_in, lat_grid_index, lon_grid_index,
                                   lat_true, lon_true, interp_order,
                                   verbosity);

    // If lon grid is cyclic, the data values at 0 and 360 must match
    const ConstVectorView& in_grid0 = gfraw_in.get_numeric_grid(0);
    const ConstVectorView& in_grid1 = gfraw_in.get_numeric_grid(1);
    const ConstVectorView& in_lat_grid = gfraw_in.get_numeric_grid(lat_grid_index);
    const ConstVectorView& in_lon_grid = gfraw_in.get_numeric_grid(lon_grid_index);

    if (is_lon_cyclic(in_lon_grid))
    {
        for (Index g0 = 0; g0 < in_grid0.nelem(); g0++)
            for (Index g1 = 0; g1 < in_grid1.nelem(); g1++)
                for (Index lat = 0; lat < in_lat_grid.nelem(); lat++)
                {
                    if (!is_same_within_epsilon(gfraw_in.data(g0, g1, lat, 0),
                                                gfraw_in.data(g0, g1, lat, in_lon_grid.nelem()-1),
                                                EPSILON_LON_CYCLIC))
                    {
                        ostringstream os;
                        os << "Data values at 0 and 360 degrees for a cyclic longitude grid must match: \n"
                        << "Mismatch at 1st grid index    : "
                        << g0 << " (" << in_grid0[g0] << ")\n"
                        << "         at 2nd grid index    : "
                        << g1 << " (" << in_grid1[g1] << ")\n"
                        << "         at latitude index    : "
                        << lat << " (" << in_lat_grid[lat] << " degrees)\n"
                        << "Value at 0 degrees longitude  : "
                        << gfraw_in.data(g0, g1, lat, 0) << "\n"
                        << "Value at 360 degrees longitude: "
                        << gfraw_in.data(g0, g1, lat, in_lon_grid.nelem()-1) << "\n"
                        << "Difference                    : "
                        << gfraw_in.data(g0, g1, lat, in_lon_grid.nelem()-1) - gfraw_in.data(g0, g1, lat, 0) << "\n"
                        << "Allowed difference            : " << EPSILON_LON_CYCLIC;
                        throw runtime_error(os.str());
                    }
            }
    }

    // Interpolate:
    for (Index i = 0; i < gfraw_in.data.nbooks(); i++)
        for (Index j = 0; j < gfraw_in.data.npages(); j++)
            interp(gfraw_out.data(i, j, joker, joker),
                   itw, gfraw_in.data(i, j, joker, joker), gp_lat, gp_lon);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldLatLonRegrid(// WS Generic Output:
                              ArrayOfGriddedField3& agfraw_out,
                              // WS Input:
                              const Vector& lat_true,
                              const Vector& lon_true,
                              // WS Generic Input:
                              const ArrayOfGriddedField3& agfraw_in,
                              const Index& interp_order,
                              const Verbosity& verbosity)
{
    agfraw_out.resize(agfraw_in.nelem());

    for (Index i = 0; i < agfraw_in.nelem(); i++)
    {
        GriddedFieldLatLonRegrid(agfraw_out[i], lat_true, lon_true, agfraw_in[i],
                                 interp_order, verbosity);
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
                                  ArrayOfGridPosPoly& gp_p,
                                  Matrix& itw,
                                  const GriddedField& gfraw_in,
                                  const Index z_grid_index,
                                  ConstVectorView z_grid,
                                  const Index& interp_order,
                                  const Index& zeropadding,
                                  const Verbosity& verbosity)
{
    CREATE_OUT2;

    chk_griddedfield_gridname(gfraw_in, z_grid_index, "Altitude");

    out2 << "  Interpolation order: " << interp_order << "\n";

    const ConstVectorView in_z_grid = gfraw_in.get_numeric_grid(z_grid_index);

    if (zeropadding)
    {
        if (in_z_grid[0] > z_grid[z_grid.nelem()-1] ||
            in_z_grid[in_z_grid.nelem()-1] < z_grid[0])
        {
            ing_min = 0;
            ing_max = ing_min-1;
        }
        else
            chk_interpolation_grids_loose_no_data_check(ing_min, ing_max,
                                                         "Raw field to z_grid",
                                                         in_z_grid,
                                                         z_grid,
                                                         interp_order);
    }
    else
    {
        ing_min = 0;
        ing_max = z_grid.nelem()-1;
        chk_interpolation_pgrids("Raw field to p_grid",
                                 in_z_grid,
                                 z_grid,
                                 interp_order);
    }

    Index nelem_in_range = ing_max - ing_min + 1;

    // Calculate grid positions:
    if (nelem_in_range>0)
    {
        gp_p.resize(nelem_in_range);
        gridpos_poly(gp_p, in_z_grid, z_grid[Range(ing_min, nelem_in_range)],
                     interp_order);

        // Interpolation weights:
        itw.resize(nelem_in_range, interp_order+1);
        interpweights(itw, gp_p);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void GriddedFieldZToPRegrid(// WS Generic Output:
                            GriddedField3& gfraw_out, //grid in P
                            // WS Input:
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid,
                            const Tensor3& z_field,
                            // WS Generic Input:
                            const GriddedField3& gfraw_in_orig,//grid in Z
                            const Index& interp_order, // Only linear interpolation allowed
                            const Index& zeropadding,
                            const Verbosity& verbosity)
{
    
    // z_field must be of the same size as its grids
    if(!((z_field.npages()==p_grid.nelem() && z_field.nrows() == lat_grid.nelem()) && z_field.ncols() == lon_grid.nelem()))
        throw std::runtime_error("*z_field* must be of the same size as *p_grid*, *lat_grid*, and *lon_grid* in *GriddedFieldZToPRegrid*.");
    
    // Must name the dimension "Altitude" to ensure user is aware of what they are doing.
    chk_griddedfield_gridname(gfraw_in_orig, 0, "Altitude");
    
    // Each lat and lon grid must be identical between z_field and in field
    const Vector& lat_in = gfraw_in_orig.get_numeric_grid(1);
    const Vector& lon_in = gfraw_in_orig.get_numeric_grid(2);
    
    if(lat_grid.nelem()!=lat_in.nelem() || lon_grid.nelem()!=lon_in.nelem())
        throw std::runtime_error("Gridding of field to regrid is bad.\n*GriddedFieldZToPRegrid* requires latitude and longitude to be on the same grid as *z_field*.");
    else
    {
        for(Index ii=0;ii<lat_grid.nelem();ii++)
            if(lat_grid[ii]!=lat_in[ii])
                throw std::runtime_error("Gridding of field to regrid is bad.\n*GriddedFieldZToPRegrid* requires latitude and longitude to be the same as for *z_field*.");
        for(Index ii=0;ii<lon_grid.nelem();ii++)
            if(lon_grid[ii]!=lon_in[ii])
                throw std::runtime_error("Gridding of field to regrid is bad.\n*GriddedFieldZToPRegrid* requires latitude and longitude to be the same as for *z_field*.");
    }
    
    // Pointer in case output is input variable (same memory allocated)
    const GriddedField3* gfraw_in_pnt;
    GriddedField3 gfraw_in_copy;
    
    if (&gfraw_in_orig == &gfraw_out)
    {
        gfraw_in_copy = gfraw_in_orig;
        gfraw_in_pnt = &gfraw_in_copy;
    }
    else
        gfraw_in_pnt = &gfraw_in_orig;
    
    // Now output and input are separate variables (not allocating the same memory)
    const GriddedField3 &gfraw_in = *gfraw_in_pnt;
    
    // Right size and order
    gfraw_out.resize(p_grid.nelem(),lat_grid.nelem(),lon_grid.nelem());
    gfraw_out.set_grid(0, p_grid);
    gfraw_out.set_grid_name(0, "Pressure");
    gfraw_out.set_grid(1, lat_grid);
    gfraw_out.set_grid_name(1, gfraw_in.get_grid_name(1));
    gfraw_out.set_grid(2, lon_grid);
    gfraw_out.set_grid_name(2, gfraw_in.get_grid_name(2));
    gfraw_out.data = 0.;

    ArrayOfGridPosPoly gp_p;
    Matrix itw;

    Index ing_min, ing_max;

    for(Index lat_index = 0;lat_index<lat_grid.nelem();lat_index++)
    {
        for(Index lon_index = 0;lon_index<lon_grid.nelem();lon_index++)
        {
            const Vector z_out = z_field(joker,lat_index,lon_index);

            GriddedFieldZToPRegridHelper(ing_min, ing_max,
                                         gp_p, itw,
                                         gfraw_in, 0,
                                         z_out, interp_order, zeropadding,
                                         verbosity);

            if (ing_max - ing_min >= 0)
            {
                Range r = joker;

                if (ing_max - ing_min + 1 != z_out.nelem())
                {
                    r = Range(ing_min, ing_max-ing_min+1);
                }

                interp(gfraw_out.data(r, lat_index, lon_index),
                       itw, gfraw_in.data(joker, lat_index, lon_index), gp_p);
            }
        }
    }
}


// Workspace method, doxygen header will be auto-generated.
// 2007-07-25 Stefan Buehler
void atm_fields_compactFromMatrix(// WS Output:
                                  GriddedField4& af, // atm_fields_compact
                                  // WS Input:
                                  const Index& atmosphere_dim,
                                  // WS Generic Input:
                                  const Matrix& im,
                                  // Control Parameters:
                                  const ArrayOfString& field_names,
                                  const Verbosity&)
{
  if (1!=atmosphere_dim)
    {
      ostringstream os; 
      os << "Atmospheric dimension must be 1.";
      throw runtime_error( os.str() );
    }

  const Index np = im.nrows();   // Number of pressure levels.
  const Index nf = im.ncols()-1; // Total number of fields.
  ArrayOfIndex f_1; // indices of non-ignored fields 
                   // All fields called "ignore" will be ignored.
  String fn_upper; // Temporary variable to hold upper case field_names.

  if (field_names.nelem()!=nf)
    {
      ostringstream os; 
      os << "Cannot extract fields from Matrix.\n"
         << "*field_names* must have one element less than there are\n"
         << "matrix columns.";
      throw runtime_error( os.str() );
    }


  // Remove additional fields from the field_names. All fields that
  // are flagged by 'ignore' in the field names, small or large letters,
  // are removed.
  for(Index f = 0; f < field_names.nelem(); f++)
    {
      fn_upper = field_names[f];
      fn_upper.toupper();
      //cout << "fieldname[" << f << "]: " << fn_upper;
      if(fn_upper != "IGNORE")
      {
        f_1.push_back(f);
        //cout << " put as element " << f_1.size()-1 << " into selection\n";
      }
    }

  // Copy required field_names to a new variable called field_names_1
  Index nf_1=f_1.size(); // Number of required fields. 
  ArrayOfString field_names_1(nf_1);
  for (Index f=0; f<nf_1; f++)
  {
    field_names_1[f] = field_names[f_1[f]];
    //cout << "new fieldname[" << f << "] (formerly [" << f_1[f] << "]): "
    //     << field_names_1[f] << "\n";
  }


  //  out3 << "Copying *" << im_name << "* to *atm_fields_compact*.\n";
  
  af.set_grid(GFIELD4_FIELD_NAMES, field_names_1);

  af.set_grid(GFIELD4_P_GRID, im(Range(joker),0));
  
  af.set_grid(GFIELD4_LAT_GRID, Vector());
  af.set_grid(GFIELD4_LON_GRID, Vector());
  
  af.resize(nf_1,np,1,1); // Resize it according to the required fields
  for (Index f = 0; f < nf_1; f++)
      af.data(f,Range(joker),0,0) = im(Range(joker),f_1[f]+1);
}


// Workspace method, doxygen header is auto-generated.
// 2007-07-31 Stefan Buehler
// 2011-05-04 Adapted by Gerrit Holl
void atm_fields_compactAddConstant(// WS Output:
                                   GriddedField4& af,
                                   // Control Parameters:
                                   const String& name,
                                   const Numeric& value,
                                   const Index& prepend,
                                   const ArrayOfString& condensibles,
                                   const Verbosity& verbosity)
{
    Index nf; // Will hold new size

    // Add book
    atm_fields_compactExpand(af, nf, name, prepend, verbosity);

    if (condensibles.nelem())
    {
        const Tensor4& vmrs = af.data;
        const ArrayOfString& species = af.get_string_grid(GFIELD4_FIELD_NAMES);
        Tensor3 condensible_sum(vmrs.npages(), vmrs.nrows(), vmrs.ncols(), 1.);
        for (Index c = 0; c < condensibles.nelem(); c++)
        {
            bool species_found = false;
            for (Index i = 0; !species_found && i < species.nelem(); i++)
            {
                if (species[i] == condensibles[c])
                {
                    condensible_sum -= vmrs(i, joker, joker, joker);
                    species_found = true;
                }
            }
            if (!species_found)
            {
                std::ostringstream os;
                os << "Condensible species \"" << condensibles[c] << "\" not found ";
                os << "in input data.";
                throw runtime_error(os.str());
            }
        }
        condensible_sum *= value;
        // Add the constant value:
        if (prepend)
            af.data(0, joker, joker, joker) = condensible_sum;
        else
            af.data(nf-1, joker, joker, joker) = condensible_sum;
    }
    else
    {
        // Add the constant value:
        if (prepend)
            af.data(0, joker, joker, joker) = value;
        else
            af.data(nf-1, joker, joker, joker) = value;
    }
}

// Workspace method, doxygen header is auto-generated
// 2011-05-02 Gerrit Holl
void atm_fields_compactAddSpecies(// WS Output:
                                  GriddedField4& atm_fields_compact,
                                  // WS Generic Input:
                                  const String& name,
                                  const GriddedField3& species,
                                  const Index& prepend,
                                  const Verbosity& verbosity)
{

  assert(atm_fields_compact.checksize());
  assert(species.checksize());

  ConstVectorView af_p_grid = atm_fields_compact.get_numeric_grid(GFIELD4_P_GRID);
  ConstVectorView af_lat_grid = atm_fields_compact.get_numeric_grid(GFIELD4_LAT_GRID);
  ConstVectorView af_lon_grid = atm_fields_compact.get_numeric_grid(GFIELD4_LON_GRID);
  ConstVectorView sp_p_grid = species.get_numeric_grid(GFIELD3_P_GRID);
  ConstVectorView sp_lat_grid = species.get_numeric_grid(GFIELD3_LAT_GRID);
  ConstVectorView sp_lon_grid = species.get_numeric_grid(GFIELD3_LON_GRID);

  Index new_n_fields; // To be set in next line
  atm_fields_compactExpand(atm_fields_compact, new_n_fields, name, prepend, verbosity);

  const Index insert_pos = (prepend)?0:new_n_fields-1;

  // Interpolate species to atm_fields_compact

  // Common for all dim
  chk_interpolation_grids("species p_grid to atm_fields_compact p_grid",
                          sp_p_grid,
                          af_p_grid);
  ArrayOfGridPos p_gridpos(af_p_grid.nelem());
  // gridpos(p_gridpos, sp_p_grid, af_p_grid);
  p2gridpos(p_gridpos, sp_p_grid, af_p_grid);

  if (sp_lat_grid.nelem() > 1)
  {
      // Common for all dim>=2
      chk_interpolation_grids("species lat_grid to atm_fields_compact lat_grid",
                              sp_lat_grid,
                              af_lat_grid);
      ArrayOfGridPos lat_gridpos(af_lat_grid.nelem());
      gridpos(lat_gridpos, sp_lat_grid, af_lat_grid);

      if (sp_lon_grid.nelem() > 1)
      { // 3D-case
          chk_interpolation_grids("species lon_grid to atm_fields_compact lon_grid",
                                  sp_lon_grid,
                                  af_lon_grid);
          ArrayOfGridPos lon_gridpos(af_lon_grid.nelem());
          gridpos(lon_gridpos, sp_lon_grid, af_lon_grid);

          Tensor4 itw(p_gridpos.nelem(), lat_gridpos.nelem(), lon_gridpos.nelem(), 8);
          interpweights(itw, p_gridpos, lat_gridpos, lon_gridpos);

          Tensor3 newfield(af_p_grid.nelem(), af_lat_grid.nelem(), af_lon_grid.nelem());
          interp(newfield, itw, species.data, p_gridpos, lat_gridpos, lon_gridpos);

          atm_fields_compact.data(insert_pos, joker, joker, joker) = newfield;
      } else { // 2D-case
          
          Tensor3 itw(p_gridpos.nelem(), lat_gridpos.nelem(), 4);
          interpweights(itw, p_gridpos, lat_gridpos);

          Matrix newfield(af_p_grid.nelem(), af_lat_grid.nelem());
          interp(newfield, itw, species.data(joker, joker, 0), p_gridpos, lat_gridpos);

          atm_fields_compact.data(insert_pos, joker, joker, 0) = newfield;
      }
  } else { // 1D-case
      Matrix itw(p_gridpos.nelem(), 2);
      interpweights(itw, p_gridpos);

      Vector newfield(af_p_grid.nelem());
      interp(newfield, itw, species.data(joker, 0, 0), p_gridpos);

      atm_fields_compact.data(insert_pos, joker, 0, 0) = newfield;
  }

}


// Workspace method, doxygen header is auto-generated
// 2015-06-30 Jana Mendrok
void atm_fields_compactCleanup (//WS Output:
                          GriddedField4& atm_fields_compact,
                          //WS Input:
                          const Numeric& threshold,
                          const Verbosity&)
{
  assert(atm_fields_compact.checksize());
  Tensor4View afd=atm_fields_compact.data;

  // Check that the data tensor does not contain realistically low (e.g.
  // negative) values. Values smaller than threshold will be set to 0).
  // Ignore T and z, though.
  for ( Index i=0; i<afd.nbooks(); i++ )
    if( atm_fields_compact.get_string_grid(GFIELD4_FIELD_NAMES)[i] != "T" &&
        atm_fields_compact.get_string_grid(GFIELD4_FIELD_NAMES)[i] != "z" )
      for ( Index j=0; j<afd.npages(); j++ )
        for ( Index k=0; k<afd.nrows(); k++ )
          for ( Index l=0; l<afd.ncols(); l++ )
            if ( afd ( i,j,k,l ) < threshold ) 
              afd ( i,j,k,l ) = 0.0;
}


// Workspace method, doxygen header is auto-generated
void atm_fields_compactCreateFromField( // WS Output:
                                        GriddedField4& atm_fields_compact,
                                        // WS Generic Input:
                                        const String& name,
                                        const GriddedField3& field,
                                        const Verbosity&)
{
    assert(field.checksize());
    
    ConstVectorView sp_p_grid = field.get_numeric_grid(GFIELD3_P_GRID);
    ConstVectorView sp_lat_grid = field.get_numeric_grid(GFIELD3_LAT_GRID);
    ConstVectorView sp_lon_grid = field.get_numeric_grid(GFIELD3_LON_GRID);
    ArrayOfString   sp_name_grid(1);
    sp_name_grid[0] = name;
    
    atm_fields_compact.set_grid(0,sp_name_grid);
    atm_fields_compact.set_grid(1,sp_p_grid);
    atm_fields_compact.set_grid(2,sp_lat_grid);
    atm_fields_compact.set_grid(3,sp_lon_grid);
    
    atm_fields_compact.data.resize(1,sp_p_grid.nelem(),sp_lat_grid.nelem(),sp_lon_grid.nelem());
    
    atm_fields_compact.data(0,joker,joker,joker) = field.data;
}



// Workspace method, doxygen header is auto-generated
// 2011-05-11 Gerrit Holl
void batch_atm_fields_compactAddConstant(// WS Output:
                                         ArrayOfGriddedField4& batch_atm_fields_compact,
                                         // WS Generic Input:
                                         const String& name,
                                         const Numeric& value,
                                         const Index& prepend,
                                         const ArrayOfString& condensibles,
                                         const Verbosity& verbosity)
{
    for (Index i=0; i<batch_atm_fields_compact.nelem(); i++)
    {
        atm_fields_compactAddConstant(batch_atm_fields_compact[i],
                                      name, value, prepend, condensibles, verbosity);
    }

}

// Workspace method, doxygen header is auto-generated
// 2011-05-09 Gerrit Holl
void batch_atm_fields_compactAddSpecies(// WS Output:
                                        ArrayOfGriddedField4& batch_atm_fields_compact,
                                        // WS Generic Input:
                                        const String& name,
                                        const GriddedField3& species,
                                        const Index& prepend,
                                        const Verbosity& verbosity)
{
    const Index nelem = batch_atm_fields_compact.nelem();

    String fail_msg;
    bool failed = false;

    // Parallelise this for-loop (some interpolation is being done, so it may
    // be beneficial)
#pragma omp parallel for      \
  if (!arts_omp_in_parallel()  \
      && nelem >= arts_omp_get_max_threads())
    for (Index i=0; i<nelem; i++)
    {
        try {
            atm_fields_compactAddSpecies(batch_atm_fields_compact[i], name, species, prepend, verbosity);
        }
        catch (runtime_error e)
        {
#pragma omp critical (batch_atm_fields_compactAddSpecies_fail)
            { fail_msg = e.what(); failed = true; }
        }
    }

    if (failed) throw runtime_error(fail_msg);
}

// Workspace method, doxygen header is auto-generated
// 2015-06-30 Jana Mendrok
void batch_atm_fields_compactCleanup (//WS Output:
                          ArrayOfGriddedField4& batch_atm_fields_compact,
                          //WS Input:
                          const Numeric& threshold,
                          const Verbosity& verbosity)
{
    for (Index i=0; i<batch_atm_fields_compact.nelem(); i++)
    {
        atm_fields_compactCleanup(batch_atm_fields_compact[i], threshold, verbosity );
    }
}



// Workspace method, doxygen header is auto-generated.
void batch_atm_fields_compactFromArrayOfMatrix(// WS Output:
                                               ArrayOfGriddedField4& batch_atm_fields_compact,
                                               // WS Input:
                                               const Index& atmosphere_dim,
                                               // WS Generic Input:
                                               const ArrayOfMatrix& am,
                                               // Control Parameters:
                                               const ArrayOfString& field_names,
                                               const Verbosity& verbosity)
{
  const Index amnelem = am.nelem();

  if (amnelem == 0)
    {
      ostringstream os; 
      os << "No elements in atmospheric scenario batch.\n"
         << "Check, whether any batch atmosphere file has been read!";
      throw runtime_error( os.str() );
    }

  // We use the existing WSMs atm_fields_compactFromMatrix and
  // atm_fields_compactAddConstant to do most of the work.

  // Make output variable the proper size:
  batch_atm_fields_compact.resize(amnelem);

  String fail_msg;
  bool failed = false;

  // Loop the batch cases:
#pragma omp parallel for      \
  if (!arts_omp_in_parallel()  \
      && amnelem >= arts_omp_get_max_threads())
  for (Index i=0; i<amnelem; ++i)
    {
      // Skip remaining iterations if an error occurred
      if (failed) continue;

      // All the input variables are visible here, despite the
      // "default(none)". The reason is that they are return by
      // reference arguments of this function, which are shared
      // automatically. 

      // The try block here is necessary to correctly handle
      // exceptions inside the parallel region. 
      try
        {
          atm_fields_compactFromMatrix(batch_atm_fields_compact[i],
                                       atmosphere_dim,
                                       am[i],
                                       field_names,
                                       verbosity);
        }
      catch (runtime_error e)
        {
#pragma omp critical (batch_atm_fields_compactFromArrayOfMatrix_fail)
            { fail_msg = e.what(); failed = true; }
        }
    }    

    if (failed) throw runtime_error(fail_msg);
}


void AtmFieldsFromCompact(// WS Output:
                          Vector& p_grid,
                          Vector& lat_grid,
                          Vector& lon_grid,
                          Tensor3& t_field,
                          Tensor3& z_field,
                          Tensor4& vmr_field,
                          Tensor4& scat_species_mass_density_field,
                          Tensor4& scat_species_mass_flux_field,
                          Tensor4& scat_species_number_density_field,
                          Tensor4& scat_species_mean_mass_field,
                          // WS Input:
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const ArrayOfString& scat_species,
                          const GriddedField4& atm_fields_compact,
                          const Index&  atmosphere_dim,
                          const String& delim,
                          const Numeric& p_min,
                          // Control parameters:
                          const Index& check_gridnames,
                          const Verbosity&)
{
  // Make a handle on atm_fields_compact to save typing:
  const GriddedField4& c = atm_fields_compact;
  
  // Check if the grids in our data match atmosphere_dim
  // (throws an error if the dimensionality is not correct):
  chk_atm_grids(atmosphere_dim,
                c.get_numeric_grid(GFIELD4_P_GRID),
                c.get_numeric_grid(GFIELD4_LAT_GRID),
                c.get_numeric_grid(GFIELD4_LON_GRID));

  // Optional check for gridnames.
  if (check_gridnames == 1)
  {
    chk_griddedfield_gridname(c, 1, "Pressure");
    chk_griddedfield_gridname(c, 2, "Latitude");
    chk_griddedfield_gridname(c, 3, "Longitude");
  }

  const Index nf   = c.get_grid_size(GFIELD4_FIELD_NAMES);
  const Index np   = c.get_grid_size(GFIELD4_P_GRID);
//  const Index nlat = c.get_grid_size(GFIELD4_LAT_GRID);
//  const Index nlon = c.get_grid_size(GFIELD4_LON_GRID);
  Index nlat = c.get_grid_size(GFIELD4_LAT_GRID);
  Index nlon = c.get_grid_size(GFIELD4_LON_GRID);
  if (nlat==0) nlat++;
  if (nlon==0) nlon++;

  // Grids:
  p_grid = c.get_numeric_grid(GFIELD4_P_GRID);
  lat_grid = c.get_numeric_grid(GFIELD4_LAT_GRID);
  lon_grid = c.get_numeric_grid(GFIELD4_LON_GRID);
  
  // Reduce p_grid to region below p_min
  Index l=np-1;
  bool search_toa=true;
  while( search_toa && l>0 )
    {
      if( p_grid[l-1] < p_min )
        l--;
      else
        search_toa=false;
    }
  if( search_toa )
    {
      ostringstream os;
      os << "At least one atmospheric level with pressure larger p_min (="
         << p_min << ")\n" << "is needed, but none is found.";
      throw runtime_error( os.str() );      
    }
  const Index npn = l+1;
  p_grid = p_grid[Range(0,npn)];

  const Index nsa = abs_species.nelem();
  const Index nsp = scat_species.nelem();
  
  // Check that there is at least one VMR species:
  if (nsa<1)
    {
      ostringstream os;
      os << "There must be at least one absorption species.";
      throw runtime_error( os.str() );
    }

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
  const String as_type="abs_species";
  const String ss_type="scat_species";

  // Find temperature field:
  found = false;
  t_field.resize(npn,nlat,nlon);
  for (Index i=0; i<nf; ++i)
    {
      if (c.get_string_grid(GFIELD4_FIELD_NAMES)[i] == "T")
        {
          if (found)
            {
               ostringstream os;
               os << "Only one temperature ('T') field allowed,\n"
                  << "but found at least 2.";
               throw runtime_error( os.str() );
            }
          else
            {
              found = true;
              t_field = c.data(i,Range(0,npn),Range(joker),Range(joker));
            }
        }
    }
  if (!found)
    {
      ostringstream os;
      os << "One temperature ('T') field required, but none found";
      throw runtime_error( os.str() );
    }

  // Find Altitude field:
  found = false;
  z_field.resize(npn,nlat,nlon);
  for (Index i=0; i<nf; ++i)
    {
      if (c.get_string_grid(GFIELD4_FIELD_NAMES)[i] == "z")
        {
          if (found)
            {
               ostringstream os;
               os << "Only one altitude ('z') field allowed,\n"
                  << "but found at least 2.";
               throw runtime_error( os.str() );
            }
          else
            {
              found = true;
              z_field = c.data(i,Range(0,npn),Range(joker),Range(joker));
            }
        }
    }
  if (!found)
    {
      ostringstream os;
      os << "One altitude ('z') field required, but none found";
      throw runtime_error( os.str() );
    }

  // Extracting the required abs_species fields:
  vmr_field.resize(nsa,npn,nlat,nlon);
  for (Index j=0; j<nsa; ++j)
    {
      using global_data::species_data;  // The species lookup data:
      const String as_name = species_data[abs_species[j][0].Species()].Name();
      found = false;
      Index i = 0;
      String species_type;
      String species_name;
      while (!found && i<nf)
        {
          parse_atmcompact_speciestype(
            species_type, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);
          // do we have an abs_species type field?
          if (species_type==as_type)
            {
              parse_atmcompact_speciesname(
                species_name, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);
              if (species_name==as_name)
                {
                  found=true;
                  vmr_field(j,Range(joker),Range(joker),Range(joker)) =
                    c.data(i,Range(0,npn),Range(joker),Range(joker));
                }
            }
          i++;
        }
      if (!found)
        {
          ostringstream os;
          os << "No field for absorption species '" << as_name << "' found.";
          throw runtime_error( os.str() );
        }
    }
 
  // Extracting the required scattering species fields:
  scat_species_mass_density_field.resize(nsp,npn,nlat,nlon);
  scat_species_mass_density_field = NAN;
  scat_species_mass_flux_field.resize(nsp,npn,nlat,nlon);
  scat_species_mass_flux_field = NAN;
  scat_species_number_density_field.resize(nsp,npn,nlat,nlon);
  scat_species_number_density_field = NAN;
  scat_species_mean_mass_field.resize(nsp,npn,nlat,nlon);
  scat_species_mean_mass_field = NAN;

  // to be on safe said (avoiding overwriting of fields), we separately monitor
  // found status of the individual fields possible per scat_species (we check
  // each of them whether present, and if sort it into the scat_species**fields.
  // but for each of them, we only take the first match.).
  bool found_md;
  bool found_mf;
  bool found_nd;
  bool found_mm;

  const String md="mass_density";
  const String mf="mass_flux";
  const String nd="number_density";
  const String mm="mean_mass";

  for (Index j=0; j<nsp; ++j)
    {
      String ss_name;
      parse_partfield_name(ss_name, scat_species[j], delim);
      found = false;
      found_md = false;
      found_mf = false;
      found_nd = false;
      found_mm = false;
      Index i = 0;
      String species_type;
      String species_name;
      String scat_type;
      while (!found && i<nf)
        {
          parse_atmcompact_speciestype(
            species_type, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);
          // do we have an scat_species type field?
          if (species_type==ss_type)
            {
              parse_atmcompact_speciesname(
                species_name, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);
              if (species_name==ss_name)
                {
                  parse_atmcompact_scattype(
                    scat_type, c.get_string_grid(GFIELD4_FIELD_NAMES)[i], delim);
                  if (!found_md && scat_type==md)
                    {
                      found_md=true;
                      scat_species_mass_density_field(j,
                        Range(joker),Range(joker),Range(joker))
                        = c.data(i,Range(0,npn),Range(joker),Range(joker));
                    }
                  else if (!found_mf && scat_type==mf)
                    {
                      found_mf=true;
                      scat_species_mass_flux_field(j,
                        Range(joker),Range(joker),Range(joker))
                        = c.data(i,Range(0,npn),Range(joker),Range(joker));
                    }
                  else if (!found_nd && scat_type==nd)
                    {
                      found_nd=true;
                      scat_species_number_density_field(j,
                        Range(joker),Range(joker),Range(joker))
                        = c.data(i,Range(0,npn),Range(joker),Range(joker));
                    }
                  else if (!found_mm && scat_type==mm)
                    {
                      found_mm=true;
                      scat_species_mean_mass_field(j,
                        Range(joker),Range(joker),Range(joker))
                        = c.data(i,Range(0,npn),Range(joker),Range(joker));
                    }
                  else
                    {
                      ostringstream os;
                      os << "scat_species field type " << scat_type
                         << "is unknown, hence can't be handled.";
                      throw runtime_error( os.str() );
                    }
                }
            }
          // if all possible scat_species fields have been filled, we don't need
          // to check the rest any longer.
          if (found_md && found_mf && found_nd && found_mm) found=true;
          i++;
        }
      // did we find at least one of the possible scat_species fields for this
      // scat species?
      found = (found_md || found_mf || found_nd || found_mm);
      if (!found)
        {
          ostringstream os;
          os << "No field for scattering species '" << ss_name << "' found.";
          throw runtime_error( os.str() );
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmosphereSet1D(// WS Output:
                     Index&    atmosphere_dim,
                     Vector&   lat_grid,
                     Vector&   lon_grid,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  out2 << "  Sets the atmospheric dimensionality to 1.\n";
  out3 << "    atmosphere_dim = 1\n";
  out3 << "    lat_grid is set to be an empty vector\n";
  out3 << "    lon_grid is set to be an empty vector\n";
  
  atmosphere_dim = 1;
  lat_grid.resize(0);
  lon_grid.resize(0);
 
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmosphereSet2D(// WS Output:
                     Index&    atmosphere_dim,
                     Vector&   lon_grid,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  out2 << "  Sets the atmospheric dimensionality to 2.\n";
  out3 << "    atmosphere_dim = 2\n";
  out3 << "    lon_grid is set to be an empty vector\n";
  
  atmosphere_dim = 2;
  lon_grid.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmosphereSet3D(// WS Output:
                     Index&    atmosphere_dim,
                     Vector&   lat_true,
                     Vector&   lon_true,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  out2 << "  Sets the atmospheric dimensionality to 3.\n";
  out3 << "    atmosphere_dim = 3\n";
  
  atmosphere_dim = 3;
  lat_true.resize(0);
  lon_true.resize(0);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsCalc(//WS Output:
                   Tensor3& t_field,
                   Tensor3& z_field,
                   Tensor4& vmr_field,
                   Tensor4& t_nlte_field,
                   //WS Input
                   const Vector&         p_grid,
                   const Vector&         lat_grid,
                   const Vector&         lon_grid,
                   const GriddedField3&        t_field_raw,
                   const GriddedField3&        z_field_raw,
                   const ArrayOfGriddedField3& vmr_field_raw,
                   const ArrayOfGriddedField3& t_nlte_field_raw,
                   const Index&          atmosphere_dim,
                   // WS Generic Input:
                   const Index& interp_order,
                   const Index& vmr_zeropadding,
                   const Index& vmr_nonegative,
                   const Index& nlte_when_negative,
                   const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  const ConstVectorView tfr_p_grid =
    t_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView tfr_lat_grid =
    t_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView tfr_lon_grid =
    t_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
  const ConstVectorView zfr_p_grid =
    z_field_raw.get_numeric_grid(GFIELD3_P_GRID);
  const ConstVectorView zfr_lat_grid =
    z_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
  const ConstVectorView zfr_lon_grid =
    z_field_raw.get_numeric_grid(GFIELD3_LON_GRID);

  out2 << "  Interpolation order: " << interp_order << "\n";
  
  // Basic checks of input variables
  //
  // Atmosphere
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  

  //==========================================================================
  if ( atmosphere_dim == 1)
    {
      if( !( tfr_lat_grid.nelem() == 1 &&
             tfr_lon_grid.nelem() == 1 ))
        throw runtime_error(
                            "Temperature data (T_field) has wrong dimension "
                            "(2D or 3D).\n"
                            );

      if( !( zfr_lat_grid.nelem() == 1 &&
             zfr_lon_grid.nelem() == 1 ))
        throw runtime_error(
                            "Altitude data (z_field) has wrong dimension "
                            "(2D or 3D).\n"
                            );

      GriddedField3 temp_gfield3;

      GriddedFieldPRegrid(temp_gfield3, p_grid, t_field_raw,
                          interp_order, 0, verbosity);
      t_field = temp_gfield3.data;

      GriddedFieldPRegrid(temp_gfield3, p_grid, z_field_raw,
                          interp_order, 0, verbosity);
      z_field = temp_gfield3.data;

      ArrayOfGriddedField3 temp_agfield3;
      try {
        GriddedFieldPRegrid(temp_agfield3, p_grid, vmr_field_raw, interp_order,
                            vmr_zeropadding, verbosity);
      } catch (runtime_error e) {
          ostringstream os;
          os << e.what() << "\n"
             << "Note that you can explicitly set vmr_zeropadding "
             << "to 1 in the method call.";
          throw runtime_error(os.str());

      }
      FieldFromGriddedField(vmr_field, p_grid, lat_grid, lon_grid, 
                            temp_agfield3, verbosity);
      
      // Non-LTE interpolation
      if(t_nlte_field_raw.nelem())
      {
        GriddedFieldPRegrid(temp_agfield3, p_grid, t_nlte_field_raw, 
                            interp_order, 0, verbosity);
        FieldFromGriddedField(t_nlte_field, p_grid, lat_grid, lon_grid,
                              temp_agfield3, verbosity);
      }
      
    }

  //=========================================================================
  else if(atmosphere_dim == 2)
    {
      if( tfr_lat_grid.nelem() == 1 &&
          tfr_lon_grid.nelem() == 1 )
        throw runtime_error(
                            "Raw data has wrong dimension (1D). "
                            "You have to use \n"
                            "AtmFieldsCalcExpand1D instead of AtmFieldsCalc."
                            );
      
      //Resize variables
      t_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
      z_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
      vmr_field.resize(vmr_field_raw.nelem(), p_grid.nelem(), lat_grid.nelem(),
                       1);
      if(t_nlte_field_raw.nelem())
        t_nlte_field.resize(t_nlte_field_raw.nelem(),
                            p_grid.nelem(), lat_grid.nelem(),1);
      else 
        t_nlte_field.resize(0,0,0,0);
      
      
      // Gridpositions:
      ArrayOfGridPosPoly gp_p(p_grid.nelem());
      ArrayOfGridPosPoly gp_lat(lat_grid.nelem());
      
      // Interpolate t_field:
      
      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_pgrids("Raw temperature to p_grid, 2D case",
                               tfr_p_grid,
                               p_grid,
                               interp_order);
      chk_interpolation_grids("Raw temperature to lat_grid, 2D case",
                              tfr_lat_grid,
                              lat_grid,
                              interp_order);

      // Calculate grid positions:
      p2gridpos_poly( gp_p, tfr_p_grid, p_grid, interp_order );
      gridpos_poly( gp_lat, tfr_lat_grid, lat_grid, interp_order );
            
      // Interpolation weights:
      Tensor3 itw(p_grid.nelem(), lat_grid.nelem(),
                  (interp_order+1)*(interp_order+1));
      // (4 interpolation weights are required for example for linear 2D
      // interpolation)
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( t_field(joker, joker, 0 ), itw,
              t_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
      
      
      // Interpolate z_field:

      // Check that interpolation grids are ok (and throw a detailed
      // error message if not):
      chk_interpolation_pgrids("Raw z to p_grid, 2D case",
                               zfr_p_grid,
                               p_grid,
                               interp_order);
      chk_interpolation_grids("Raw z to lat_grid, 2D case",
                              zfr_lat_grid,
                              lat_grid,
                              interp_order);

      // Calculate grid positions:
      p2gridpos_poly( gp_p, zfr_p_grid, p_grid, interp_order );
      gridpos_poly( gp_lat, zfr_lat_grid, lat_grid, interp_order );
            
      // Interpolation weights:
      interpweights( itw, gp_p, gp_lat);
      
      // Interpolate:
      interp( z_field(joker, joker, 0), itw, 
              z_field_raw.data(joker, joker, 0), gp_p, gp_lat);
      
      
      // Interpolate vmr_field. 
      // Loop over the gaseous species:
      for (Index gas_i = 0; gas_i < vmr_field_raw.nelem(); gas_i++)
        {
          ostringstream os; 

          if( !( 
            vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID).nelem() != 1 &&
            vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LON_GRID).nelem() == 1 ))
            {
              os << "VMR data of the " << gas_i << " the species has "
                 << "wrong dimension (1D or 3D). \n";
              throw runtime_error( os.str() );
            }

          // Check that interpolation grids are ok (and throw a detailed
          // error message if not):
          os << "Raw VMR[" << gas_i << "] to p_grid, 2D case";
          chk_interpolation_pgrids(os.str(),
            vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID),
            p_grid, interp_order);
          os.str("");
          os << "Raw VMR[" << gas_i << "] to lat_grid, 2D case";
          chk_interpolation_grids(os.str(),
            vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID),
            lat_grid, interp_order);

          // Calculate grid positions:
          p2gridpos_poly(gp_p,
                         vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_P_GRID), 
                         p_grid, interp_order);
          gridpos_poly(gp_lat,
                       vmr_field_raw[gas_i].get_numeric_grid(GFIELD3_LAT_GRID), 
                       lat_grid, interp_order);
                  
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat);
          
          // Interpolate:
          interp( vmr_field(gas_i, joker, joker, 0),
                  itw, vmr_field_raw[gas_i].data(joker, joker, 0),
                  gp_p, gp_lat);
        }
        
        // Interpolat Non-LTE
        for (Index qi_i = 0; qi_i < t_nlte_field_raw.nelem(); qi_i++)
        {
          ostringstream os; 

          if( !( 
            t_nlte_field_raw[qi_i].get_numeric_grid(GFIELD3_LAT_GRID).nelem() != 1 &&
            t_nlte_field_raw[qi_i].get_numeric_grid(GFIELD3_LON_GRID).nelem() == 1 ))
            {
              os << "NLTE data of the " << qi_i << " temperature field has "
                 << "wrong dimension (1D or 3D). \n";
              throw std::runtime_error( os.str() );
            }

          // Check that interpolation grids are ok (and throw a detailed
          // error message if not):
          os << "Raw NLTE[" << qi_i << "] to p_grid, 2D case";
          chk_interpolation_pgrids(os.str(),
            t_nlte_field_raw[qi_i].get_numeric_grid(GFIELD3_P_GRID),
            p_grid, interp_order);
          os.str("");
          os << "Raw NLTE[" << qi_i << "] to lat_grid, 2D case";
          chk_interpolation_grids(os.str(),
            t_nlte_field_raw[qi_i].get_numeric_grid(GFIELD3_LAT_GRID),
            lat_grid, interp_order);

          // Calculate grid positions:
          p2gridpos_poly(gp_p, 
                         t_nlte_field_raw[qi_i].get_numeric_grid(GFIELD3_P_GRID), 
                         p_grid, interp_order);
          gridpos_poly(gp_lat,
                       t_nlte_field_raw[qi_i].get_numeric_grid(GFIELD3_LAT_GRID), 
                       lat_grid, interp_order);
                  
          // Interpolation weights:
          interpweights( itw, gp_p, gp_lat);
          
          // Interpolate:
          interp( t_nlte_field(qi_i, joker, joker, 0),
                  itw, t_nlte_field_raw[qi_i].data(joker, joker, 0),
                  gp_p, gp_lat);
        }
    }

  //================================================================
  // atmosphere_dim = 3    
  else if(atmosphere_dim == 3)
    {
      if( tfr_lat_grid.nelem() == 1 &&
          tfr_lon_grid.nelem() == 1 )
        throw runtime_error(
                            "Raw data has wrong dimension. You have to use \n"
                            "AtmFieldsCalcExpand1D instead of AtmFieldsCalc."
                            );

      GriddedField3 temp_gfield3;

      GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                               t_field_raw, interp_order, verbosity);
      GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3,
                          interp_order, 0, verbosity);
      t_field = temp_gfield3.data;

      GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                               z_field_raw, interp_order, verbosity);
      GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3,
                          interp_order, 0, verbosity);
      z_field = temp_gfield3.data;

      ArrayOfGriddedField3 temp_agfield3;
      GriddedFieldLatLonRegrid(temp_agfield3, lat_grid, lon_grid,
                               vmr_field_raw, interp_order, verbosity);
      try {
        GriddedFieldPRegrid(temp_agfield3, p_grid, temp_agfield3, interp_order,
                            vmr_zeropadding, verbosity);
      } catch (runtime_error e) {
          ostringstream os;
          os << e.what() << "\n"
             << "Note that you can explicitly set vmr_zeropadding "
             << "to 1 in the method call.";
          throw runtime_error(os.str());

      }
      FieldFromGriddedField(vmr_field, p_grid, lat_grid, lon_grid,
                           temp_agfield3, verbosity);
      
      if(t_nlte_field_raw.nelem())
      {
        GriddedFieldLatLonRegrid(temp_agfield3, lat_grid, lon_grid,
                                 t_nlte_field_raw, interp_order, verbosity);
        GriddedFieldPRegrid(temp_agfield3, p_grid, temp_agfield3,
                            interp_order, 0, verbosity);
      }
    }
  else
  {
    // We can never get here, since there was a runtime 
    // error check for atmosphere_dim at the beginning.
    assert(false);
  }

  // remove negatives?
  if( vmr_nonegative )
    {
      for( Index ib=0; ib<vmr_field.nbooks(); ib++ )
        {
          for( Index ip=0; ip<vmr_field.npages(); ip++ )
            {
              for( Index ir=0; ir<vmr_field.nrows(); ir++ )
                {
                  for( Index ic=0; ic<vmr_field.ncols(); ic++ )
                    {
                      if( vmr_field(ib,ip,ir,ic) < 0 )
                        { vmr_field(ib,ip,ir,ic) = 0; }
        }   }   }   } 
    }
    
  // what to do with negative nlte temperatures?
  if(nlte_when_negative!=-1) // Nothing
    if(t_nlte_field_raw.nelem())
      for( Index ib=0; ib<t_nlte_field.nbooks(); ib++ )
        for( Index ip=0; ip<t_nlte_field.npages(); ip++ )
          for( Index ir=0; ir<t_nlte_field.nrows(); ir++ )
            for( Index ic=0; ic<t_nlte_field.ncols(); ic++ )
              if( t_nlte_field(ib,ip,ir,ic) < 0 )
                // Set to atmospheric temperature or to nil.
                t_nlte_field(ib,ip,ir,ic) = 
                  nlte_when_negative==1?t_field(ip,ir,ic):0; 
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MagFieldsCalc( //WS Output:
                    Tensor3& mag_u_field,
                    Tensor3& mag_v_field,
                    Tensor3& mag_w_field,
                    //WS Input
                    const Vector&         p_grid,
                    const Vector&         lat_grid,
                    const Vector&         lon_grid,
                    const GriddedField3&        mag_u_field_raw,
                    const GriddedField3&        mag_v_field_raw,
                    const GriddedField3&        mag_w_field_raw,
                    const Index&          atmosphere_dim,
                    // WS Generic Input:
                    const Index& interp_order,
                    const Verbosity& verbosity)
{
    CREATE_OUT2;
    
    const ConstVectorView ufr_p_grid =
      mag_u_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    const ConstVectorView ufr_lat_grid =
      mag_u_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
    const ConstVectorView ufr_lon_grid =
      mag_u_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
    const ConstVectorView vfr_p_grid =
      mag_v_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    const ConstVectorView vfr_lat_grid =
      mag_v_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
    const ConstVectorView vfr_lon_grid =
      mag_v_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
    const ConstVectorView wfr_p_grid =
      mag_w_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    const ConstVectorView wfr_lat_grid =
      mag_w_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
    const ConstVectorView wfr_lon_grid =
      mag_w_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
    
    out2 << "  Interpolation order: " << interp_order << "\n";
    
    // Basic checks of input variables
    //
    // Atmosphere
    chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
    chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
    
    
    //==========================================================================
    if ( atmosphere_dim == 1)
    {
        if( !( ufr_lat_grid.nelem() == 1 && ufr_lon_grid.nelem() == 1 ))
            throw std::runtime_error(
            "Magnetic u field data has wrong dimension (2D or 3D).\n");
        if( !( vfr_lat_grid.nelem() == 1 && vfr_lon_grid.nelem() == 1 ))
            throw std::runtime_error(
            "Magnetic v field data has wrong dimension (2D or 3D).\n");
        if( !( wfr_lat_grid.nelem() == 1 && wfr_lon_grid.nelem() == 1 ))
            throw std::runtime_error(
            "Magnetic w field data has wrong dimension  (2D or 3D).\n");
        
        GriddedField3 temp_gfield3;
        
        GriddedFieldPRegrid(temp_gfield3, p_grid, mag_u_field_raw,
                            interp_order, 0, verbosity);
        mag_u_field = temp_gfield3.data;
        
        GriddedFieldPRegrid(temp_gfield3, p_grid, mag_v_field_raw,
                            interp_order, 0, verbosity);
        mag_v_field = temp_gfield3.data;
        
        GriddedFieldPRegrid(temp_gfield3, p_grid, mag_w_field_raw,
                            interp_order, 0, verbosity);
        mag_w_field = temp_gfield3.data;
    }
    
    //=========================================================================
    else if(atmosphere_dim == 2)
    {
        if( ufr_lat_grid.nelem() == 1 && ufr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension (1D). You have to use \n"
            "MagFieldsCalcExpand1D instead of MagFieldsCalc.");
        if( vfr_lat_grid.nelem() == 1 && vfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension (1D). You have to use \n"
            "MagFieldsCalcExpand1D instead of MagFieldsCalc.");
        if( wfr_lat_grid.nelem() == 1 && wfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension (1D). You have to use \n"
            "MagFieldsCalcExpand1D instead of MagFieldsCalc.");
        
        //Resize variables
        mag_u_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
        mag_v_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
        mag_w_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
        
        // Gridpositions:
        ArrayOfGridPosPoly gp_p(p_grid.nelem());
        ArrayOfGridPosPoly gp_lat(lat_grid.nelem());
        
        // Interpolate mag_u_field:
        
        // Check that interpolation grids are ok (and throw a detailed
        // error message if not):
        chk_interpolation_pgrids("Raw u field to p_grid, 2D case",
                                 ufr_p_grid,p_grid,interp_order);
        chk_interpolation_grids("Raw u field to lat_grid, 2D case",
                                ufr_lat_grid,lat_grid,interp_order);
        
        // Calculate grid positions:
        p2gridpos_poly( gp_p, ufr_p_grid, p_grid, interp_order );
        gridpos_poly( gp_lat, ufr_lat_grid, lat_grid, interp_order );
        
        // Interpolation weights:
        Tensor3 itwu(p_grid.nelem(), lat_grid.nelem(),
                     (interp_order+1)*(interp_order+1));
        // (4 interpolation weights are required for example for linear 2D
        // interpolation)
        interpweights( itwu, gp_p, gp_lat);
        
        // Interpolate:
        interp( mag_u_field(joker, joker, 0 ), itwu,
                mag_u_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
        
        // Interpolate mag_v_field:
        
        // Check that interpolation grids are ok (and throw a detailed
        // error message if not):
        chk_interpolation_pgrids("Raw v field to p_grid, 2D case",
                                 vfr_p_grid,p_grid,interp_order);
        chk_interpolation_grids("Raw v field to lat_grid, 2D case",
                                vfr_lat_grid,lat_grid,interp_order);
        
        // Calculate grid positions:
        p2gridpos_poly( gp_p, vfr_p_grid, p_grid, interp_order );
        gridpos_poly( gp_lat, vfr_lat_grid, lat_grid, interp_order );
        
        // Interpolation weights:
        Tensor3 itwv(p_grid.nelem(), lat_grid.nelem(),
                     (interp_order+1)*(interp_order+1));
        // (4 interpolation weights are required for example for linear 2D
        // interpolation)
        interpweights( itwv, gp_p, gp_lat);
        
        // Interpolate:
        interp( mag_v_field(joker, joker, 0 ), itwv,
                mag_v_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
        
        // Interpolate mag_w_field:
        
        // Check that interpolation grids are ok (and throw a detailed
        // error message if not):
        chk_interpolation_pgrids("Raw w field to p_grid, 2D case",
                                 wfr_p_grid,p_grid,interp_order);
        chk_interpolation_grids("Raw w field to lat_grid, 2D case",
                                wfr_lat_grid,lat_grid,interp_order);
        
        // Calculate grid positions:
        p2gridpos_poly( gp_p, wfr_p_grid, p_grid, interp_order );
        gridpos_poly( gp_lat, wfr_lat_grid, lat_grid, interp_order );
        
        // Interpolation weights:
        Tensor3 itww(p_grid.nelem(), lat_grid.nelem(),
                     (interp_order+1)*(interp_order+1));
        // (4 interpolation weights are required for example for linear 2D
        // interpolation)
        interpweights( itww, gp_p, gp_lat);
        
        // Interpolate:
        interp( mag_w_field(joker, joker, 0 ), itww,
                mag_w_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
    }
    
    //================================================================
    // atmosphere_dim = 3    
    else if(atmosphere_dim == 3)
    {
        if( ufr_lat_grid.nelem() == 1 && ufr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension. You have to use \n"
            "MagFieldsCalcExpand1D instead of MagFieldsCalc.");
        if( vfr_lat_grid.nelem() == 1 && vfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension. You have to use \n"
            "MagFieldsCalcExpand1D instead of MagFieldsCalc.");
        if( wfr_lat_grid.nelem() == 1 && wfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension. You have to use \n"
            "MagFieldsCalcExpand1D instead of MagFieldsCalc.");
        
        GriddedField3 temp_gfield3;
        
        GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                                 mag_u_field_raw, interp_order, verbosity);
        GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3,
                            interp_order, 0, verbosity);
        mag_u_field = temp_gfield3.data;
        
        GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                                 mag_v_field_raw, interp_order, verbosity);
        GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3,
                            interp_order, 0, verbosity);
        mag_v_field = temp_gfield3.data;
        
        GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                                 mag_w_field_raw, interp_order, verbosity);
        GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3,
                            interp_order, 0, verbosity);
        mag_w_field = temp_gfield3.data;
        
        
    }
    else
    {
        // We can never get here, since there was a runtime 
        // error check for atmosphere_dim at the beginning.
        assert(false);
    }
    
}


/* Workspace method: Doxygen documentation will be auto-generated */
void WindFieldsCalc(//WS Output:
                    Tensor3& wind_u_field,
                    Tensor3& wind_v_field,
                    Tensor3& wind_w_field,
                    //WS Input
                    const Vector&         p_grid,
                    const Vector&         lat_grid,
                    const Vector&         lon_grid,
                    const GriddedField3&        wind_u_field_raw,
                    const GriddedField3&        wind_v_field_raw,
                    const GriddedField3&        wind_w_field_raw,
                    const Index&          atmosphere_dim,
                    // WS Generic Input:
                    const Index& interp_order,
                    const Verbosity& verbosity)
{
    CREATE_OUT2;
    
    const ConstVectorView ufr_p_grid =
      wind_u_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    const ConstVectorView ufr_lat_grid = 
      wind_u_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
    const ConstVectorView ufr_lon_grid =
      wind_u_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
    const ConstVectorView vfr_p_grid =
      wind_v_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    const ConstVectorView vfr_lat_grid =
      wind_v_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
    const ConstVectorView vfr_lon_grid =
      wind_v_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
    const ConstVectorView wfr_p_grid =
      wind_w_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    const ConstVectorView wfr_lat_grid =
       wind_w_field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
    const ConstVectorView wfr_lon_grid =
       wind_w_field_raw.get_numeric_grid(GFIELD3_LON_GRID);
    
    out2 << "  Interpolation order: " << interp_order << "\n";
    
    // Basic checks of input variables
    //
    // Atmosphere
    chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
    chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
    
    
    //==========================================================================
    if ( atmosphere_dim == 1)
    {
        if( !( ufr_lat_grid.nelem() == 1 && ufr_lon_grid.nelem() == 1 ))
            throw std::runtime_error(
            "Wind u field data has wrong dimension (2D or 3D).\n");
        if( !( vfr_lat_grid.nelem() == 1 && vfr_lon_grid.nelem() == 1 ))
            throw std::runtime_error(
            "Wind v field data has wrong dimension (2D or 3D).\n");
        if( !( wfr_lat_grid.nelem() == 1 && wfr_lon_grid.nelem() == 1 ))
            throw std::runtime_error(
            "Wind w field data has wrong dimension  (2D or 3D).\n");
        
        GriddedField3 temp_gfield3;
        
        GriddedFieldPRegrid(temp_gfield3, p_grid, wind_u_field_raw,
                            interp_order, 0, verbosity);
        wind_u_field = temp_gfield3.data;
        
        GriddedFieldPRegrid(temp_gfield3, p_grid, wind_v_field_raw,
                            interp_order, 0, verbosity);
        wind_v_field = temp_gfield3.data;
        
        GriddedFieldPRegrid(temp_gfield3, p_grid, wind_w_field_raw,
                            interp_order, 0, verbosity);
        wind_w_field = temp_gfield3.data;
    }
    
    //=========================================================================
    else if(atmosphere_dim == 2)
    {
        if( ufr_lat_grid.nelem() == 1 && ufr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension (1D). You have to use \n"
            "WindFieldsCalcExpand1D instead of WindFieldsCalc.");
        if( vfr_lat_grid.nelem() == 1 && vfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension (1D). You have to use \n"
            "WindFieldsCalcExpand1D instead of WindFieldsCalc.");
        if( wfr_lat_grid.nelem() == 1 && wfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension (1D). You have to use \n"
            "WindFieldsCalcExpand1D instead of WindFieldsCalc.");
        
        //Resize variables
        wind_u_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
        wind_v_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
        wind_w_field.resize(p_grid.nelem(), lat_grid.nelem(), 1);
        
        // Gridpositions:
        ArrayOfGridPosPoly gp_p(p_grid.nelem());
        ArrayOfGridPosPoly gp_lat(lat_grid.nelem());
        
        // Interpolate wind_u_field:
        
        // Check that interpolation grids are ok (and throw a detailed
        // error message if not):
        chk_interpolation_pgrids("Raw u field to p_grid, 2D case",
                                 ufr_p_grid,p_grid,interp_order);
        chk_interpolation_grids("Raw u field to lat_grid, 2D case",
                                ufr_lat_grid,lat_grid,interp_order);
        
        // Calculate grid positions:
        p2gridpos_poly( gp_p, ufr_p_grid, p_grid, interp_order );
        gridpos_poly( gp_lat, ufr_lat_grid, lat_grid, interp_order );
        
        // Interpolation weights:
        Tensor3 itwu(p_grid.nelem(), lat_grid.nelem(),
                     (interp_order+1)*(interp_order+1));
        // (4 interpolation weights are required for example for linear 2D
        // interpolation)
        interpweights( itwu, gp_p, gp_lat);
        
        // Interpolate:
        interp( wind_u_field(joker, joker, 0 ), itwu,
                wind_u_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
        
        // Interpolate wind_v_field:
        
        // Check that interpolation grids are ok (and throw a detailed
        // error message if not):
        chk_interpolation_pgrids("Raw v field to p_grid, 2D case",
                                 vfr_p_grid,p_grid,interp_order);
        chk_interpolation_grids("Raw v field to lat_grid, 2D case",
                                vfr_lat_grid,lat_grid,interp_order);
        
        // Calculate grid positions:
        p2gridpos_poly( gp_p, vfr_p_grid, p_grid, interp_order );
        gridpos_poly( gp_lat, vfr_lat_grid, lat_grid, interp_order );
        
        // Interpolation weights:
        Tensor3 itwv(p_grid.nelem(), lat_grid.nelem(),
                     (interp_order+1)*(interp_order+1));
        // (4 interpolation weights are required for example for linear 2D
        // interpolation)
        interpweights( itwv, gp_p, gp_lat);
        
        // Interpolate:
        interp( wind_v_field(joker, joker, 0 ), itwv,
                wind_v_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
        
        // Interpolate wind_w_field:
        
        // Check that interpolation grids are ok (and throw a detailed
        // error message if not):
        chk_interpolation_pgrids("Raw w field to p_grid, 2D case",
                                 wfr_p_grid,p_grid,interp_order);
        chk_interpolation_grids("Raw w field to lat_grid, 2D case",
                                wfr_lat_grid,lat_grid,interp_order);
        
        // Calculate grid positions:
        p2gridpos_poly( gp_p, wfr_p_grid, p_grid, interp_order );
        gridpos_poly( gp_lat, wfr_lat_grid, lat_grid, interp_order );
        
        // Interpolation weights:
        Tensor3 itww(p_grid.nelem(), lat_grid.nelem(),
                     (interp_order+1)*(interp_order+1));
        // (4 interpolation weights are required for example for linear 2D
        // interpolation)
        interpweights( itww, gp_p, gp_lat);
        
        // Interpolate:
        interp( wind_w_field(joker, joker, 0 ), itww,
                wind_w_field_raw.data(joker, joker, 0),  gp_p, gp_lat);
    }
    
    //================================================================
    // atmosphere_dim = 3    
    else if(atmosphere_dim == 3)
    {
        if( ufr_lat_grid.nelem() == 1 && ufr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension. You have to use \n"
            "WindFieldsCalcExpand1D instead of WindFieldsCalc.");
        if( vfr_lat_grid.nelem() == 1 && vfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension. You have to use \n"
            "WindFieldsCalcExpand1D instead of WindFieldsCalc.");
        if( wfr_lat_grid.nelem() == 1 && wfr_lon_grid.nelem() == 1 )
            throw std::runtime_error(
            "Raw data has wrong dimension. You have to use \n"
            "WindFieldsCalcExpand1D instead of WindFieldsCalc.");
        
        GriddedField3 temp_gfield3;
        
        GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                                 wind_u_field_raw, interp_order, verbosity);
        GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3,
                            interp_order, 0, verbosity);
        wind_u_field = temp_gfield3.data;
        
        GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                                 wind_v_field_raw, interp_order, verbosity);
        GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3, 
                            interp_order, 0, verbosity);
        wind_v_field = temp_gfield3.data;
        
        GriddedFieldLatLonRegrid(temp_gfield3, lat_grid, lon_grid,
                                 wind_w_field_raw, interp_order, verbosity);
        GriddedFieldPRegrid(temp_gfield3, p_grid, temp_gfield3,
                            interp_order, 0, verbosity);
        wind_w_field = temp_gfield3.data;
        
        
    }
    else
    {
        // We can never get here, since there was a runtime 
        // error check for atmosphere_dim at the beginning.
        assert(false);
    }
    
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsCalcExpand1D(Tensor3&              t_field,
                           Tensor3&              z_field,
                           Tensor4&              vmr_field,
                           Tensor4&              t_nlte_field,
                           const Vector&         p_grid,
                           const Vector&         lat_grid,
                           const Vector&         lon_grid,
                           const GriddedField3&  t_field_raw,
                           const GriddedField3&  z_field_raw,
                           const ArrayOfGriddedField3& vmr_field_raw,
                           const ArrayOfGriddedField3& t_nlte_field_raw,
                           const Index&          atmosphere_dim,
                           const Index&          interp_order,
                           const Index&          vmr_zeropadding,
                           const Index&          vmr_nonegative,
                           const Index&          nlte_when_negative,
                           const Verbosity&      verbosity)
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  if( atmosphere_dim == 1 )
    { throw runtime_error( 
     "This function is intended for 2D and 3D. For 1D, use *AtmFieldsCalc*.");}

  // Make 1D interpolation using some dummy variables
  Vector    vempty(0);
  Tensor3   t_temp, z_temp;
  Tensor4   vmr_temp, nlte_temp;
  AtmFieldsCalc(t_temp, z_temp, vmr_temp, nlte_temp, p_grid, vempty, vempty, 
                t_field_raw, z_field_raw, vmr_field_raw, t_nlte_field_raw, 1, interp_order,
                vmr_zeropadding, vmr_nonegative, nlte_when_negative, verbosity);

  // Move values from the temporary tensors to the return arguments
  const Index   np = p_grid.nelem();
  const Index   nlat = lat_grid.nelem();
        Index   nlon = lon_grid.nelem();
  if( atmosphere_dim == 2 )
    { nlon = 1; }
  const Index   nspecies = vmr_temp.nbooks();
  //
  assert( t_temp.npages() == np );
  //
  t_field.resize( np, nlat, nlon );
  z_field.resize( np, nlat, nlon );
  vmr_field.resize( nspecies, np, nlat, nlon );
  if(t_nlte_field_raw.nelem())
    t_nlte_field.resize(t_nlte_field_raw.nelem(),np,nlat,nlon);
  else
    t_nlte_field.resize(0,0,0,0);
  //
  for( Index ilon=0; ilon<nlon; ilon++ )
    {
      for( Index ilat=0; ilat<nlat; ilat++ )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              t_field(ip,ilat,ilon) = t_temp(ip,0,0);
              z_field(ip,ilat,ilon) = z_temp(ip,0,0);
              for( Index is=0; is<nspecies; is++ )
                { vmr_field(is,ip,ilat,ilon) = vmr_temp(is,ip,0,0); }
              for( Index is=0; is<t_nlte_field_raw.nelem(); is++ )
                t_nlte_field(is,ip,ilat,ilon) = nlte_temp(is,ip,0,0);
            }
        }
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MagFieldsCalcExpand1D(Tensor3&              mag_u_field,
                           Tensor3&              mag_v_field,
                           Tensor3&              mag_w_field,
                           const Vector&         p_grid,
                           const Vector&         lat_grid,
                           const Vector&         lon_grid,
                           const GriddedField3&  mag_u_field_raw,
                           const GriddedField3&  mag_v_field_raw,
                           const GriddedField3&  mag_w_field_raw,
                           const Index&          atmosphere_dim,
                           const Index&          interp_order,
                           const Verbosity&      verbosity)
{
    chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
    chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
    
    if( atmosphere_dim == 1 )
        throw std::runtime_error( "This function is intended for 2D and 3D. For 1D, use *MagFieldsCalc*.");
    
    // Make 1D interpolation using some dummy variables
    Vector    vempty(0);
    Tensor3   mag_u_field_temp, mag_v_field_temp, mag_w_field_temp;
    MagFieldsCalc( mag_u_field_temp, mag_v_field_temp, mag_w_field_temp,
                    p_grid, vempty, vempty, 
                    mag_u_field_raw, mag_v_field_raw, mag_w_field_raw,
                    /*atmosphere_dim = */ 1, 
                    interp_order, verbosity);
    
    // Move values from the temporary tensors to the return arguments
    const Index   np = p_grid.nelem();
    const Index   nlat = lat_grid.nelem();
    Index   nlon = lon_grid.nelem();
    if( atmosphere_dim == 2 )
    { nlon = 1; }
    //
    assert( mag_u_field_temp.npages() == np );
    assert( mag_v_field_temp.npages() == np );
    assert( mag_w_field_temp.npages() == np );
    //
    mag_u_field.resize( np, nlat, nlon );
    mag_v_field.resize( np, nlat, nlon );
    mag_w_field.resize( np, nlat, nlon );
    
    for( Index ilon=0; ilon<nlon; ilon++ )
    {
        for( Index ilat=0; ilat<nlat; ilat++ )
        {
            for( Index ip=0; ip<np; ip++ )
            {
                mag_u_field(ip,ilat,ilon) = mag_u_field_temp(ip,0,0);
                mag_v_field(ip,ilat,ilon) = mag_v_field_temp(ip,0,0);
                mag_w_field(ip,ilat,ilon) = mag_w_field_temp(ip,0,0);
            }
        }
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void WindFieldsCalcExpand1D(Tensor3&              wind_u_field,
                            Tensor3&              wind_v_field,
                            Tensor3&              wind_w_field,
                            const Vector&         p_grid,
                            const Vector&         lat_grid,
                            const Vector&         lon_grid,
                            const GriddedField3&  wind_u_field_raw,
                            const GriddedField3&  wind_v_field_raw,
                            const GriddedField3&  wind_w_field_raw,
                            const Index&          atmosphere_dim,
                            const Index&          interp_order,
                            const Verbosity&      verbosity)
{
    chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
    chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
    
    if( atmosphere_dim == 1 )
        throw std::runtime_error( "This function is intended for 2D and 3D. For 1D, use *WindFieldsCalc*.");
    
    // Make 1D interpolation using some dummy variables
    Vector    vempty(0);
    Tensor3   wind_u_field_temp, wind_v_field_temp, wind_w_field_temp;
    MagFieldsCalc( wind_u_field_temp, wind_v_field_temp, wind_w_field_temp,
                   p_grid, vempty, vempty, 
                   wind_u_field_raw, wind_v_field_raw, wind_w_field_raw,
                   /*atmosphere_dim = */ 1, 
                   interp_order, verbosity);
    
    // Move values from the temporary tensors to the return arguments
    const Index   np = p_grid.nelem();
    const Index   nlat = lat_grid.nelem();
    Index   nlon = lon_grid.nelem();
    if( atmosphere_dim == 2 )
    { nlon = 1; }
    //
    assert( wind_u_field_temp.npages() == np );
    assert( wind_v_field_temp.npages() == np );
    assert( wind_w_field_temp.npages() == np );
    //
    wind_u_field.resize( np, nlat, nlon );
    wind_v_field.resize( np, nlat, nlon );
    wind_w_field.resize( np, nlat, nlon );
    
    for( Index ilon=0; ilon<nlon; ilon++ )
    {
        for( Index ilat=0; ilat<nlat; ilat++ )
        {
            for( Index ip=0; ip<np; ip++ )
            {
                wind_u_field(ip,ilat,ilon) = wind_u_field_temp(ip,0,0);
                wind_v_field(ip,ilat,ilon) = wind_v_field_temp(ip,0,0);
                wind_w_field(ip,ilat,ilon) = wind_w_field_temp(ip,0,0);
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsExpand1D(Tensor3&         t_field,
                       Tensor3&         z_field,
                       Tensor4&         vmr_field,
                       const Vector&    p_grid,
                       const Vector&    lat_grid,
                       const Vector&    lon_grid,
                       const Index&     atmosphere_dim,
                       const Index&     chk_vmr_nan,
                       const Verbosity&)
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );

  // Sizes
  const Index   np = p_grid.nelem();
  const Index   nlat = lat_grid.nelem();
  const Index   nlon = max( Index(1), lon_grid.nelem() );
  const Index   nspecies = vmr_field.nbooks();

  const bool chknan = chk_vmr_nan;

  if( atmosphere_dim == 1 )
    { throw runtime_error( "No use in calling this method for 1D.");}
  chk_atm_field( "t_field", t_field, 1, p_grid, Vector(0), Vector(0) );
  chk_atm_field( "z_field", z_field, 1, p_grid, Vector(0), Vector(0) );
  if( nspecies )
    chk_atm_field( "vmr_field", vmr_field, 1, nspecies,
                   p_grid, Vector(0), Vector(0),
                   chknan );

  // Temporary containers
  Tensor3 t_temp = t_field, z_temp = z_field;
  Tensor4 vmr_temp = vmr_field;

  // Resize and fill
  t_field.resize( np, nlat, nlon );
  z_field.resize( np, nlat, nlon );
  vmr_field.resize( nspecies, np, nlat, nlon );
  //
  for( Index ilon=0; ilon<nlon; ilon++ )
    {
      for( Index ilat=0; ilat<nlat; ilat++ )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              t_field(ip,ilat,ilon) = t_temp(ip,0,0);
              z_field(ip,ilat,ilon) = z_temp(ip,0,0);
              for( Index is=0; is<nspecies; is++ )
                { vmr_field(is,ip,ilat,ilon) = vmr_temp(is,ip,0,0); }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsExtract1D(
                     Index&    atmosphere_dim,
                     Vector&   lat_grid,
                     Vector&   lon_grid,
                     Tensor3&  t_field,
                     Tensor3&  z_field,
                     Tensor4&  vmr_field,
              const Index&     ilat,
              const Index&     ilon,
              const Verbosity& verbosity )
{
  if( atmosphere_dim == 1 )
    { return; }

  if( ilat < 0  ||  ilat >= lat_grid.nelem() )
    throw runtime_error( "Invalid of *ilat*. It must be >= 0 and less than "
                         "length of *lat_grid*." );
  
  if( atmosphere_dim == 2 )
    {
      Vector vtmp;
      vtmp = t_field(joker,ilat,0);
      t_field.resize( t_field.npages(), 1, 1 );
      t_field(joker,0,0) = vtmp; 
      vtmp = z_field(joker,ilat,0);
      z_field.resize( z_field.npages(), 1, 1 );
      z_field(joker,0,0) = vtmp; 
      Matrix mtmp;
      mtmp = vmr_field(joker,joker,ilat,0);
      vmr_field.resize( vmr_field.nbooks(), vmr_field.npages(), 1, 1 );      
      vmr_field(joker,joker,0,0) = mtmp; 
    }
  else if( atmosphere_dim == 3 )
    {
      if( ilat < 0  ||  ilon >= lon_grid.nelem() )
        throw runtime_error( "Invalid of *ilon*. It must be >= 0 and less than "
                             "length of *lon_grid*." );
      Vector vtmp;
      vtmp = t_field(joker,ilat,ilon);
      t_field.resize( t_field.npages(), 1, 1 );
      t_field(joker,0,0) = vtmp; 
      vtmp = z_field(joker,ilat,ilon);
      z_field.resize( z_field.npages(), 1, 1 );
      z_field(joker,0,0) = vtmp; 
      Matrix mtmp;
      mtmp = vmr_field(joker,joker,ilat,ilon);
      vmr_field.resize( vmr_field.nbooks(), vmr_field.npages(), 1, 1 );      
      vmr_field(joker,joker,0,0) = mtmp; 
    }
  
  else 
    {
      throw runtime_error( "Invalid of *atmosphere_dim*. It must be 1-3." );
    }
  
  AtmosphereSet1D( atmosphere_dim, lat_grid, lon_grid, verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AtmFieldsRefinePgrid(// WS Output:
                          Vector& p_grid,
                          Tensor3& t_field,
                          Tensor3& z_field,
                          Tensor4& vmr_field,
                          Index& atmfields_checked,
                          Index& atmgeom_checked,
                          Index& cloudbox_checked,
                          // WS Input:
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Index& atmosphere_dim,
                          // Control Parameters:
                          const Numeric& p_step,
                          const Index& interp_order,
                          const Verbosity& verbosity)
{
  // Checks on input parameters:
  
  // We don't actually need lat_grid and lon_grid, but we have them as
  // input variables, so that we can use the standard functions to
  // check atmospheric fields and grids. A bit cheesy, but I don't
  // want to program all the checks explicitly.

  // Check grids:
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);

  // Check T field:
  chk_atm_field("t_field", t_field, atmosphere_dim,
                p_grid, lat_grid, lon_grid);
 
  // Check z field:
  chk_atm_field("z_field", z_field, atmosphere_dim,
                p_grid, lat_grid, lon_grid);
 
  // Check VMR field (and abs_species):
  chk_atm_field("vmr_field", vmr_field, atmosphere_dim,
                vmr_field.nbooks(), p_grid, lat_grid, lon_grid);

  // Move original p_grid to p_old, freeing p_grid for the refined one.
  Vector p_old;
  p_old = p_grid;

  p_gridRefine(p_grid,
               atmfields_checked, atmgeom_checked, cloudbox_checked,
               p_old, p_step, verbosity);

  AtmFieldPRegrid(z_field, z_field, p_grid, p_old, interp_order, verbosity);
  AtmFieldPRegrid(t_field, t_field, p_grid, p_old, interp_order, verbosity);
  AtmFieldPRegrid(vmr_field, vmr_field, p_grid, p_old, interp_order, verbosity);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void AtmRawRead(//WS Output:
                GriddedField3&        t_field_raw,
                GriddedField3&        z_field_raw,
                ArrayOfGriddedField3& vmr_field_raw,
                ArrayOfGriddedField3& t_nlte_field_raw,
                ArrayOfQuantumIdentifier& nlte_quantum_identifiers,
                //WS Input:
                const ArrayOfArrayOfSpeciesTag& abs_species,
                //Keyword:
                const String&    basename,
                const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  String tmp_basename = basename;
  if (basename.length() && basename[basename.length()-1] != '/')
    tmp_basename += ".";

      // Read the temperature field:
  String file_name = tmp_basename + "t.xml";
  xml_read_from_file( file_name, t_field_raw, verbosity);
  
  out3 << "Temperature field read from file: " << file_name << "\n";  

  // Read geometrical altitude field:
  file_name = tmp_basename + "z.xml";
  xml_read_from_file( file_name, z_field_raw, verbosity);

  out3 << "Altitude field read from file: " << file_name << "\n";  

  // Clear out vmr_field_raw
  vmr_field_raw.resize(0);

  // The species lookup data:
  using global_data::species_data;
  
  // We need to read one profile for each tag group.
  for ( Index i=0; i<abs_species.nelem(); i ++)
    {
      // Determine the name.
      file_name = tmp_basename +
        species_data[abs_species[i][0].Species()].Name() + ".xml";
      
      // Add an element for this tag group to the vmr profiles:
      GriddedField3 vmr_field_data;
      vmr_field_raw.push_back(vmr_field_data);
      
      // Read the VMR:
      xml_read_from_file( file_name, vmr_field_raw[vmr_field_raw.nelem()-1], verbosity);
      
      // state the source of profile.
      out3 << "  " << species_data[abs_species[i][0].Species()].Name()
           << " profile read from file: " << file_name << "\n";
    }
    
    // NLTE is ignored by doing this
    t_nlte_field_raw.resize(0);
    nlte_quantum_identifiers.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MagRawRead(//WS Output:
                GriddedField3&  mag_u_field_raw,
                GriddedField3&  mag_v_field_raw,
                GriddedField3&  mag_w_field_raw,
                //Keyword:
                const String&    basename,
                const Verbosity& verbosity)
{
    CREATE_OUT3;
    
    String tmp_basename = basename;
    if (basename.length() && basename[basename.length()-1] != '/')
        tmp_basename += ".";
    
    // Read magfield u component:
    String file_name = tmp_basename + "mag_u.xml";
    xml_read_from_file( file_name, mag_u_field_raw, verbosity);
    
    out3 << "Bu field read from file: " << file_name << "\n";  
    
    // Read magfield v component:
    file_name = tmp_basename + "mag_v.xml";
    xml_read_from_file( file_name, mag_v_field_raw, verbosity);
    
    out3 << "Bv field read from file: " << file_name << "\n"; 
    
    // Read magfield w component:
    file_name = tmp_basename + "mag_w.xml";
    xml_read_from_file( file_name, mag_w_field_raw, verbosity);
    
    out3 << "Bw field read from file: " << file_name << "\n"; 
}


/* Workspace method: Doxygen documentation will be auto-generated */
void WindRawRead(//WS Output:
                 GriddedField3&  wind_u_field_raw,
                 GriddedField3&  wind_v_field_raw,
                 GriddedField3&  wind_w_field_raw,
                 //Keyword:
                 const String&    basename,
                 const Verbosity& verbosity)
{
    CREATE_OUT3;
    
    String tmp_basename = basename;
    if (basename.length() && basename[basename.length()-1] != '/')
        tmp_basename += ".";
    
    // Read wind field u component:
    String file_name = tmp_basename + "wind_u.xml";
    xml_read_from_file( file_name, wind_u_field_raw, verbosity);
    
    out3 << "Wind u field read from file: " << file_name << "\n";  
    
    // Read wind field u component:
    file_name = tmp_basename + "wind_v.xml";
    xml_read_from_file( file_name, wind_v_field_raw, verbosity);
    
    out3 << "Wind v field read from file: " << file_name << "\n"; 
    
    // Read wind field u component:
    file_name = tmp_basename + "wind_w.xml";
    xml_read_from_file( file_name, wind_w_field_raw, verbosity);
    
    out3 << "Wind w field read from file: " << file_name << "\n"; 
}
  

/* Workspace method: Doxygen documentation will be auto-generated */
void AtmWithNLTERawRead(//WS Output:
                        GriddedField3&        t_field_raw,
                        GriddedField3&        z_field_raw,
                        ArrayOfGriddedField3& vmr_field_raw,
                        ArrayOfGriddedField3& t_nlte_field_raw,
                        ArrayOfQuantumIdentifier& nlte_quantum_identifiers,
                        //WS Input:
                        const ArrayOfArrayOfSpeciesTag& abs_species,
                        //Keyword:
                        const String&    basename,
                        const Verbosity& verbosity)
{
    CREATE_OUT3;
  
  String tmp_basename = basename;
  if (basename.length() && basename[basename.length()-1] != '/')
    tmp_basename += ".";

      // Read the temperature field:
  String file_name = tmp_basename + "t.xml";
  xml_read_from_file( file_name, t_field_raw, verbosity);
  
  out3 << "Temperature field read from file: " << file_name << "\n";  

  // Read geometrical altitude field:
  file_name = tmp_basename + "z.xml";
  xml_read_from_file( file_name, z_field_raw, verbosity);

  out3 << "Altitude field read from file: " << file_name << "\n";  

  // Clear out vmr_field_raw
  vmr_field_raw.resize(0);

  // The species lookup data:
  using global_data::species_data;
  
  // We need to read one profile for each tag group.
  for ( Index i=0; i<abs_species.nelem(); i ++)
    {
      // Determine the name.
      file_name = tmp_basename +
        species_data[abs_species[i][0].Species()].Name() + ".xml";
      
      // Add an element for this tag group to the vmr profiles:
      GriddedField3 vmr_field_data;
      vmr_field_raw.push_back(vmr_field_data);
      
      // Read the VMR:
      xml_read_from_file( file_name, vmr_field_raw[vmr_field_raw.nelem()-1], verbosity);
      
      // state the source of profile.
      out3 << "  " << species_data[abs_species[i][0].Species()].Name()
           << " profile read from file: " << file_name << "\n";
    }

  // Read each nlte temperature field:
  file_name = tmp_basename + "t_nlte.xml";
  xml_read_from_file( file_name, t_nlte_field_raw, verbosity);
  
  out3 << "NLTE temperature fieldarray read from file: " << file_name << "\n";  
  
  // Read each nlte identifier field:
  file_name = tmp_basename + "qi.xml";
  xml_read_from_file( file_name, nlte_quantum_identifiers, verbosity);
  
  out3 << "NLTE identifier array read from file: " << file_name << "\n";
  
  if(t_nlte_field_raw.nelem()!=nlte_quantum_identifiers.nelem())
  {
    ostringstream os;
    os << "The quantum identifers and the NLTE temperature fields\n"
       << "are of different lengths.  This should not be the case.\n"
       << "please check the qi.xml and t_nlte.xml files under\n"
       << basename << std::endl << "to correct this error.\n";
       throw std::runtime_error( os.str() );
  }
}
  

/* Workspace method: Doxygen documentation will be auto-generated */
void InterpAtmFieldToPosition(
          Numeric&   outvalue,
    const Index&     atmosphere_dim,
    const Vector&    p_grid,
    const Vector&    lat_grid,
    const Vector&    lon_grid,
    const Tensor3&   z_field,
    const Vector&    rtp_pos,
    const Tensor3&   field,
    const Verbosity& verbosity)
{
  // Input checks
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "input argument *field*", field, atmosphere_dim, 
                                                  p_grid, lat_grid, lon_grid );
  chk_rte_pos( atmosphere_dim, rtp_pos );

  // Determine grid positions
  GridPos gp_p, gp_lat, gp_lon;
  rte_pos2gridpos( gp_p, gp_lat, gp_lon, atmosphere_dim, 
                   p_grid, lat_grid, lon_grid, z_field, rtp_pos );

  // Interpolate
  outvalue = interp_atmfield_by_gp( atmosphere_dim, field, gp_p, gp_lat, 
                                                                 gp_lon );

  CREATE_OUT3;
  out3 << "    Result = " << outvalue << "\n";
}



/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridDensify(// WS Output:
                   Vector& p_grid,
                   Index& atmfields_checked,
                   Index& atmgeom_checked,
                   Index& cloudbox_checked,
                   // WS Input:
                   const Vector& p_grid_old,
                   // Control Parameters:
                   const Index&  nfill,
                   const Verbosity& verbosity )
{
  // Check that p_grid and p_grid_old are not the same variable (pointing to the
  // same memory space). this as p_grid will be overwritten, but we will need
  // both data later on for data regridding.
  if (&p_grid == &p_grid_old)
    {
      ostringstream os;
      os << "The old and new grids (p_grid and p_grid_old) are not allowed\n"
         << "to be identical (pointing to same memory space).\n"
         << "But they are doing in your case.";
      throw runtime_error( os.str() );
    }

  // as we manipoulate the overall vertical grid (but not simultaneously the
  // atmospheric fields), we reset all atmfields related checked WSV to
  // unchecked, forcing the user to do the checks again.
  atmfields_checked = 0;
  atmgeom_checked = 0;
  cloudbox_checked = 0;

  // Check the keyword argument:
  if( nfill < 0 ) 
    { throw runtime_error( "Argument *nfill* must be >= 0." ); }

  // Nothing to do if nfill=0
  if( nfill > 0 )
    {
      // Allocate new size for p_grid
      const Index n0 = p_grid_old.nelem();
      p_grid.resize( (n0-1)*(1+nfill) + 1 );

      Index iout = 0;
      p_grid[0] = p_grid_old[0];

      for( Index i=1; i<n0; i++ )
        {
          Vector pnew;
          VectorNLogSpace( pnew, 2+nfill, p_grid_old[i-1], p_grid_old[i], verbosity );
          for( Index j=1; j<nfill+2; j++ )
            {
              iout += 1;
              p_grid[iout] = pnew[j];
            }
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridRefine(// WS Output:
                  Vector& p_grid,
                  Index& atmfields_checked,
                  Index& atmgeom_checked,
                  Index& cloudbox_checked,
                  // WS Input:
                  const Vector& p_grid_old,
                  // Control Parameters:
                  const Numeric& p_step10,
                  const Verbosity&)
{
  // Check that p_grid and p_grid_old are not the same variable (pointing to the
  // same memory space). this as p_grid will be overwritten, but we will need
  // both data later on for data regridding.
  if (&p_grid == &p_grid_old)
    {
      ostringstream os;
      os << "The old and new grids (p_grid and p_grid_old) are not allowed\n"
         << "to be identical (pointing to same memory space).\n"
         << "But they are doing in your case.";
      throw runtime_error( os.str() );
    }

  // as we manipoulate the overall vertical grid (but not simultaneously the
  // atmospheric fields), we reset all atmfields related checked WSV to
  // unchecked, forcing the user to do the checks again.
  atmfields_checked = 0;
  atmgeom_checked = 0;
  cloudbox_checked = 0;

  // Check the keyword argument:
  if ( p_step10 <= 0  )
    {
      ostringstream os;
      os << "The keyword argument p_step must be >0.";
      throw runtime_error( os.str() );
    }

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

  for (Index i=1; i<log_p_old.nelem(); ++i)
    {
      const Numeric dp =  log_p_old[i-1] - log_p_old[i]; // The grid is descending.

      const Numeric dp_by_p_step = dp/p_step;
      //          cout << "dp_by_p_step: " << dp_by_p_step << "\n";

      // How many times does p_step fit into dp?
      const Index n = (Index) ceil(dp_by_p_step); 
      // n is the number of intervals that we want to have in the
      // new grid. The number of additional points to insert is
      // n-1. But we have to insert the original point as well.
      //          cout << n << "\n";

      const Numeric ddp = dp/(Numeric)n;
      //          cout << "ddp: " << ddp << "\n";

      for (Index j=1; j<=n; ++j)
        log_p_new.push_back(log_p_old[i-1] - (Numeric)j*ddp);          
    }

  // Copy ArrayOfNumeric to proper vector. 
  Vector log_p(log_p_new.nelem());
  for (Index i=0; i<log_p_new.nelem(); ++i)
    log_p[i] = log_p_new[i];

  // Copy the new grid to abs_p, removing the log:
  p_grid.resize(log_p.nelem());
  transform(p_grid, exp, log_p);
} 



/* Workspace method: Doxygen documentation will be auto-generated */
void p_gridFromZRaw(//WS Output
                      Vector& p_grid,
                      //WS Input
                      const GriddedField3& z_field_raw,
                      const Index& no_negZ,
                      const Verbosity&)
{
  // original version excludes negative z. not clear, why this is. maybe is
  // currently a convention somehwere in ARTS (DOIT?). negative z seem, however,
  // to work fine for clear-sky cases. so we make the negative z exclude an
  // option (for consistency until unclarities solved, default: do exclude)
  Vector p_grid_raw=z_field_raw.get_numeric_grid(GFIELD3_P_GRID);

  Index i;
  if (is_increasing(z_field_raw.data(joker,0,0)))
    {
      i=0;
      if (no_negZ)
        {
          while ( z_field_raw.data(i,0,0)< 0.0 ) i++;
        }
      p_grid=p_grid_raw[Range(i,joker)];
    }
  else if (is_decreasing(z_field_raw.data(joker,0,0)))
    {
      i=z_field_raw.data.npages()-1;
      if (no_negZ)
        {
          while ( z_field_raw.data(i,0,0)< 0.0 ) i--;
        }
      p_grid=p_grid_raw[Range(i,joker,-1)];
    }
  else
    {
       ostringstream os;
       os << "z_field_raw needs to be monotonous, but this is not the case.\n";
       throw runtime_error( os.str() );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void lat_gridFromRawField(//WS Output
                      Vector& lat_grid,
                      //WS Input
                      const GriddedField3& field_raw,
                      const Verbosity&)
{
  lat_grid = field_raw.get_numeric_grid(GFIELD3_LAT_GRID);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void lon_gridFromRawField(//WS Output
                      Vector& lon_grid,
                      //WS Input
                      const GriddedField3& field_raw,
                      const Verbosity&)
{
  lon_grid = field_raw.get_numeric_grid(GFIELD3_LON_GRID);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void wind_u_fieldIncludePlanetRotation(
         Tensor3&   wind_u_field,
   const Index&     atmosphere_dim,
   const Vector&    p_grid,
   const Vector&    lat_grid,
   const Vector&    lon_grid,
   const Vector&    refellipsoid,
   const Tensor3&   z_field,
   const Numeric&   planet_rotation_period,
   const Verbosity&)
{
  if( atmosphere_dim < 3 )
    throw runtime_error( "No need to use this method for 1D and 2D." );

  const Index  np = p_grid.nelem();
  const Index  na = lat_grid.nelem();
  const Index  no = lon_grid.nelem();
  
  chk_atm_field( "z_field", z_field, atmosphere_dim, 
                                                  p_grid, lat_grid, lon_grid );
  if( wind_u_field.npages() > 0 ) 
    { chk_atm_field( "wind_u_field", wind_u_field, atmosphere_dim, 
                                                  p_grid, lat_grid, lon_grid );}
  else
    {
        wind_u_field.resize( np, na, no );
        wind_u_field = 0.;
    }

  const Numeric k1 = 2 * PI / planet_rotation_period;

  
  for( Index a=0; a<na; a++ )
    {
      const Numeric k2 = k1 * cos( DEG2RAD*lat_grid[a] );
      const Numeric re = refell2r( refellipsoid, lat_grid[a] );

      for( Index o=0; o<no; o++ )
        {
          for( Index p=0; p<np; p++ )
            {
              wind_u_field(p,a,o) += k2 * ( re + z_field(p,a,o) );
            }
        }
    }
}




// A small help function
void z2g(
               Numeric& g,
         const Numeric& r,
         const Numeric& g0,
         const Numeric& z )
{
  const Numeric x = r/(r+z);
  g = g0 * x*x;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void z_fieldFromHSE(
         Workspace&   ws,
         Tensor3&     z_field,
   const Index&       atmosphere_dim,
   const Vector&      p_grid,
   const Vector&      lat_grid,
   const Vector&      lon_grid,
   const Vector&      lat_true,
   const Vector&      lon_true,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Tensor3&     t_field,
   const Tensor4&     vmr_field,
   const Vector&      refellipsoid,
   const Matrix&      z_surface,
   const Index&       atmfields_checked,
   const Agenda&      g0_agenda,
   const Numeric&     molarmass_dry_air,
   const Numeric&     p_hse,
   const Numeric&     z_hse_accuracy,
   const Verbosity&   verbosity)
{
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );

  // Some general variables
  const Index np   = p_grid.nelem();
  const Index nlat = t_field.nrows();
  const Index nlon = t_field.ncols();
  //
  const Index firstH2O = find_first_species_tg( abs_species,
                                      species_index_from_species_name("H2O") );

  if( firstH2O < 0 )
    {
      ostringstream os;
      os << "No water vapour tag group in *abs_species*.\n"
         << "Be aware that this leads to significant variations in atmospheres\n"
         << "that contain considerable amounts of water vapour (e.g. Earth)!\n";
      CREATE_OUT1;
      out1 << os.str();
      //throw runtime_error(
      //   "Water vapour is required (must be a tag group in *abs_species*)." );
    }
  //
  if( p_hse>p_grid[0]  ||  p_hse < p_grid[np-1] )
    {
      ostringstream os;
      os << "The value of *p_hse* must be inside the range of *p_grid*:"
         << "  p_hse  = " << p_hse << " Pa\n"
         << "  p_grid = << p_grid[np-1]" << " - " << p_grid[0] << " Pa\n";
      throw runtime_error( os.str() );
    }
  //
  if( z_hse_accuracy <= 0 )
    { throw runtime_error( "The value of *z_hse_accuracy* must be > 0." ); }
  //
  chk_latlon_true( atmosphere_dim, lat_grid, lat_true, lon_true );


  // Determine interpolation weights for p_hse
  //
  ArrayOfGridPos gp(1);
  Matrix itw( 1, 2);
  p2gridpos( gp, p_grid, Vector(1,p_hse) );
  interpweights ( itw, gp );


  // // Molecular weight of water vapour
  const Numeric mw = 18.016;   

  // mw/molarmass_dry_air matches eps in Eq. 3.14 in Wallace&Hobbs:
  const Numeric k  = 1 - mw/molarmass_dry_air; 

  // Gas constant for 1kg dry air:
  const Numeric rd = 1e3*GAS_CONSTANT/molarmass_dry_air; 

  // The calculations
  //
  for( Index ilat=0; ilat<nlat; ilat++ )
    {
      // The reference ellipsoid is already adjusted to internal 1D or 2D
      // views, and lat_grid is the relevant grid for *refellipsoid*, also
      // for 2D. On the other hand, extraction of g0 requires that the true
      // position is determined.

      // Radius of reference ellipsoid
      Numeric re;
      if( atmosphere_dim == 1 )
        { re = refellipsoid[0]; }
      else
        { re = refell2r( refellipsoid, lat_grid[ilat] ); }

      for( Index ilon=0; ilon<nlon; ilon++ )
        {
          // Determine true latitude and longitude
          Numeric lat, lon;
          Vector  pos(atmosphere_dim);    // pos[0] can be a dummy value
          if( atmosphere_dim >= 2 )
            {
              pos[1] = lat_grid[ilat];
              if( atmosphere_dim == 3 )
                { pos[2] = lon_grid[ilon]; }
            }
          pos2true_latlon( lat, lon, atmosphere_dim, lat_grid, lat_true, 
                                                               lon_true, pos );

          // Get g0
          Numeric g0;
          g0_agendaExecute( ws, g0, lat, lon, g0_agenda );

          // Determine altitude for p_hse
          Vector z_hse(1);
          interp( z_hse, itw, z_field(joker,ilat,ilon), gp );

          Numeric z_acc = 2 * z_hse_accuracy;

          while( z_acc > z_hse_accuracy )
            {
              z_acc = 0;
              Numeric g2=g0;

              for( Index ip=0; ip<np-1; ip++ )
                {
                  // Calculate average g
                  if( ip == 0 )
                    { z2g( g2, re, g0, z_field(ip,ilat,ilon) ); }
                  const Numeric g1 = g2;
                  z2g( g2, re, g0, z_field(ip+1,ilat,ilon) );
                  //
                  const Numeric g = ( g1 + g2 ) / 2.0;

                  //Average water vapour VMR
                  Numeric hm;
                  if( firstH2O < 0 )
                    {
                      hm = 0.0;
                    }
                  else
                    {
                      hm = 0.5* ( vmr_field(firstH2O,ip,ilat,ilon)
                                + vmr_field(firstH2O,ip+1,ilat,ilon) );
                    }
  
                  // Average virtual temperature (no liquid water)
                  // (cf. e.g. Eq. 3.16 in Wallace&Hobbs)
                  const Numeric tv = ( 1 / ( 2 * (1-hm*k) ) ) * 
                           ( t_field(ip,ilat,ilon) + t_field(ip+1,ilat,ilon) );
  
                  // Change in vertical altitude from ip to ip+1 
                  //  (cf. e.g. Eq. 3.24 in Wallace&Hobbs)
                  const Numeric dz = rd*(tv/g)*log( p_grid[ip]/p_grid[ip+1] );

                  // New altitude
                  Numeric znew = z_field(ip,ilat,ilon) + dz;
                  z_acc = max( z_acc, fabs(znew-z_field(ip+1,ilat,ilon)) );
                  z_field(ip+1,ilat,ilon) = znew;
                }

              // Adjust to z_hse
              Vector z_tmp(1);
              interp( z_tmp, itw, z_field(joker,ilat,ilon), gp );
              z_field(joker,ilat,ilon) -= z_tmp[0] - z_hse[0];
            }
        }
    }

  // Check that there is no gap between the surface and lowest pressure 
  // level
  // (This code is copied from *basics_checkedCalc*. Make this to an internal
  // function if used in more places.)
  for( Index row=0; row<z_surface.nrows(); row++ )
    {
      for( Index col=0; col<z_surface.ncols(); col++ )
        {
          if( z_surface(row,col)<z_field(0,row,col) ||
                   z_surface(row,col)>=z_field(z_field.npages()-1,row,col) )
            {
              ostringstream os;
              os << "The surface altitude (*z_surface*) cannot be outside "
                 << "of the altitudes in *z_field*.";
              throw runtime_error( os.str() );
            }
        }
    }

}




/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldSetConstant(
   Tensor4&         vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const String&    species,
   const Numeric&   vmr_value, 
   const Verbosity&  )
{
  // Check input
  chk_if_in_range( "vmr_value", vmr_value, 0, 1 ); 
  //
  if( abs_species.nelem() != vmr_field.nbooks() )
    throw runtime_error( 
          "Size of *vmr_field* and length of *abs_species* do not agree." );

  // Find position for this species. 
  ArrayOfSpeciesTag tag;
  array_species_tag_from_string( tag, species );
  Index si = chk_contains( "species", abs_species, tag );

  // Apply value
  vmr_field(si,joker,joker,joker) = vmr_value;
}


