/* Copyright (C) 2003-2008 Oliver Lemke <olemke@core-dump.info>

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


////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   nc_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-26

  \brief This file contains functions to handle NetCDF data files.

*/

#include "arts.h"
#include "nc_io.h"
#include "nc_io_types.h"


//=== Matrix ==========================================================

//! Reads a Matrix from a NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param m       Matrix
*/
void
nc_read_from_file (NcFile &ncf,
                   Matrix& m)
{
  NcVar *data;
  nc_read_var (ncf, &data, 2, "Matrix");

  NcDim *rows = data->get_dim(0);
  NcDim *cols = data->get_dim(1);
  m.resize (rows->size(), cols->size());

  Index i = 0;
  for (Index r = 0; r < m.nrows(); r++)
    for (Index c = 0; c < m.ncols(); c++, i++)
      m(r, c) = data->values()->as_double(i);
}


//! Writes a Matrix to a NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param m       Matrix
*/
void
nc_write_to_file (NcFile &ncf,
                  const Matrix& m)
{
  NcDim *nrows = ncf.add_dim ("nrows", m.nrows());
  NcDim *ncols = ncf.add_dim ("ncols", m.ncols());

  NcVar *data = ncf.add_var ("Matrix", ncDouble, nrows, ncols);

  const Numeric *np = m.get_c_array();
  data->put (np, m.nrows(), m.ncols());

}

//=== Tensor3 ==========================================================

//! Reads a Tensor3 from a NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param t       Tensor3
*/
void
nc_read_from_file (NcFile &ncf,
                   Tensor3& t)
{
  NcVar *data;
  nc_read_var (ncf, &data, 3, "Tensor3");

  NcDim *npages = data->get_dim(0);
  NcDim *nrows = data->get_dim(1);
  NcDim *ncols = data->get_dim(2);
  t.resize (npages->size(), nrows->size(), ncols->size());

  Index i = 0;
  for (Index p = 0; p < t.npages(); p++)
    for (Index r = 0; r < t.nrows(); r++)
      for (Index c = 0; c < t.ncols(); c++, i++)
        t(p, r, c) = data->values()->as_double(i);
}


//! Writes a Tensor3 to a NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param t       Tensor3
*/
void
nc_write_to_file (NcFile &ncf,
                  const Tensor3& t)
{
  NcDim *npages = ncf.add_dim ("npages", t.npages());
  NcDim *nrows  = ncf.add_dim ("nrows", t.nrows());
  NcDim *ncols  = ncf.add_dim ("ncols", t.ncols());

  NcVar *data = ncf.add_var ("Tensor3", ncDouble, npages, nrows, ncols);

  const Numeric *np = t.get_c_array();
  data->put (np, t.npages(), t.nrows(), t.ncols());
}

//=== Vector ==========================================================

//! Reads a Vector from a NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param v       Vector
*/
void
nc_read_from_file (NcFile &ncf,
                   Vector& v)
{
  NcVar *data;
  nc_read_var (ncf, &data, 1, "Vector");

  NcDim *nelem = data->get_dim(0);
  v.resize (nelem->size());

  for (Index i = 0; i < v.nelem(); i++)
    v[i] = data->values()->as_double(i);
}


//! Writes a Vector to a NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param v       Vector
*/
void
nc_write_to_file (NcFile &ncf,
                  const Vector& v)
{
  NcDim *nelem = ncf.add_dim ("nelem", v.nelem());

  NcVar *data = ncf.add_var ("Vector", ncDouble, nelem);

  const Numeric *np = v.get_c_array();
  data->put (np, v.nelem());
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what) \
  void nc_write_to_file (NcFile&, const what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  } \
  void nc_read_from_file (NcFile&, what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  }

//=== Basic Types ==========================================================

TMPL_NC_READ_WRITE_FILE_DUMMY( Index )
TMPL_NC_READ_WRITE_FILE_DUMMY( Numeric )
TMPL_NC_READ_WRITE_FILE_DUMMY( Sparse )
TMPL_NC_READ_WRITE_FILE_DUMMY( String )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor5 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor7 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Timer )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY


