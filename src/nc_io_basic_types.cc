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

  \brief This file contains basic functions to handle NetCDF data files.

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
nc_read_from_file (NcFile &ncf _U_,
                   Matrix& m _U_)
{
  if (ncf.num_dims() != 2 || ncf.num_vars() != 1)
    {
      ostringstream os;
      os << "Error reading NetCDF file." << endl;
      throw runtime_error (os.str());
    }

  NcVar *data = ncf.get_var("Matrix");
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
  NcDim *cols = ncf.add_dim ("cols", m.ncols());
  NcDim *rows = ncf.add_dim ("rows", m.nrows());

  NcVar *data = ncf.add_var ("Matrix", ncDouble, rows, cols);

  const Numeric *np = m.get_c_array();
  data->put (np, m.nrows(), m.ncols());

}

//=== Vector ==========================================================

//! Reads a Vector from a NetCDF file
/*!
  \param ncf     NetCDF file discriptor
  \param v       Vector
*/
void
nc_read_from_file (NcFile &ncf,
                   Vector& v _U_)
{
  if (ncf.num_dims() != 1 || ncf.num_vars() != 1)
    {
      ostringstream os;
      os << "Error reading NetCDF file." << endl;
      throw runtime_error (os.str());
    }

  NcVar *data = ncf.get_var("Vector");
  v.resize (data->values()->num());
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
  NcDim *dim = ncf.add_dim ("nelem", v.nelem());

  NcVar *data = ncf.add_var ("Vector", ncDouble, dim);

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
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor5 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Tensor7 )
TMPL_NC_READ_WRITE_FILE_DUMMY( Timer )

//=== Compound Types =======================================================

TMPL_NC_READ_WRITE_FILE_DUMMY( Agenda )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GasAbsLookup )
TMPL_NC_READ_WRITE_FILE_DUMMY( GridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( IsotopeRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( MCAntenna )
TMPL_NC_READ_WRITE_FILE_DUMMY( Ppath )
TMPL_NC_READ_WRITE_FILE_DUMMY( RetrievalQuantity )
TMPL_NC_READ_WRITE_FILE_DUMMY( SLIData2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( SingleScatteringData )
TMPL_NC_READ_WRITE_FILE_DUMMY( SpeciesRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( SpeciesTag )

//=== Array Types ==========================================================

TMPL_NC_READ_WRITE_FILE_DUMMY( Array<IsotopeRecord> )
TMPL_NC_READ_WRITE_FILE_DUMMY( Array<SpeciesRecord> )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfIndex )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfLineRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfMatrix )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfSpeciesTag )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfTensor3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfTensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfIndex )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfLineRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfLineshapeSpec )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfMatrix )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfPpath )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfRetrievalQuantity )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfSingleScatteringData )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfSpeciesTag )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfString )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor7 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfVector )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY


