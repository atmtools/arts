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
  \file   nc_io_array_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-10-01

  \brief This file contains functions to handle NetCDF data files.

*/

#include "arts.h"

#ifdef ENABLE_NETCDF

#include "nc_io.h"
#include "nc_io_types.h"

//=== ArrayOfMatrix ==========================================================

//! Reads an ArrayOfMatrix from a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param aom     ArrayOfMatrix
*/
void
nc_read_from_file (const int ncid,
                   ArrayOfMatrix& aom)
{
  Index nelem, nelem_total;
  nelem = nc_get_dim (ncid, "nelem");
  nelem_total = nc_get_dim (ncid, "nelem_total");

  long *vnrows = new long[nelem];
  long *vncols = new long[nelem];
  aom.resize (nelem);
  nc_get_data_long (ncid, "Matrix_nrows", vnrows);
  nc_get_data_long (ncid, "Matrix_ncols", vncols);
  size_t pos = 0;
  for (Index i = 0; i < nelem; i++)
    {
      aom[i].resize (vnrows[i], vncols[i]);
      nc_get_dataa_double (ncid, "ArrayOfMatrix", pos, vnrows[i] * vncols[i],
                           aom[i].get_c_array());
      pos += vnrows[i] * vncols[i];
    }

  delete [] vnrows;
  delete [] vncols;
}


//! Writes an ArrayOfMatrix to a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param aom     ArrayOfMatrix
*/
void
nc_write_to_file (const int ncid,
                  const ArrayOfMatrix& aom)
{
  int retval;
  int ncdim, varid_nrows, varid_ncols;
  int ncdim_total, varid;
  long nelem_total = 0;
  long *vncols = new long[aom.nelem()];
  long *vnrows = new long[aom.nelem()];
  for (Index i = 0; i < aom.nelem(); i++)
    {
      vnrows[i] = aom[i].nrows();
      vncols[i] = aom[i].ncols();
      nelem_total += vnrows[i] * vncols[i];
    }

  if ((retval = nc_def_dim (ncid, "nelem", aom.nelem(), &ncdim)))
    ncerror (retval, "nc_def_dim");
  if ((retval = nc_def_dim (ncid, "nelem_total", nelem_total, &ncdim_total)))
    ncerror (retval, "nc_def_dim");

  if ((retval = nc_def_var (ncid, "Matrix_nrows", NC_LONG, 1,
                            &ncdim, &varid_nrows)))
    ncerror (retval, "nc_def_var");
  if ((retval = nc_def_var (ncid, "Matrix_ncols", NC_LONG, 1,
                            &ncdim, &varid_ncols)))
    ncerror (retval, "nc_def_var");
  if ((retval = nc_def_var (ncid, "ArrayOfMatrix", NC_DOUBLE, 1,
                            &ncdim_total, &varid)))
    ncerror (retval, "nc_def_var");

  if ((retval = nc_enddef (ncid))) ncerror (retval, "nc_enddef");

  if ((retval = nc_put_var_long (ncid, varid_nrows, vnrows)))
    ncerror (retval, "nc_put_var");
  if ((retval = nc_put_var_long (ncid, varid_ncols, vncols)))
    ncerror (retval, "nc_put_var");

  size_t pos = 0;
  for (Index i = 0; i < aom.nelem(); i++)
    {
      size_t count = aom[i].nrows() * aom[i].ncols();
      if ((retval = nc_put_vara_double (ncid, varid, &pos, &count,
                                        aom[i].get_c_array())))
        ncerror (retval, "nc_put_var");
      pos += count;
    }

  delete [] vnrows;
  delete [] vncols;
}


//=== ArrayOfVector ==========================================================

//! Reads an ArrayOfVector from a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param aov     ArrayOfVector
*/
void
nc_read_from_file (const int ncid,
                   ArrayOfVector& aov)
{
  Index nelem, nelem_total;
  nelem = nc_get_dim (ncid, "nelem");
  nelem_total = nc_get_dim (ncid, "nelem_total");

  long *vnelem = new long[nelem];
  aov.resize (nelem);
  nc_get_data_long (ncid, "Vector_nelem", vnelem);
  size_t pos = 0;
  for (Index i = 0; i < nelem; i++)
    {
      aov[i].resize (vnelem[i]);
      nc_get_dataa_double (ncid, "ArrayOfVector", pos, vnelem[i],
                           aov[i].get_c_array());
      pos += vnelem[i];
    }

  delete [] vnelem;
}


//! Writes an ArrayOfVector to a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param aov     ArrayOfVector
*/
void
nc_write_to_file (const int ncid,
                  const ArrayOfVector& aov)
{
  int retval;
  int ncdim, varid_nelem;
  int ncdim_total, varid;
  long nelem_total = 0;
  long *velems = new long[aov.nelem()];
  for (Index i = 0; i < aov.nelem(); i++)
    {
      velems[i] = aov[i].nelem();
      nelem_total += velems[i];
    }

  if ((retval = nc_def_dim (ncid, "nelem", aov.nelem(), &ncdim)))
    ncerror (retval, "nc_def_dim");
  if ((retval = nc_def_dim (ncid, "nelem_total", nelem_total, &ncdim_total)))
    ncerror (retval, "nc_def_dim");

  if ((retval = nc_def_var (ncid, "Vector_nelem", NC_LONG, 1,
                            &ncdim, &varid_nelem)))
    ncerror (retval, "nc_def_var");
  if ((retval = nc_def_var (ncid, "ArrayOfVector", NC_DOUBLE, 1,
                            &ncdim_total, &varid)))
    ncerror (retval, "nc_def_var");

  if ((retval = nc_enddef (ncid))) ncerror (retval, "nc_enddef");

  if ((retval = nc_put_var_long (ncid, varid_nelem, velems)))
    ncerror (retval, "nc_put_var");

  size_t pos = 0;
  for (Index i = 0; i < aov.nelem(); i++)
    {
      size_t count = aov[i].nelem();
      if ((retval = nc_put_vara_double (ncid, varid, &pos, &count,
                                        aov[i].get_c_array())))
        ncerror (retval, "nc_put_var");
      pos += count;
    }

  delete [] velems;
}


////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what) \
  void nc_write_to_file (const int, const what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  } \
  void nc_read_from_file (const int, what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  }

//=== Array Types ==========================================================

TMPL_NC_READ_WRITE_FILE_DUMMY( Array<IsotopeRecord> )
TMPL_NC_READ_WRITE_FILE_DUMMY( Array<SpeciesRecord> )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGriddedField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGriddedField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfIndex )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfLineRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfMatrix )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfSpeciesTag )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfTensor3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfTensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGriddedField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGriddedField2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGriddedField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGriddedField4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfIndex )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfLineRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfLineshapeSpec )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfPpath )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfRetrievalQuantity )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfSingleScatteringData )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfSpeciesTag )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfString )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor7 )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY

#endif /* ENABLE_NETCDF */


