/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

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
void nca_read_from_file(const int ncid, ArrayOfMatrix& aom, const Verbosity&) {
  Index nelem;
  nelem = nc_get_dim(ncid, "nelem");

  long* vnrows = new long[nelem];
  long* vncols = new long[nelem];
  aom.resize(nelem);
  nca_get_data(ncid, "Matrix_nrows", vnrows);
  nca_get_data(ncid, "Matrix_ncols", vncols);
  size_t pos = 0;
  for (Index i = 0; i < nelem; i++) {
    aom[i].resize(vnrows[i], vncols[i]);
    nca_get_data(ncid,
                         "ArrayOfMatrix",
                         pos,
                         vnrows[i] * vncols[i],
                         aom[i].get_c_array());
    pos += vnrows[i] * vncols[i];
  }

  delete[] vnrows;
  delete[] vncols;
}

//! Writes an ArrayOfMatrix to a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param aom     ArrayOfMatrix
*/
void nca_write_to_file(const int ncid,
                       const ArrayOfMatrix& aom,
                       const Verbosity&) {
  int retval;
  int ncdim, varid_nrows, varid_ncols;
  int ncdim_total, varid;
  long nelem_total = 0;
  long* vncols = new long[aom.nelem()];
  long* vnrows = new long[aom.nelem()];
  for (Index i = 0; i < aom.nelem(); i++) {
    vnrows[i] = aom[i].nrows();
    vncols[i] = aom[i].ncols();
    nelem_total += vnrows[i] * vncols[i];
  }

  if ((retval = nc_def_dim(ncid, "nelem", aom.nelem(), &ncdim)))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "nelem_total", nelem_total, &ncdim_total)))
    nca_error(retval, "nc_def_dim");

  if ((retval =
           nc_def_var(ncid, "Matrix_nrows", NC_INT64, 1, &ncdim, &varid_nrows)))
    nca_error(retval, "nc_def_var");
  if ((retval =
           nc_def_var(ncid, "Matrix_ncols", NC_INT64, 1, &ncdim, &varid_ncols)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_def_var(
           ncid, "ArrayOfMatrix", NC_DOUBLE, 1, &ncdim_total, &varid)))
    nca_error(retval, "nc_def_var");

  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");

  if ((retval = nc_put_var_long(ncid, varid_nrows, vnrows)))
    nca_error(retval, "nc_put_var");
  if ((retval = nc_put_var_long(ncid, varid_ncols, vncols)))
    nca_error(retval, "nc_put_var");

  size_t pos = 0;
  for (Index i = 0; i < aom.nelem(); i++) {
    size_t count = aom[i].nrows() * aom[i].ncols();
    if ((retval = nc_put_vara_double(
             ncid, varid, &pos, &count, aom[i].get_c_array())))
      nca_error(retval, "nc_put_var");
    pos += count;
  }

  delete[] vnrows;
  delete[] vncols;
}

//=== ArrayOfVector ==========================================================

//! Reads an ArrayOfVector from a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param aov     ArrayOfVector
*/
void nca_read_from_file(const int ncid, ArrayOfVector& aov, const Verbosity&) {
  Index nelem;
  nelem = nc_get_dim(ncid, "nelem");

  long* vnelem = new long[nelem];
  aov.resize(nelem);
  nca_get_data(ncid, "Vector_nelem", vnelem);
  size_t pos = 0;
  for (Index i = 0; i < nelem; i++) {
    aov[i].resize(vnelem[i]);
    nca_get_data(
        ncid, "ArrayOfVector", pos, vnelem[i], aov[i].get_c_array());
    pos += vnelem[i];
  }

  delete[] vnelem;
}

//! Writes an ArrayOfVector to a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param aov     ArrayOfVector
*/
void nca_write_to_file(const int ncid,
                       const ArrayOfVector& aov,
                       const Verbosity&) {
  int retval;
  int ncdim, varid_nelem;
  int ncdim_total, varid;
  long nelem_total = 0;
  long* velems = new long[aov.nelem()];
  for (Index i = 0; i < aov.nelem(); i++) {
    velems[i] = aov[i].nelem();
    nelem_total += velems[i];
  }

  if ((retval = nc_def_dim(ncid, "nelem", aov.nelem(), &ncdim)))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "nelem_total", nelem_total, &ncdim_total)))
    nca_error(retval, "nc_def_dim");

  if ((retval =
           nc_def_var(ncid, "Vector_nelem", NC_INT64, 1, &ncdim, &varid_nelem)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_def_var(
           ncid, "ArrayOfVector", NC_DOUBLE, 1, &ncdim_total, &varid)))
    nca_error(retval, "nc_def_var");

  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");

  if ((retval = nc_put_var_long(ncid, varid_nelem, velems)))
    nca_error(retval, "nc_put_var");

  size_t pos = 0;
  for (Index i = 0; i < aov.nelem(); i++) {
    size_t count = aov[i].nelem();
    if ((retval = nc_put_vara_double(
             ncid, varid, &pos, &count, aov[i].get_c_array())))
      nca_error(retval, "nc_put_var");
    pos += count;
  }

  delete[] velems;
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what)                                   \
  void nca_write_to_file(const int, const what&, const Verbosity&) {          \
    throw runtime_error("NetCDF support not yet implemented for this type!"); \
  }                                                                           \
  void nca_read_from_file(const int, what&, const Verbosity&) {               \
    throw runtime_error("NetCDF support not yet implemented for this type!"); \
  }

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY

#endif /* ENABLE_NETCDF */
