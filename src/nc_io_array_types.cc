////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   nc_io_array_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-10-01

  \brief This file contains functions to handle NetCDF data files.

*/

#include "nc_io.h"
#include "nc_io_types.h"

//=== ArrayOfIndex ==========================================================

//! Reads a ArrayOfIndex from a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param v       ArrayOfIndex
*/
void nca_read_from_file(const int ncid, ArrayOfIndex& v) {
  Index size;
  size = nca_get_dim(ncid, "size");

  v.resize(size);
  nca_get_data(ncid, "ArrayOfIndex", v.data());
}

//! Writes a ArrayOfIndex to a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param v       ArrayOfIndex
*/
void nca_write_to_file(const int ncid, const ArrayOfIndex& v) {
  int retval;
  int ncdim, varid;
  if ((retval = nc_def_dim(ncid, "size", v.size(), &ncdim)))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_var(ncid, "ArrayOfIndex", NC_INT64, 1, &ncdim, &varid)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");
  if ((retval = nc_put_var(ncid, varid, v.data())))
    nca_error(retval, "nc_put_var");
}

//=== ArrayOfMatrix ==========================================================

//! Reads an ArrayOfMatrix from a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param aom     ArrayOfMatrix
*/
void nca_read_from_file(const int ncid, ArrayOfMatrix& aom) {
  Index size;
  size = nca_get_dim(ncid, "size");

  long* vnrows = new long[size];
  long* vncols = new long[size];
  aom.resize(size);
  nca_get_data(ncid, "Matrix_nrows", vnrows);
  nca_get_data(ncid, "Matrix_ncols", vncols);
  size_t pos = 0;
  for (Index i = 0; i < size; i++) {
    aom[i].resize(vnrows[i], vncols[i]);
    nca_get_data(ncid,
                         "ArrayOfMatrix",
                         pos,
                         vnrows[i] * vncols[i],
                         aom[i].unsafe_data_handle());
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
                       const ArrayOfMatrix& aom) {
  int retval;
  int ncdim, varid_nrows, varid_ncols;
  int ncdim_total, varid;
  long size_total = 0;
  long* vncols = new long[aom.size()];
  long* vnrows = new long[aom.size()];
  for (Index i = 0; i < aom.size(); i++) {
    vnrows[i] = aom[i].nrows();
    vncols[i] = aom[i].ncols();
    size_total += vnrows[i] * vncols[i];
  }

  if ((retval = nc_def_dim(ncid, "size", aom.size(), &ncdim)))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "size_total", size_total, &ncdim_total)))
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
  for (Index i = 0; i < aom.size(); i++) {
    size_t count = aom[i].nrows() * aom[i].ncols();
    if ((retval = nc_put_vara_double(
             ncid, varid, &pos, &count, aom[i].unsafe_data_handle())))
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
void nca_read_from_file(const int ncid, ArrayOfVector& aov) {
  Index size;
  size = nca_get_dim(ncid, "size");

  long* vsize = new long[size];
  aov.resize(size);
  nca_get_data(ncid, "Vector_size", vsize);
  size_t pos = 0;
  for (Index i = 0; i < size; i++) {
    aov[i].resize(vsize[i]);
    nca_get_data(
        ncid, "ArrayOfVector", pos, vsize[i], aov[i].unsafe_data_handle());
    pos += vsize[i];
  }

  delete[] vsize;
}

//! Writes an ArrayOfVector to a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param aov     ArrayOfVector
*/
void nca_write_to_file(const int ncid,
                       const ArrayOfVector& aov) {
  int retval;
  int ncdim, varid_size;
  int ncdim_total, varid;
  long size_total = 0;
  long* velems = new long[aov.size()];
  for (Index i = 0; i < aov.size(); i++) {
    velems[i] = aov[i].size();
    size_total += velems[i];
  }

  if ((retval = nc_def_dim(ncid, "size", aov.size(), &ncdim)))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "size_total", size_total, &ncdim_total)))
    nca_error(retval, "nc_def_dim");

  if ((retval =
           nc_def_var(ncid, "Vector_size", NC_INT64, 1, &ncdim, &varid_size)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_def_var(
           ncid, "ArrayOfVector", NC_DOUBLE, 1, &ncdim_total, &varid)))
    nca_error(retval, "nc_def_var");

  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");

  if ((retval = nc_put_var_long(ncid, varid_size, velems)))
    nca_error(retval, "nc_put_var");

  size_t pos = 0;
  for (Index i = 0; i < aov.size(); i++) {
    size_t count = aov[i].size();
    if ((retval = nc_put_vara_double(
             ncid, varid, &pos, &count, aov[i].unsafe_data_handle())))
      nca_error(retval, "nc_put_var");
    pos += count;
  }

  delete[] velems;
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what)                                    \
  void nca_write_to_file(const int, const what &) {         \
    ARTS_USER_ERROR("NetCDF support not yet implemented for this type!");      \
  }                                                                            \
  void nca_read_from_file(const int, what &) {              \
    ARTS_USER_ERROR("NetCDF support not yet implemented for this type!");      \
  }

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY
