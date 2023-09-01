////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   nc_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-26

  \brief This file contains functions to handle NetCDF data files.

*/

#include "nc_io.h"
#include "nc_io_types.h"

//=== Matrix ==========================================================

//! Reads a Matrix from a NetCDF file
/*!
  \param ncdi    NetCDF file descriptor
  \param m       Matrix
*/
void nca_read_from_file(const int ncid, Matrix& m) {
  Index nrows, ncols;
  nrows = nca_get_dim(ncid, "nrows");
  ncols = nca_get_dim(ncid, "ncols");

  m.resize(nrows, ncols);
  nca_get_data(ncid, "Matrix", m.unsafe_data_handle());
}

//! Writes a Matrix to a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param m       Matrix
*/
void nca_write_to_file(const int ncid, const Matrix& m) {
  int retval;
  int ncdims[2], varid;
  if ((retval = nc_def_dim(ncid, "nrows", m.nrows(), &ncdims[0])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "ncols", m.ncols(), &ncdims[1])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_var(ncid, "Matrix", NC_DOUBLE, 2, &ncdims[0], &varid)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");
  if ((retval = nc_put_var_double(ncid, varid, m.unsafe_data_handle())))
    nca_error(retval, "nc_put_var");
}

//=== Tensor3 ==========================================================

//! Reads a Tensor3 from a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param t       Tensor3
*/
void nca_read_from_file(const int ncid, Tensor3& t) {
  Index npages, nrows, ncols;
  npages = nca_get_dim(ncid, "npages");
  nrows = nca_get_dim(ncid, "nrows");
  ncols = nca_get_dim(ncid, "ncols");

  t.resize(npages, nrows, ncols);
  nca_get_data(ncid, "Tensor3", t.unsafe_data_handle());
}

//! Writes a Tensor3 to a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param t       Tensor3
*/
void nca_write_to_file(const int ncid, const Tensor3& t) {
  int retval;
  int ncdims[3], varid;
  if ((retval = nc_def_dim(ncid, "npages", t.npages(), &ncdims[0])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "nrows", t.nrows(), &ncdims[1])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "ncols", t.ncols(), &ncdims[2])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_var(ncid, "Tensor3", NC_DOUBLE, 3, &ncdims[0], &varid)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");
  if ((retval = nc_put_var_double(ncid, varid, t.unsafe_data_handle())))
    nca_error(retval, "nc_put_var");
}

//=== Tensor4 ==========================================================

//! Reads a Tensor4 from a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param t       Tensor4
*/
void nca_read_from_file(const int ncid, Tensor4& t) {
  Index nbooks, npages, nrows, ncols;
  nbooks = nca_get_dim(ncid, "nbooks");
  npages = nca_get_dim(ncid, "npages");
  nrows = nca_get_dim(ncid, "nrows");
  ncols = nca_get_dim(ncid, "ncols");

  t.resize(nbooks, npages, nrows, ncols);
  nca_get_data(ncid, "Tensor4", t.unsafe_data_handle());
}

//! Writes a Tensor4 to a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param t       Tensor4
*/
void nca_write_to_file(const int ncid, const Tensor4& t) {
  int retval;
  int ncdims[4], varid;
  if ((retval = nc_def_dim(ncid, "nbooks", t.nbooks(), &ncdims[0])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "npages", t.npages(), &ncdims[1])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "nrows", t.nrows(), &ncdims[2])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "ncols", t.ncols(), &ncdims[3])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_var(ncid, "Tensor4", NC_DOUBLE, 4, &ncdims[0], &varid)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");
  if ((retval = nc_put_var_double(ncid, varid, t.unsafe_data_handle())))
    nca_error(retval, "nc_put_var");
}

//=== Tensor5 ==========================================================

//! Reads a Tensor5 from a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param t       Tensor5
*/
void nca_read_from_file(const int ncid, Tensor5& t) {
  Index nshelves, nbooks, npages, nrows, ncols;
  nshelves = nca_get_dim(ncid, "nshelves");
  nbooks = nca_get_dim(ncid, "nbooks");
  npages = nca_get_dim(ncid, "npages");
  nrows = nca_get_dim(ncid, "nrows");
  ncols = nca_get_dim(ncid, "ncols");

  t.resize(nshelves, nbooks, npages, nrows, ncols);
  nca_get_data(ncid, "Tensor5", t.unsafe_data_handle());
}

//! Writes a Tensor5 to a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param t       Tensor5
*/
void nca_write_to_file(const int ncid, const Tensor5& t) {
  int retval;
  int ncdims[5], varid;
  if ((retval = nc_def_dim(ncid, "nshelves", t.nshelves(), &ncdims[0])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "nbooks", t.nbooks(), &ncdims[1])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "npages", t.npages(), &ncdims[2])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "nrows", t.nrows(), &ncdims[3])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_dim(ncid, "ncols", t.ncols(), &ncdims[4])))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_var(ncid, "Tensor5", NC_DOUBLE, 5, &ncdims[0], &varid)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");
  if ((retval = nc_put_var_double(ncid, varid, t.unsafe_data_handle())))
    nca_error(retval, "nc_put_var");
}

//=== Vector ==========================================================

//! Reads a Vector from a NetCDF file
/*!
  \param ncf     NetCDF file descriptor
  \param v       Vector
*/
void nca_read_from_file(const int ncid, Vector& v) {
  Index nelem;
  nelem = nca_get_dim(ncid, "nelem");

  v.resize(nelem);
  nca_get_data(ncid, "Vector", v.unsafe_data_handle());
}

//! Writes a Vector to a NetCDF file
/*!
  \param ncid    NetCDF file descriptor
  \param v       Vector
*/
void nca_write_to_file(const int ncid, const Vector& v) {
  int retval;
  int ncdim, varid;
  if ((retval = nc_def_dim(ncid, "nelem", v.nelem(), &ncdim)))
    nca_error(retval, "nc_def_dim");
  if ((retval = nc_def_var(ncid, "Vector", NC_DOUBLE, 1, &ncdim, &varid)))
    nca_error(retval, "nc_def_var");
  if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");
  if ((retval = nc_put_var_double(ncid, varid, v.unsafe_data_handle())))
    nca_error(retval, "nc_put_var");
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what)                               \
  void nca_write_to_file(const int, const what&) {                        \
    ARTS_USER_ERROR("NetCDF support not yet implemented for this type!"); \
  }                                                                       \
  void nca_read_from_file(const int, what&) {                             \
    ARTS_USER_ERROR("NetCDF support not yet implemented for this type!"); \
  }

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY
