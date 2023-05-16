////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   nc_io.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-12

  \brief This file contains basic functions to handle NetCDF data files.
*/

#include "arts.h"

#ifdef ENABLE_NETCDF

#ifndef nc_io_h
#define nc_io_h

#include <netcdf.h>

#include "exceptions.h"
#include "messages.h"
#include "mystring.h"
#include "species_tags.h"

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void nca_filename(String& filename, const String& varname);

void nca_filename_with_index(String& filename,
                             const Index& file_index,
                             const String& varname);

////////////////////////////////////////////////////////////////////////////
//   Generic IO routines for XML files
////////////////////////////////////////////////////////////////////////////

template <typename T>
void nca_read_from_file(const String& filename,
                        T& type,
                        const Verbosity& verbosity);

template <typename T>
void nca_write_to_file(const String& filename,
                       const T& type,
                       const Verbosity& verbosity);

/*void nc_read_var(const int ncf, const int **ncvar,
                  const Index dims, const String& name);*/

void nca_def_dim(const int ncid,
                 const String& name,
                 const Index nelem,
                 int* ncdim);
void nca_def_var(const int ncid,
                 const String& name,
                 const nc_type type,
                 const int ndims,
                 const int* dims,
                 int* varid);

int nca_def_ArrayOfIndex(const int ncid,
                         const String& name,
                         const ArrayOfIndex& a);

int nca_def_Vector(const int ncid, const String& name, const Vector& v);

int nca_def_Matrix(const int ncid, const String& name, const Matrix& m);

int nca_def_Tensor4(const int ncid, const String& name, const Tensor4& t);

Index nca_get_dim(const int ncid, const String &name,
                  const bool noerror = false);

void nca_get_data(const int ncid, const String& name, int* data);

void nca_get_data(const int ncid, const String &name, long *data);

void nca_get_data(const int ncid, const String &name, long long *data);

void nca_get_data(const int ncid, const String &name, Numeric *data);

void nca_get_data(const int ncid, const String &name, size_t start,
                  size_t count, Numeric *data);

void nca_get_data(const int ncid, const String &name, char *data);

void nca_get_data(const int ncid, const String &name,
                               ArrayOfIndex &aoi, const bool noerror);

void nca_get_data(const int ncid, const String &name,
                  ArrayOfArrayOfSpeciesTag &aast, const bool noerror);

void nca_get_data(const int ncid, const String &name, Vector &v,
                  const bool noerror = false);

void nca_get_data(const int ncid, const String &name, Matrix &m,
                  const bool noerror = false);

void nca_get_data(const int ncid, const String &name, Tensor4 &m,
                  const bool noerror = false);

void nca_put_var(const int ncid, const int varid, const long *ind_arr);

void nca_put_var(const int ncid, const int varid, const long long*ind_arr);

bool nca_put_var(const int ncid, const int varid, const ArrayOfIndex &a);

bool nca_put_var(const int ncid, const int varid, const Vector& v);

bool nca_put_var(const int ncid, const int varid, const Matrix& m);

bool nca_put_var(const int ncid, const int varid, const Tensor4& t);

void nca_error(const int err, const std::string_view msg);

#endif /* nc_io_h */

#endif /* ENABLE_NETCDF */
