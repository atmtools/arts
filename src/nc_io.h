/* Copyright (C) 2008 Oliver Lemke <olemke@core-dump.info>

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
#include "mystring.h"
#include "exceptions.h"
#include "messages.h"
#include "abs_species_tags.h"


////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_nc(String& filename, const String&  varname);

void filename_nc_with_index(String& filename, const Index& file_index, const String& varname);


////////////////////////////////////////////////////////////////////////////
//   Generic IO routines for XML files
////////////////////////////////////////////////////////////////////////////

template<typename T>
void nc_read_from_file(const String& filename, T& type, const Verbosity& verbosity);


template<typename T>
void nc_write_to_file(const String&  filename, const T&  type, const Verbosity& verbosity);


/*void nc_read_var(const int ncf, const int **ncvar,
                  const Index dims, const String& name);*/

void nc_get_data_int(const int ncid, const String &name, int *data);

void nc_get_data_long(const int ncid, const String &name, long *data);

void nc_get_data_double(const int ncid, const String &name, Numeric *data);

void nc_get_dataa_double(const int ncid, const String &name,
                         size_t start, size_t count, Numeric *data);

void nc_get_data_text(const int ncid, const String &name, char *data);

Index nc_get_dim(const int ncid, const String &name, const bool noerror = false);

void nc_get_data_ArrayOfIndex(const int ncid, const String &name, ArrayOfIndex &aoi,
                              const bool noerror);

void nc_get_data_ArrayOfArrayOfSpeciesTag(const int ncid, const String &name,
                                          ArrayOfArrayOfSpeciesTag &aast,
                                          const bool noerror);

void nc_get_data_Vector(const int ncid, const String &name, Vector &v, const bool noerror = false);

void nc_get_data_Matrix(const int ncid, const String &name, Matrix &m, const bool noerror = false);

void nc_get_data_Tensor4(const int ncid, const String &name, Tensor4 &m, const bool noerror = false);

void ncerror(const int err, const String msg);

#endif /* nc_io_h */

#endif /* ENABLE_NETCDF */

