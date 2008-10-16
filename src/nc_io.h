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


#ifndef nc_io_h
#define nc_io_h

#include <netcdf.h>
#include "mystring.h"
#include "exceptions.h"
#include "messages.h"


////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void
filename_nc (      String&  filename,
             const String&  varname);

void
filename_nc_with_index (       String&  filename,
                         const Index&   file_index,
                         const String&  varname );


////////////////////////////////////////////////////////////////////////////
//   Generic IO routines for XML files
////////////////////////////////////////////////////////////////////////////

template<typename T> void
nc_read_from_file (const String& filename _U_,
                         T&      type _U_);


template<typename T> void
nc_write_to_file (const String&  filename,
                  const      T&  type);


/*void nc_read_var (const int ncf, const int **ncvar,
                  const Index dims, const String& name);*/

void nc_get_data_int (const int ncid, const String &name, int *data);

void nc_get_data_double (const int ncid, const String &name, Numeric *data);

void nc_get_dataa_double (const int ncid, const String &name,
                          size_t start, size_t count, Numeric *data);

Index nc_get_dim (const int ncid, const String &name);

void ncerror (const int err, const String msg);

#endif /* nc_io_h */

