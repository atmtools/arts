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
  \file   nc_io.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-12

  \brief This file contains basic functions to handle NetCDF data files.

*/

#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Please run ./configure in the top arts directory before compiling."
#endif

#ifdef ENABLE_NETCDF

#include "arts.h"
#include "nc_io.h"
#include "nc_io_types.h"
#include "nc_io_instantiation.h"
#include "file.h"
#include "messages.h"
#include "exceptions.h"


////////////////////////////////////////////////////////////////////////////
//   Default file name
////////////////////////////////////////////////////////////////////////////

//! Gives the default filename for the NetCDF formats.
/*!
  The default name is only used if the filename is empty.

  \param filename filename
  \param varname variable name
*/
void
filename_nc (      String&  filename,
             const String&  varname )
{
  if ("" == filename)
    {
      extern const String out_basename;
      filename = out_basename + "." + varname + ".nc";
    }
}



//! Gives the default filename, with file index, for the NetCDF formats.
/*!
  The default name is only used if the filename is empty.

  \param[out] filename   filename
  \param[in]  file_index Index appended to the filename
  \param[in]  varname    variable name
*/
void
filename_nc_with_index (
                    String&  filename,
              const Index&   file_index,
              const String&  varname )
{
  Index dummy = file_index;
  dummy += 1;

  if ("" == filename)
    {
      extern const String out_basename;
      ostringstream os;
      os << out_basename << "." << varname << "." << file_index << ".nc";
      filename = os.str();
    }
  else
    {
      ostringstream os;
      os << filename << "." << file_index << ".nc";
      filename = os.str();
    }
}


template<typename T> void
nc_read_from_file (const String& filename,
                         T&      type)
{
  String efilename = expand_path(filename);
  
  out2 << "  Reading " << efilename << '\n';

  int ncid;
  if (nc_open (efilename.c_str(), NC_NOWRITE, &ncid))
    {
      ostringstream os;
      os << "Error reading file: " << efilename << endl;
      throw runtime_error (os.str());
    }

  nc_read_from_file (ncid, type);

  nc_close (ncid);
}


template<typename T> void
nc_write_to_file (const String&  filename,
                  const      T&  type)
{
  String efilename = expand_path(filename);
  
  out2 << "  Writing " << efilename << '\n';

  int ncid;
  if (nc_create (efilename.c_str(), NC_CLOBBER | NC_FORMAT_NETCDF4, &ncid))
    {
      ostringstream os;
      os << "Error writing file: " << efilename << endl;
      throw runtime_error (os.str());
    }

  nc_write_to_file (ncid, type);

  nc_close (ncid);
}


void nc_get_data_int (const int ncid, const String &name, int *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_var_int (ncid, varid, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


void nc_get_data_long (const int ncid, const String &name, long *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_var_long (ncid, varid, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


void nc_get_data_double (const int ncid, const String &name, Numeric *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_var_double (ncid, varid, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


void nc_get_dataa_double (const int ncid, const String &name,
                          size_t start, size_t count, Numeric *data)
{
  int retval, varid;
  if ((retval = nc_inq_varid (ncid, name.c_str(), &varid)))
    ncerror (retval, "nc_inq_varid("+name+")");
  if ((retval = nc_get_vara_double (ncid, varid, &start, &count, data)))
    ncerror (retval, "nc_get_var("+name+")");
}


Index nc_get_dim (const int ncid, const String &name)
{
  int retval, dimid;
  size_t ndim;
  if ((retval = nc_inq_dimid (ncid, name.c_str(), &dimid)))
    ncerror (retval, "nc_inq_ndims("+name+")");
  if ((retval = nc_inq_dimlen (ncid, dimid, &ndim)))
    ncerror (retval, "nc_inq_dimlen("+name+")");

  return (Index)ndim;
}


void ncerror (const int e, const String s)
{
  ostringstream os;
  os << "NetCDF error: " << s << ", " << e;
  throw runtime_error (os.str());
}

#endif /* ENABLE_NETCDF */

