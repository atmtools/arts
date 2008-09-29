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

#include <netcdfcpp.h>
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
                         T&      type _U_)
{
  throw runtime_error("NetCDF reading not implemented for this type yet.");
}


template<typename T> void
nc_write_to_file (const String&  filename,
                  const      T&  type)
{
  out2 << "  Writing " << filename << '\n';

  NcFile ncf(filename.c_str(), NcFile::Replace);

  if (!ncf.is_valid())
    {
      ostringstream os;
      os << "Error writing file: " << filename << endl;
      throw runtime_error (os.str());
    }

  nc_write_to_file (ncf, type);
}


template<typename T> void
nc_write_to_file (const NcFile&,
                  const T&)
{
  throw runtime_error("NetCDF writing not implemented for this type yet.");
}

#endif /* nc_io_h */
