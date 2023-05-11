/* Copyright (C) 2012 Oliver Lemke <olemke@core-dump.info>

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

/*!
  \file   m_nc.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-19

  \brief  Workspace methods and template functions for supergeneric NetCDF IO.

*/

#ifndef m_nc_h
#define m_nc_h

#include "arts.h"

#include "agenda_class.h"
#include "exceptions.h"
#include "workspace_ng.h"

#ifdef ENABLE_NETCDF

#include "nc_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadNetCDF(  // WS Generic Input:
    T& v,
    const String& v_name _U_,
    const String& f,
    // WS Generic Input Names:
    const String& f_name _U_)

{
  nca_read_from_file(f, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDF(  // WS Generic Input:
    const T& v,
    const String& f,
    // WS Generic Input Names:
    const String& v_name,
    const String& f_name _U_)

{
  String filename = f;

  // Create default filename if empty
  nca_filename(filename, v_name);

  nca_write_to_file(filename, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDFIndexed(  //WS Input:
    const Index& file_index,
    // WS Generic Input:
    const T& v,
    const String& f,
    // WS Generic Input Names:
    const String& v_name,
    const String& f_name)

{
  String filename = f;

  // Create default filename if empty
  nca_filename_with_index(filename, file_index, v_name);

  WriteNetCDF(v, filename, v_name, f_name);
}

#else  // NetCDF not enabled

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadNetCDF(  // WS Generic Input:
    T&,
    const String&,
    const String&,
    // WS Generic Input Names:
    const String&)

{
  throw runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDF(  // WS Generic Input:
    const T&,
    const String&,
    // WS Generic Input Names:
    const String&,
    const String&)

{
  throw runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDFIndexed(  //WS Input:
    const Index&,
    // WS Generic Input:
    const T&,
    const String&,
    // WS Generic Input Names:
    const String&,
    const String&)

{
  throw runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
}

#endif  // ENABLE_NETCDF

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadNetCDF(Workspace& ws _U_,
                // WS Generic Input:
                T& v,
                const String& v_name,
                const String& f,
                // WS Generic Input Names:
                const String& f_name)

{
  ReadNetCDF(v, f, v_name, f_name);
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void WriteNetCDF(Workspace& ws _U_,
                 // WS Generic Input:
                 const Agenda& v,
                 const String& f,
                 // WS Generic Input Names:
                 const String& v_name,
                 const String& f_name) {
  WriteNetCDF(v, f, v_name, f_name);
}

#endif  // m_nc_h
