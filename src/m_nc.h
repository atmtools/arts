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

/*!
  \file   m_nc.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-19

  \brief  Workspace methods and template functions for supergeneric NetCDF IO.

*/

#ifndef m_nc_h
#define m_nc_h

#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Please run ./configure in the top arts directory before compiling."
#endif

#include "exceptions.h"
#include "workspace_ng.h"
#include "agenda_class.h"


#ifdef ENABLE_NETCDF

#include "nc_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
ReadNetCDF (// WS Generic Input:
            T&            v,
            const String& v_name _U_,
            const String& f,
            // WS Generic Input Names:
            const String& f_name _U_)

{
  nc_read_from_file (f, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
WriteNetCDF (// WS Generic Input:
             const T&      v,
             const String& f,
             // WS Generic Input Names:
             const String& v_name,
             const String& f_name _U_)

{
  String filename = f;

  // Create default filename if empty
  filename_nc (filename, v_name);

  nc_write_to_file (filename, v);
}

#else // NetCDF not enabled

/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
ReadNetCDF (// WS Generic Input:
            T&            v _U_,
            const String& v_name _U_,
            const String& f _U_,
            // WS Generic Input Names:
            const String& f_name _U_)

{
  throw runtime_error("This version of arts was compiled without NetCDF support.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
WriteNetCDF (// WS Generic Input:
             const T&      v _U_,
             const String& f _U_,
             // WS Generic Input Names:
             const String& v_name _U_,
             const String& f_name _U_)

{
  throw runtime_error("This version of arts was compiled without NetCDF support.");
}

#endif // ENABLE_NETCDF

/* Workspace method: Doxygen documentation will be auto-generated */
template<typename T> void
ReadNetCDF (Workspace& ws _U_,
            // WS Generic Input:
            T&            v,
            const String& v_name,
            const String& f,
            // WS Generic Input Names:
            const String& f_name)

{
  ReadNetCDF (v, f, v_name, f_name);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void
WriteNetCDF (Workspace& ws _U_,
             // WS Generic Input:
             const Agenda& v,
             const String& f,
             // WS Generic Input Names:
             const String& v_name,
             const String& f_name)
{
  WriteNetCDF (v, f, v_name, f_name);
}

#endif // m_nc_h

