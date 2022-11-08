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
  \file   nc_io_types.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle NetCDF data files.

*/

#ifndef nc_io_instantiation_h
#define nc_io_instantiation_h

#include <cfloat>
#include <stdexcept>
#include "nc_io.h"
#include "nc_io_types.h"

#define TMPL_NC_READ_WRITE_FILE(what)                \
  template void nca_write_to_file<what>(             \
      const String&, const what&, const Verbosity&); \
  template void nca_read_from_file<what>(            \
      const String&, what&, const Verbosity&);

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for NetCDF streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_NC_READ_WRITE_FILE(Matrix)
TMPL_NC_READ_WRITE_FILE(Tensor3)
TMPL_NC_READ_WRITE_FILE(Tensor4)
TMPL_NC_READ_WRITE_FILE(Tensor5)
TMPL_NC_READ_WRITE_FILE(Vector)

//=== Compound Types =======================================================

TMPL_NC_READ_WRITE_FILE(Agenda)
TMPL_NC_READ_WRITE_FILE(GasAbsLookup)

//=== Array Types ==========================================================

TMPL_NC_READ_WRITE_FILE(ArrayOfIndex)
TMPL_NC_READ_WRITE_FILE(ArrayOfMatrix)
TMPL_NC_READ_WRITE_FILE(ArrayOfVector)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE

/*void
xml_parse_from_stream (istream&, Vector&, bifstream *, ArtsXMLTag&, const Verbosity&);

void
xml_read_from_stream (istream&, ArrayOfLineRecord&,
                      const Numeric, const Numeric, bifstream * = NULL,
                      const Verbosity&);

void
xml_parse_from_stream (istream&, ArrayOfString&, bifstream *, ArtsXMLTag&, const Verbosity&);*/

#endif /* nc_io_instantiation_h */
