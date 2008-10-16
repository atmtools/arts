/* Copyright (C) 2003-2008 Oliver Lemke <olemke@core-dump.info>

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
  \file   nc_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-26

  \brief This file contains functions to handle NetCDF data files.

*/

#include "arts.h"
#include "nc_io.h"
#include "nc_io_types.h"

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

#define TMPL_NC_READ_WRITE_FILE_DUMMY(what) \
  void nc_write_to_file (const int, const what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  } \
  void nc_read_from_file (const int, what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  }

//=== Compound Types =======================================================

TMPL_NC_READ_WRITE_FILE_DUMMY( Agenda )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GField4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( GasAbsLookup )
TMPL_NC_READ_WRITE_FILE_DUMMY( GridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( IsotopeRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( MCAntenna )
TMPL_NC_READ_WRITE_FILE_DUMMY( Ppath )
TMPL_NC_READ_WRITE_FILE_DUMMY( RetrievalQuantity )
TMPL_NC_READ_WRITE_FILE_DUMMY( SLIData2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( SingleScatteringData )
TMPL_NC_READ_WRITE_FILE_DUMMY( SpeciesRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( SpeciesTag )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY


