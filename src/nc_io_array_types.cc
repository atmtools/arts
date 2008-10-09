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
  \file   nc_io_array_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-10-01

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
  void nc_write_to_file (NcFile&, const what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  } \
  void nc_read_from_file (NcFile&, what&) \
  { \
    throw runtime_error ("NetCDF support not yet implemented for this type!"); \
  }

//=== Array Types ==========================================================

TMPL_NC_READ_WRITE_FILE_DUMMY( Array<IsotopeRecord> )
TMPL_NC_READ_WRITE_FILE_DUMMY( Array<SpeciesRecord> )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfIndex )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfLineRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfMatrix )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfSpeciesTag )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfTensor3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfArrayOfTensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField1 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField2 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGField4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfGridPos )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfIndex )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfLineRecord )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfLineshapeSpec )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfMatrix )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfPpath )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfRetrievalQuantity )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfSingleScatteringData )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfSpeciesTag )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfString )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor3 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor4 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor6 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfTensor7 )
TMPL_NC_READ_WRITE_FILE_DUMMY( ArrayOfVector )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE_DUMMY


