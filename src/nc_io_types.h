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

#include "arts.h"

#ifdef ENABLE_NETCDF

#ifndef nc_io_types_h
#define nc_io_types_h

#include <stdexcept>
#include <cfloat>
#include "nc_io.h"
#include "absorption.h"
#include "agenda_class.h"
#include "array.h"
#include "bifstream.h"
#include "bofstream.h"
#include "gas_abs_lookup.h"
#include "gridded_fields.h"
#include "jacobian.h"
#include "m_general.h"
#include "mc_antenna.h"
#include "mc_interp.h"
#include "matpackII.h"
#include "matpackVII.h"
#include "messages.h"
#include "optproperties.h"
#include "ppath.h"


#define TMPL_NC_READ_WRITE_FILE(what) \
  void nca_write_to_file(const int, const what&, const Verbosity&); \
  void nca_read_from_file(const int, what&, const Verbosity&);


////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for NetCDF streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_NC_READ_WRITE_FILE(Index)
TMPL_NC_READ_WRITE_FILE(Matrix)
TMPL_NC_READ_WRITE_FILE(Numeric)
TMPL_NC_READ_WRITE_FILE(Sparse)
TMPL_NC_READ_WRITE_FILE(String)
TMPL_NC_READ_WRITE_FILE(Tensor3)
TMPL_NC_READ_WRITE_FILE(Tensor4)
TMPL_NC_READ_WRITE_FILE(Tensor5)
TMPL_NC_READ_WRITE_FILE(Tensor6)
TMPL_NC_READ_WRITE_FILE(Tensor7)
TMPL_NC_READ_WRITE_FILE(Timer)
TMPL_NC_READ_WRITE_FILE(Vector)

//=== Compound Types =======================================================

TMPL_NC_READ_WRITE_FILE(Agenda)
TMPL_NC_READ_WRITE_FILE(GriddedField1)
TMPL_NC_READ_WRITE_FILE(GriddedField2)
TMPL_NC_READ_WRITE_FILE(GriddedField3)
TMPL_NC_READ_WRITE_FILE(GriddedField4)
TMPL_NC_READ_WRITE_FILE(GasAbsLookup)
TMPL_NC_READ_WRITE_FILE(GridPos)
TMPL_NC_READ_WRITE_FILE(IsotopologueRecord)
TMPL_NC_READ_WRITE_FILE(MCAntenna)
TMPL_NC_READ_WRITE_FILE(Ppath)
TMPL_NC_READ_WRITE_FILE(RetrievalQuantity)
TMPL_NC_READ_WRITE_FILE(SLIData2)
TMPL_NC_READ_WRITE_FILE(SingleScatteringData)
TMPL_NC_READ_WRITE_FILE(SpeciesAuxData)
TMPL_NC_READ_WRITE_FILE(SpeciesRecord)
TMPL_NC_READ_WRITE_FILE(SpeciesTag)

//=== Array Types ==========================================================

TMPL_NC_READ_WRITE_FILE(Array<IsotopologueRecord> )
TMPL_NC_READ_WRITE_FILE(Array<SpeciesRecord> )
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfArrayOfArrayOfGridPos)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfGriddedField1)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfGriddedField3)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfGridPos)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfArrayOfGridPos)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfIndex)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfLineRecord)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfMatrix)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfSpeciesTag)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfTensor3)
TMPL_NC_READ_WRITE_FILE(ArrayOfArrayOfTensor6)
TMPL_NC_READ_WRITE_FILE(ArrayOfGriddedField1)
TMPL_NC_READ_WRITE_FILE(ArrayOfGriddedField2)
TMPL_NC_READ_WRITE_FILE(ArrayOfGriddedField3)
TMPL_NC_READ_WRITE_FILE(ArrayOfGriddedField4)
TMPL_NC_READ_WRITE_FILE(ArrayOfGridPos)
TMPL_NC_READ_WRITE_FILE(ArrayOfIndex)
TMPL_NC_READ_WRITE_FILE(ArrayOfLineRecord)
TMPL_NC_READ_WRITE_FILE(ArrayOfLineshapeSpec)
TMPL_NC_READ_WRITE_FILE(ArrayOfMatrix)
TMPL_NC_READ_WRITE_FILE(ArrayOfPpath)
TMPL_NC_READ_WRITE_FILE(ArrayOfRetrievalQuantity)
TMPL_NC_READ_WRITE_FILE(ArrayOfSingleScatteringData)
TMPL_NC_READ_WRITE_FILE(ArrayOfSpeciesTag)
TMPL_NC_READ_WRITE_FILE(ArrayOfString)
TMPL_NC_READ_WRITE_FILE(ArrayOfTensor3)
TMPL_NC_READ_WRITE_FILE(ArrayOfTensor4)
TMPL_NC_READ_WRITE_FILE(ArrayOfTensor6)
TMPL_NC_READ_WRITE_FILE(ArrayOfTensor7)
TMPL_NC_READ_WRITE_FILE(ArrayOfVector)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE

/*void
xml_parse_from_stream (istream&, Vector&, bifstream *, ArtsXMLTag&, const Verbosity&);

void
xml_read_from_stream (istream&, ArrayOfLineRecord&,
                      const Numeric, const Numeric, bifstream * = NULL, const Verbosity&);

void
xml_parse_from_stream (istream&, ArrayOfString&, bifstream *, ArtsXMLTag&, const Verbosity&);*/

#endif  /* nc_io_types_h */

#endif /* ENABLE_NETCDF */

