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
  \file   xml_io_compound_types.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_compound_types_h
#define xml_io_compound_types_h

#include "messages.h"
#include "ppath.h"
#include "agenda_class.h"
#include "absorption.h"
#include "gas_abs_lookup.h"
#include "optproperties.h"
#include "gridded_fields.h"
#include "bifstream.h"
#include "bofstream.h"
#include "jacobian.h"
#include "mc_interp.h"
#include "mc_antenna.h"

#define TMPL_XML_READ_WRITE_STREAM(what) \
  void xml_read_from_stream (istream&, what&, bifstream * = NULL); \
  void xml_write_to_stream (ostream&, const what&, bofstream * = NULL, \
                            const String& = "");


////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

//=== Compound Types =======================================================

TMPL_XML_READ_WRITE_STREAM( Agenda )
TMPL_XML_READ_WRITE_STREAM( GField1 )
TMPL_XML_READ_WRITE_STREAM( GField2 )
TMPL_XML_READ_WRITE_STREAM( GField3 )
TMPL_XML_READ_WRITE_STREAM( GField4 )
TMPL_XML_READ_WRITE_STREAM( GasAbsLookup )
TMPL_XML_READ_WRITE_STREAM( GridPos )
TMPL_XML_READ_WRITE_STREAM( IsotopeRecord )
TMPL_XML_READ_WRITE_STREAM( MCAntenna )
TMPL_XML_READ_WRITE_STREAM( Ppath )
TMPL_XML_READ_WRITE_STREAM( RetrievalQuantity )
TMPL_XML_READ_WRITE_STREAM( SLIData2 )
TMPL_XML_READ_WRITE_STREAM( SingleScatteringData )
TMPL_XML_READ_WRITE_STREAM( SpeciesRecord )
TMPL_XML_READ_WRITE_STREAM( SpeciesTag )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE_STREAM


#endif  /* xml_io_compound_types_h */
