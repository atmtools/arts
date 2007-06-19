/* Copyright (C) 2003-2007 Oliver Lemke <olemke@core-dump.info>

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

#ifndef xml_io_compound_h
#define xml_io_compound_h

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

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

void
xml_read_from_stream (istream&, Agenda&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Agenda&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, GasAbsLookup&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const GasAbsLookup&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, GriddedField3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const GriddedField3&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, GridPos&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const GridPos&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, IsotopeRecord&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const IsotopeRecord&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, MCAntenna&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const MCAntenna&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Ppath&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Ppath&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, RetrievalQuantity&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const RetrievalQuantity&,
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, SingleScatteringData&,
                      bifstream * = NULL);

void
xml_write_to_stream (ostream&, const SingleScatteringData&,
                     bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, SLIData2&,
                      bifstream * = NULL);
void
xml_write_to_stream (ostream&, const SLIData2&,
                     bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, SpeciesRecord&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const SpeciesRecord&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, SpeciesTag&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const SpeciesTag&, bofstream * = NULL,
                     const String& = "");

#endif  /* xml_io_compound_types_h */
