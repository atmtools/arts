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
  \file   xml_io_array_types.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_array_types_h
#define xml_io_array_types_h

#include "matpackVII.h"
#include "array.h"
#include "messages.h"
#include "ppath.h"
#include "agenda_class.h"
#include "absorption.h"
#include "gas_abs_lookup.h"
#include "gridded_fields.h"
#include "optproperties.h"
#include "bifstream.h"
#include "bofstream.h"
#include "jacobian.h"


////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

void
xml_read_from_stream (istream&, Array<IsotopeRecord>&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Array<IsotopeRecord>&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Array<SpeciesRecord>&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Array<SpeciesRecord>&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfSpeciesTag&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfSpeciesTag&,
                     bofstream * = NULL, const String& = "");
void
xml_read_from_stream (istream&, ArrayOfSingleScatteringData&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfSingleScatteringData&,
                     bofstream * = NULL, const String & = "");

void
xml_read_from_stream (istream&, ArrayOfGriddedField3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfGriddedField3&,
                     bofstream * = NULL, const String & = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfGriddedField3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfGriddedField3&,
                     bofstream * = NULL, const String & = "");

void
xml_read_from_stream (istream&, ArrayOfGriddedField4&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfGriddedField4&,
                     bofstream * = NULL, const String & = "");

void
xml_read_from_stream (istream&, ArrayOfLineRecord&, bifstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfLineRecord&,
                      const Numeric, const Numeric, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfLineRecord&,
                     bofstream * = NULL, const String & = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfLineRecord&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfLineRecord&,
                     bofstream * = NULL, const String & = "");

void
xml_read_from_stream (istream&, ArrayOfLineshapeSpec&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfLineshapeSpec&,
                     bofstream * = NULL, const String & = "");

void
xml_read_from_stream (istream&, ArrayOfPpath&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfPpath&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfTensor3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfTensor3&,
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfTensor6&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfTensor6&,
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfGridPos&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfGridPos&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfGridPos&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfGridPos&, 
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfArrayOfGridPos&, 
                      bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfArrayOfGridPos&, 
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfArrayOfArrayOfGridPos&,
                      bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfArrayOfArrayOfGridPos&, 
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfIndex&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfIndex&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfArrayOfIndex&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfIndex&,
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfMatrix&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfMatrix&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&,ArrayOfArrayOfMatrix&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfMatrix&,
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfRetrievalQuantity&,
                      bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfRetrievalQuantity&,
                     bofstream * = NULL, const String& = "");

void
xml_read_from_stream (istream&, ArrayOfSpeciesTag&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfSpeciesTag&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfString&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfString&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfTensor3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfTensor3&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfTensor4&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfTensor4&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfTensor6&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfTensor6&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfTensor7&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfTensor7&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, ArrayOfVector&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfVector&, bofstream * = NULL,
                     const String& = "");

#endif  /* xml_io_array_types_h */
