/* Copyright (C) 2003 Oliver Lemke <olemke@uni-bremen.de>

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
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_array_types_h
#define xml_io_array_types_h

#include "matpackI.h"
#include "matpackII.h"
#include "matpackIII.h"
#include "matpackIV.h"
#include "matpackV.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "array.h"
#include "messages.h"
#include "ppath.h"
#include "agenda_class.h"
#include "absorption.h"
#include "gas_abs_lookup.h"
#include "optproperties.h"
#include "bifstream.h"
#include "bofstream.h"


////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

void
xml_read_from_stream (istream&, Array<SpeciesRecord>&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Array<SpeciesRecord>&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfArrayOfSpeciesTag&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfSpeciesTag&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfSingleScatteringData&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfSingleScatteringData&,
                     bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfArrayOfTensor3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfTensor3&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfArrayOfTensor6&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfTensor6&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfGridPos&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfGridPos&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfIndex&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfIndex&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfMatrix&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfMatrix&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfSpeciesTag&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfSpeciesTag&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfString&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfString&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfTensor3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfTensor3&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfTensor6&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfTensor6&, bofstream * = NULL);

void
xml_read_from_stream (istream&, ArrayOfVector&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const ArrayOfVector&, bofstream * = NULL);

#endif  /* xml_io_array_types_h */
