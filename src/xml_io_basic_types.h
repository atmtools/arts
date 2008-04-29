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
  \file   xml_io_basic_types.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_basic_types_h
#define xml_io_basic_types_h

#include <stdexcept>
#include <cfloat>
#include "matpackII.h"
#include "matpackVII.h"
#include "array.h"
#include "messages.h"
#include "ppath.h"
#include "agenda_class.h"
#include "absorption.h"
#include "gas_abs_lookup.h"
#include "optproperties.h"
#include "m_general.h"
#include "bifstream.h"
#include "bofstream.h"


////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

void
xml_read_from_stream (istream&, Index&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Index&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Matrix&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Matrix&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Numeric&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Numeric&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Sparse&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Sparse&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, String&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const String&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Tensor3&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Tensor3&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Tensor4&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Tensor4&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Tensor5&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Tensor5&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Tensor6&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Tensor6&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Tensor7&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Tensor7&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Timer&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Timer&, bofstream * = NULL,
                     const String& = "");

void
xml_read_from_stream (istream&, Vector&, bifstream * = NULL);

void
xml_write_to_stream (ostream&, const Vector&, bofstream * = NULL,
                     const String& = "");


#endif  /* xml_io_basic_types_h */
