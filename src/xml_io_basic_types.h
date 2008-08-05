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
#include "xml_io_private.h"

#define TMPL_XML_READ_WRITE_STREAM(what) \
  void xml_read_from_stream (istream&, what&, bifstream * = NULL); \
  void xml_write_to_stream (ostream&, const what&, bofstream * = NULL, \
                            const String& = "");


////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_XML_READ_WRITE_STREAM( Index )
TMPL_XML_READ_WRITE_STREAM( Matrix )
TMPL_XML_READ_WRITE_STREAM( Numeric )
TMPL_XML_READ_WRITE_STREAM( Sparse )
TMPL_XML_READ_WRITE_STREAM( String )
TMPL_XML_READ_WRITE_STREAM( Tensor3 )
TMPL_XML_READ_WRITE_STREAM( Tensor4 )
TMPL_XML_READ_WRITE_STREAM( Tensor5 )
TMPL_XML_READ_WRITE_STREAM( Tensor6 )
TMPL_XML_READ_WRITE_STREAM( Tensor7 )
TMPL_XML_READ_WRITE_STREAM( Timer )
TMPL_XML_READ_WRITE_STREAM( Vector )

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE_STREAM


void
xml_parse_from_stream (istream&, Vector&, bifstream *, ArtsXMLTag&);


#endif  /* xml_io_basic_types_h */
