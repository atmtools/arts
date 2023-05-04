////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io_types_matpack.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_general_types_h
#define xml_io_general_types_h

#include <cfloat>
#include <stdexcept>
#include "array.h"
#include "bifstream.h"
#include "bofstream.h"
#include "matpack_data.h"
#include "matpack_sparse.h"
#include "messages.h"

#define TMPL_XML_READ_WRITE_STREAM(what)                  \
  void xml_read_from_stream(                              \
      istream &, what &, bifstream *, const Verbosity &); \
  void xml_write_to_stream(ostream &,                     \
                           const what &,                  \
                           bofstream *,                   \
                           const String &,                \
                           const Verbosity &);

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_XML_READ_WRITE_STREAM(Index)
TMPL_XML_READ_WRITE_STREAM(Matrix)
TMPL_XML_READ_WRITE_STREAM(Numeric)
TMPL_XML_READ_WRITE_STREAM(Sparse)
TMPL_XML_READ_WRITE_STREAM(String)
TMPL_XML_READ_WRITE_STREAM(Tensor3)
TMPL_XML_READ_WRITE_STREAM(Tensor4)
TMPL_XML_READ_WRITE_STREAM(Tensor5)
TMPL_XML_READ_WRITE_STREAM(Tensor6)
TMPL_XML_READ_WRITE_STREAM(Tensor7)
TMPL_XML_READ_WRITE_STREAM(Vector)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE_STREAM

#endif /* xml_io_general_types_h */
