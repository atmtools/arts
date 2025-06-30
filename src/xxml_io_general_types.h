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

#include <bifstream.h>
#include <bofstream.h>
#include <matpack.h>
#include <mystring.h>
#include <supergeneric.h>

#include <cfloat>
#include <stdexcept>

#define TMPL_XML_READ_WRITE_STREAM(what)                          \
  void xml_read_from_stream(std::istream &, what &, bifstream *); \
  void xml_write_to_stream(                                       \
      std::ostream &, const what &, bofstream *, const String &);

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_XML_READ_WRITE_STREAM(Any)
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
TMPL_XML_READ_WRITE_STREAM(Vector2)
TMPL_XML_READ_WRITE_STREAM(Vector3)

TMPL_XML_READ_WRITE_STREAM(ComplexVector)
TMPL_XML_READ_WRITE_STREAM(ComplexMatrix)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE_STREAM

#endif /* xml_io_general_types_h */
