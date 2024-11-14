////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io_arts_types.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_arts_types_h
#define xml_io_arts_types_h

#include <cfloat>
#include <stdexcept>

// All workspace groups are known by the agenda class through tokval.h
#include <workspace.h>

// Extras
#include "fwd_spectral_radiance.h"
#include "jacobian.h"
#include "matpack_data.h"
#include "mc_interp.h"
#include "operators.h"
#include "path_point.h"
#include "template_partfun.h"
#include "xml_io_general_types.h"

#define TMPL_XML_READ_WRITE_STREAM(what)                          \
  void xml_read_from_stream(std::istream &, what &, bifstream *); \
  void xml_write_to_stream(                                       \
      std::ostream &, const what &, bofstream *, const String &);

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_XML_READ_WRITE_STREAM(PartitionFunctionsData)
TMPL_XML_READ_WRITE_STREAM(ScatteringSpecies)

//=== Extras ===============================================================

TMPL_XML_READ_WRITE_STREAM(AtmFunctionalData)
TMPL_XML_READ_WRITE_STREAM(GridPos)

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE_STREAM

#endif /* xml_io_arts_types_h */
