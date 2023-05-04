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

#ifndef nc_io_instantiation_h
#define nc_io_instantiation_h

#include <cfloat>
#include <stdexcept>
#include "nc_io.h"
#include "nc_io_types.h"

#define TMPL_NC_READ_WRITE_FILE(what)                \
  template void nca_write_to_file<what>(             \
      const String&, const what&, const Verbosity&); \
  template void nca_read_from_file<what>(            \
      const String&, what&, const Verbosity&);

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for NetCDF streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_NC_READ_WRITE_FILE(Matrix)
TMPL_NC_READ_WRITE_FILE(Tensor3)
TMPL_NC_READ_WRITE_FILE(Tensor4)
TMPL_NC_READ_WRITE_FILE(Tensor5)
TMPL_NC_READ_WRITE_FILE(Vector)

//=== Compound Types =======================================================

TMPL_NC_READ_WRITE_FILE(Agenda)
TMPL_NC_READ_WRITE_FILE(GasAbsLookup)

//=== Array Types ==========================================================

TMPL_NC_READ_WRITE_FILE(ArrayOfIndex)
TMPL_NC_READ_WRITE_FILE(ArrayOfMatrix)
TMPL_NC_READ_WRITE_FILE(ArrayOfVector)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE

/*void
xml_parse_from_stream (istream&, Vector&, bifstream *, ArtsXMLTag&, const Verbosity&);

void
xml_read_from_stream (istream&, ArrayOfLineRecord&,
                      const Numeric, const Numeric, bifstream * = NULL,
                      const Verbosity&);

void
xml_parse_from_stream (istream&, ArrayOfString&, bifstream *, ArtsXMLTag&, const Verbosity&);*/

#endif /* nc_io_instantiation_h */
