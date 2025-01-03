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

#ifndef nc_io_types_h
#define nc_io_types_h

#include "array.h"
#include "matpack_arrays.h"

#define TMPL_NC_READ_WRITE_FILE(what)             \
  void nca_write_to_file(const int, const what&); \
  void nca_read_from_file(const int, what&);

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

//=== Array Types ==========================================================

TMPL_NC_READ_WRITE_FILE(ArrayOfIndex)
TMPL_NC_READ_WRITE_FILE(ArrayOfMatrix)
TMPL_NC_READ_WRITE_FILE(ArrayOfVector)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_NC_READ_WRITE_FILE

#endif /* nc_io_types_h */
