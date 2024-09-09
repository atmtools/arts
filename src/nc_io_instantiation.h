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

#include <workspace.h>

#include <stdexcept>

#include "nc_io.h"
#include "nc_io_types.h"

template <WorkspaceGroup T>
void nca_write_to_file(const String& filename, const T& type) {
  String efilename = add_basedir(filename);

  bool fail = false;
  String fail_msg;
#pragma omp critical(netcdf__critical_region)
  {
    int ncid;
    if (nc_create(efilename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid)) {
      fail = true;
      fail_msg = "Error opening file for writing.";
    } else {
      try {
        nca_write_to_file(ncid, type);
      } catch (const std::runtime_error &e) {
        fail = true;
        fail_msg = e.what();
      }
      nc_close(ncid);
    }
  }

  if (fail)
    ARTS_USER_ERROR("Error writing file: {}\n{}", efilename, fail_msg);
}

template <WorkspaceGroup T>
void nca_read_from_file(const String &filename, T &type) {
  String efilename = expand_path(filename);

  bool fail = false;
  String fail_msg;
#pragma omp critical(netcdf__critical_region)
  {
    int ncid;
    if (nc_open(efilename.c_str(), NC_NOWRITE, &ncid)) {
      fail = true;
      fail_msg = "Error opening file. Does it exists?";
    } else {
      try {
        nca_read_from_file(ncid, type);
      } catch (const std::runtime_error &e) {
        fail = true;
        fail_msg = e.what();
      }
      nc_close(ncid);
    }
  }

  if (fail)
    ARTS_USER_ERROR("Error reading file: {}\n{}", efilename, fail_msg);
}

#define TMPL_NC_READ_WRITE_FILE(what)                                \
  template void nca_write_to_file<what>(const String&, const what&); \
  template void nca_read_from_file<what>(const String&, what&);

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
xml_parse_from_stream (istream&, Vector&, bifstream *, ArtsXMLTag&);

void
xml_read_from_stream (istream&, ArrayOfLineRecord&,
                      const Numeric, const Numeric, bifstream * = NULL);

void
xml_parse_from_stream (istream&, ArrayOfString&, bifstream *, ArtsXMLTag&);*/

#endif /* nc_io_instantiation_h */
