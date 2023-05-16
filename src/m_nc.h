/*!
  \file   m_nc.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-09-19

  \brief  Workspace methods and template functions for supergeneric NetCDF IO.

*/

#ifndef m_nc_h
#define m_nc_h

#include "arts.h"

#include "agenda_class.h"
#include "exceptions.h"
#include "workspace_ng.h"

#ifdef ENABLE_NETCDF

#include "nc_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadNetCDF(  // WS Generic Input:
    T& v,
    const String& v_name _U_,
    const String& f,
    // WS Generic Input Names:
    const String& f_name _U_,
    const Verbosity& verbosity)

{
  nca_read_from_file(f, v, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDF(  // WS Generic Input:
    const T& v,
    const String& f,
    // WS Generic Input Names:
    const String& v_name,
    const String& f_name _U_,
    const Verbosity& verbosity)

{
  String filename = f;

  // Create default filename if empty
  nca_filename(filename, v_name);

  nca_write_to_file(filename, v, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDFIndexed(  //WS Input:
    const Index& file_index,
    // WS Generic Input:
    const T& v,
    const String& f,
    // WS Generic Input Names:
    const String& v_name,
    const String& f_name,
    const Verbosity& verbosity)

{
  String filename = f;

  // Create default filename if empty
  nca_filename_with_index(filename, file_index, v_name);

  WriteNetCDF(v, filename, v_name, f_name, verbosity);
}

#else  // NetCDF not enabled

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadNetCDF(  // WS Generic Input:
    T&,
    const String&,
    const String&,
    // WS Generic Input Names:
    const String&,
    const Verbosity&)

{
  throw runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDF(  // WS Generic Input:
    const T&,
    const String&,
    // WS Generic Input Names:
    const String&,
    const String&,
    const Verbosity&)

{
  throw runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void WriteNetCDFIndexed(  //WS Input:
    const Index&,
    // WS Generic Input:
    const T&,
    const String&,
    // WS Generic Input Names:
    const String&,
    const String&,
    const Verbosity&)

{
  throw runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
}

#endif  // ENABLE_NETCDF

/* Workspace method: Doxygen documentation will be auto-generated */
template <typename T>
void ReadNetCDF(Workspace& ws _U_,
                // WS Generic Input:
                T& v,
                const String& v_name,
                const String& f,
                // WS Generic Input Names:
                const String& f_name,
                const Verbosity& verbosity)

{
  ReadNetCDF(v, f, v_name, f_name, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void WriteNetCDF(Workspace& ws _U_,
                 // WS Generic Input:
                 const Agenda& v,
                 const String& f,
                 // WS Generic Input Names:
                 const String& v_name,
                 const String& f_name,
                 const Verbosity& verbosity) {
  WriteNetCDF(v, f, v_name, f_name, verbosity);
}

#endif  // m_nc_h
