/*!
  \file   m_xml.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   2002-06-18

  \brief  Workspace methods and template functions for supergeneric XML IO.

*/

#ifndef m_xml_h
#define m_xml_h

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include <workspace.h>
#include "exceptions.h"
#include "xml_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void ReadXML(  // WS Generic Output:
    T& v,
    // WS Generic Input:
    const String& f) {
  String filename = f;

  // Create default filename if empty
  filename_xml(filename);

  xml_read_from_file(filename, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void ReadXMLIndexed(  // WS Generic Output:
    T& v,
    // WS Input:
    const Index& file_index,
    // WS Generic Input:
    const String& f,
    const Index& digits) {
  String filename = f;

  // Create default filename if empty
  filename_xml_with_index(filename, file_index, digits);

  xml_read_from_file(filename, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void WriteXML(  //WS Input:
    const String& file_format,
    // WS Generic Input:
    const T& v,
    const String& f,
    const Index& no_clobber)

{
  // If MPI is enabled make sure only master process performs the write.
#ifdef ENABLE_MPI
  int initialized;
  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(nullptr, nullptr);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) {
    return;
  }
#endif  // ENABLE_MPI

  String filename = f;

  // Create default filename if empty
  filename_xml(filename);

  const FileType ftype = string2filetype(file_format);

  String errmsg;

#pragma omp critical(WriteXML_critical_region)
  {
    try {
      xml_write_to_file(filename, v, ftype, no_clobber);
    } catch (const std::exception& e) {
      errmsg = e.what();
    }
  }

  ARTS_USER_ERROR_IF (errmsg.length(), errmsg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void WriteXMLIndexed(  //WS Input:
    const String& file_format,
    const Index& file_index,
    // WS Generic Input:
    const T& v,
    const String& f,
    const Index& digits) {
  String filename = f;

  // Create default filename if empty
  filename_xml_with_index(filename, file_index, digits);

  WriteXML(file_format, v, filename, 0);
}

#endif  // m_xml_h
