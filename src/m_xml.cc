#include <workspace.h>
#include <xml.h>
#ifdef ENABLE_MPI
#include "mpi.h"
#endif

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadXML(  // WS Generic Output:
    AnyOutput v,
    // WS Generic Input:
    const String& f) {
  ARTS_TIME_REPORT

  String filename = f;

  // Create default filename if empty
  filename_xml(filename);

  std::visit([&](const auto& ptr) { xml_read_from_file(filename, *ptr); }, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadXMLIndexed(  // WS Generic Output:
    AnyOutput v,
    // WS Input:
    const Index& file_index,
    // WS Generic Input:
    const String& f,
    const Index&  digits) {
  ARTS_TIME_REPORT

  String filename = f;

  // Create default filename if empty
  filename_xml_with_index(filename, file_index, digits);

  std::visit([&](const auto& ptr) { xml_read_from_file(filename, *ptr); }, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void WriteXML(  //WS Input:
    const String& file_format,
    // WS Generic Input:
    const AnyInput v,
    const String&  f,
    const Index&   no_clobber)

{
  ARTS_TIME_REPORT

  // If MPI is enabled make sure only master process performs the write.
#ifdef ENABLE_MPI
  int initialized;
  MPI_Initialized(&initialized);
  if (!initialized) { MPI_Init(nullptr, nullptr); }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) { return; }
#endif  // ENABLE_MPI

  String filename = f;

  // Create default filename if empty
  filename_xml(filename);

  const FileType ftype = to<FileType>(file_format);

  String errmsg;

#pragma omp critical(WriteXML_critical_region)
  {
    try {
      std::visit([&](const auto& ptr) { xml_write_to_file(filename, *ptr, ftype, no_clobber); }, v);
    } catch (const std::exception& e) { errmsg = e.what(); }
  }

  ARTS_USER_ERROR_IF(errmsg.length(), "{}", errmsg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void WriteXMLIndexed(  //WS Input:
    const String& file_format,
    const Index&  file_index,
    // WS Generic Input:
    const AnyInput v,
    const String&  f,
    const Index&   digits) {
  ARTS_TIME_REPORT

  String filename = f;

  // Create default filename if empty
  filename_xml_with_index(filename, file_index, digits);

  WriteXML(file_format, v, filename, 0);
}
