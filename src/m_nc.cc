#include <workspace.h>

#ifdef ENABLE_NETCDF

#include "nc_io.h"

#endif

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void ReaderNetCDF(  // WS Generic Input:
    T& v [[maybe_unused]],
    const String& f [[maybe_unused]]) {
  ARTS_TIME_REPORT

#ifdef ENABLE_NETCDF
  nca_read_from_file(f, v);
#else
  throw std::runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
#endif
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void WriterNetCDF(  // WS Generic Input:
    const T& v [[maybe_unused]],
    const String& f [[maybe_unused]]) {
  ARTS_TIME_REPORT

#ifdef ENABLE_NETCDF
  String filename = f;

  // Create default filename if empty
  nca_filename(filename);

  nca_write_to_file(filename, v);
#else
  throw std::runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
#endif
}

/* Workspace method: Doxygen documentation will be auto-generated */
template <WorkspaceGroup T>
void WriterNetCDFIndexed(  //WS Input:
    const Index& file_index [[maybe_unused]],
    // WS Generic Input:
    const T& v [[maybe_unused]],
    const String& f [[maybe_unused]]) {
  ARTS_TIME_REPORT

#ifdef ENABLE_NETCDF
  String filename = f;

  // Create default filename if empty
  nca_filename_with_index(filename, file_index);

  WriteNetCDF(v, filename);
#else
  throw std::runtime_error(
      "This version of ARTS was compiled without NetCDF support.");
#endif
}

#define ARTS_NETCDF_READER(T)                                              \
  void WriteNetCDF(const T& v, const String& f) { WriterNetCDF<T>(v, f); } \
  void WriteNetCDFIndexed(                                                 \
      const Index& file_index, const T& v, const String& f) {              \
    WriterNetCDFIndexed<T>(file_index, v, f);                              \
  }                                                                        \
  void ReadNetCDF(T& v, const String& f) { ReaderNetCDF<T>(v, f); }

ARTS_NETCDF_READER(ArrayOfIndex)
ARTS_NETCDF_READER(Vector)
ARTS_NETCDF_READER(Matrix)
ARTS_NETCDF_READER(Tensor3)
ARTS_NETCDF_READER(Tensor4)
ARTS_NETCDF_READER(Tensor5)
ARTS_NETCDF_READER(ArrayOfVector)
ARTS_NETCDF_READER(ArrayOfMatrix)
