# - Find the NetCDF library
#
# NETCDF_FOUND       - system has netcdf
# NETCDF_INCLUDE_DIR - the netcdf include directory
# NETCDF_LIBRARIES   - Link these to use netcdf
#
# Copyright (c) 2018, Oliver Lemke, <olemke@core-dump.info>

if (NOT NO_FFTW)
  find_path (FFTW_INCLUDE_DIR fftw3.h)

  find_library (FFTW_LIBRARY NAMES fftw3)

  set (FFTW_LIBRARIES ${FFTW_LIBRARY})

  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (
    FFTW DEFAULT_MSG
    FFTW_LIBRARIES
    FFTW_INCLUDE_DIR
    )

  mark_as_advanced (FFTW_INCLUDE_DIR FFTW_LIBRARY)

  if (FFTW_FOUND)
    set (ENABLE_FFTW 1)
  endif (FFTW_FOUND)
endif (NOT NO_FFTW)

