# - Find the NetCDF library
#
# NETCDF_FOUND       - system has netcdf
# NETCDF_INCLUDE_DIR - the netcdf include directory
# NETCDF_LIBRARIES   - Link these to use netcdf
#
# Copyright (c) 2009, Oliver Lemke, <olemke@core-dump.info>

if (NOT NO_NETCDF)
  find_path (NETCDF_INCLUDE_DIR netcdf.h PATH_SUFFIXES netcdf)

  find_library (NETCDF_LIBRARY NAMES netcdf4 netcdf)

  set (NETCDF_LIBRARIES ${NETCDF_LIBRARY} ${NETCDFXX_LIBRARY})

  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (NetCDF DEFAULT_MSG
                                     NETCDF_LIBRARIES
                                     NETCDF_INCLUDE_DIR)

  mark_as_advanced (NETCDF_INCLUDE_DIR NETCDF_LIBRARY NETCDFXX_LIBRARY)

  if (NETCDF_FOUND)
    set (ENABLE_NETCDF 1)
  endif (NETCDF_FOUND)
endif (NOT NO_NETCDF)

