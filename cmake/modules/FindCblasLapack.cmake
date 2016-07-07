# - Find the Lapack library
#
# LAPACK_FOUND       - system has lapack
# LAPACK_LIBRARIES   - Link these to use lapack
#
# Copyright (c) 2016, Oliver Lemke, <olemke@core-dump.info>

find_library (CBLAS_LIBRARIES NAMES cblas)
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Cblas DEFAULT_MSG
                                   CBLAS_LIBRARIES)

if (CBLAS_FOUND)
  unset(LAPACK_LIBRARIES)
  find_library (LAPACK_LIBRARIES NAMES lapack)
  find_package_handle_standard_args (Lapack DEFAULT_MSG
                                     LAPACK_LIBRARIES)
  set(LAPACK_LIBRARIES ${CBLAS_LIBRARIES} ${LAPACK_LIBRARIES})
 
  if (LAPACK_FOUND)
    set (ENABLE_LAPACK)
  endif ()
endif ()

