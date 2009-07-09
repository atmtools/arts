# - Check for OpenMP support
#
# Copyright (c) 2009, Oliver Lemke, <olemke@core-dump.info>

if (NOT NO_OPENMP)
  check_cxx_compiler_flag (-fopenmp CXX_OPENMP_SUPPORTED)

  if (CXX_OPENMP_SUPPORTED)
    if (CMAKE_COMPILER_IS_GNUCXX)
      set (OPENMP_FLAGS "-fopenmp")
    endif (CMAKE_COMPILER_IS_GNUCXX)

    if (CMAKE_COMPILER_IS_INTELCXX)
      set (OPENMP_FLAGS "-openmp")
    endif (CMAKE_COMPILER_IS_INTELCXX)

    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENMP_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_FLAGS}")
    set (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
    set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
  endif (CXX_OPENMP_SUPPORTED)
endif (NOT NO_OPENMP)

