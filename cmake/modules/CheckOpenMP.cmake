# - Check for OpenMP support
#
# Copyright (c) 2010, Oliver Lemke, <olemke@core-dump.info>

if (NOT NO_OPENMP)
  check_cxx_compiler_flag (-fopenmp OPENMP_FOUND)

  if (OPENMP_FOUND)
    if (CMAKE_COMPILER_IS_GNUCXX)
      set (OPENMP_FLAGS "-fopenmp")
    endif (CMAKE_COMPILER_IS_GNUCXX)

    if (CMAKE_COMPILER_IS_INTELCXX)
      set (OPENMP_FLAGS "-openmp")
    endif (CMAKE_COMPILER_IS_INTELCXX)

   mark_as_advanced(
     OPENMP_FOUND
     OPENMP_FLAGS
   )

  endif (OPENMP_FOUND)
endif (NOT NO_OPENMP)

