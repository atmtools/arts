# Evans' polradtran: RT4 (ENABLE_RT4) and RT3 (ENABLE_RT3), see README.
# polradtran is MIT licensed (LICENSE), so both are also fine in LGPL builds.
if(ENABLE_RT4)
  message(STATUS "RT4 is enabled")
else()
  message(STATUS "RT4 is disabled")
endif()
if(ENABLE_RT3)
  message(STATUS "RT3 is enabled")
else()
  message(STATUS "RT3 is disabled")
endif()
if(NOT ENABLE_RT4 AND NOT ENABLE_RT3)
  return()
endif()

# CMake initialises CMAKE_Fortran_FLAGS from the FFLAGS environment variable
# the first time Fortran is enabled in a build tree; conda environments export
# FFLAGS (e.g. -march=nocona -O2 ...).  Configure with -DCMAKE_Fortran_FLAGS=""
# for a reproducible polradtran build.  The flags in use are reported here so
# that a leak is visible.
message(STATUS "polradtran CMAKE_Fortran_FLAGS: '${CMAKE_Fortran_FLAGS}'")

# Keep GNU's Darwin unwind settings confined to Fortran code (see tmatrix).
set(polradtran_library_type STATIC)
if(APPLE AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(polradtran_library_type SHARED)
endif()

# The legacy sources use default REAL literals and variables (e.g. LOG(2.)
# for the doubling count); ARTS 2.6 compiled them in double precision, so
# promote default REAL to 8 bytes there.  Default INTEGER is not promoted:
# the rt*_c_interface.f90 shims use explicit kinds and convert at the C
# boundary.  set_source_files_properties is directory scoped, so rt4/ and
# rt3/ apply these options to their own sources.
set(polradtran_legacy_options)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(polradtran_legacy_options -std=legacy -fdefault-real-8 -fdefault-double-8)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  set(polradtran_legacy_options -r8)
endif()

# radmat.f: Evans' general matrix routines (MCOPY ... MINVERT and LINPACK's
# DGEFA/DGEDI), used by both RT4 and RT3.  They keep no state between calls
# (no COMMON, SAVE or DATA), so one copy serves both solvers.  RT4 and RT3
# serialise their calls with separate mutexes, so a call into radmat from
# one may run concurrently with a call from the other: MINVERT's local
# arrays INDX and WORK must therefore live on the stack.  GNU does that by
# default for arrays of this size; -frecursive (GNU) and -auto (Intel, whose
# default is static local arrays) make it explicit.
set(polradtran_radmat_options ${polradtran_legacy_options})
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  list(APPEND polradtran_radmat_options -frecursive)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
  list(APPEND polradtran_radmat_options -auto)
endif()
add_library(polradtran_radmat ${polradtran_library_type} radmat.f)
set_target_properties(polradtran_radmat PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_source_files_properties(radmat.f PROPERTIES COMPILE_OPTIONS "${polradtran_radmat_options}")

add_subdirectory(rt4)
add_subdirectory(rt3)
