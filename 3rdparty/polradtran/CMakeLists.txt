# Evans' polradtran: RT4 (ENABLE_RT4) and RT3 (ENABLE_RT3), see README.ARTS.
# polradtran is MIT licensed (LICENSE), so both are also fine in LGPL builds.
# Built by default when a Fortran compiler is found (top-level CMakeLists.txt).
if(ARTS_FORTRAN_FOUND)
  set(polradtran_off_reason "turned off")
else()
  set(polradtran_off_reason "no Fortran compiler found; set FC or CMAKE_Fortran_COMPILER")
endif()
foreach(solver RT4 RT3)
  if(ENABLE_${solver})
    message(STATUS "${solver} is enabled (Fortran compiler ${CMAKE_Fortran_COMPILER})")
  else()
    message(STATUS "${solver} is disabled (${polradtran_off_reason})")
  endif()
endforeach()
if(NOT ENABLE_RT4 AND NOT ENABLE_RT3)
  return()
endif()
arts_add_bundled_license(polradtran SPDX MIT FILES LICENSE)

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
# boundary.
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
# MMULT calls BLAS DGEMM: the BLAS of ARTS's LAPACK (README.ARTS)
target_link_libraries(polradtran_radmat PRIVATE ${LAPACK_LIBRARIES})
set_source_files_properties(radmat.f PROPERTIES COMPILE_OPTIONS "${polradtran_radmat_options}")

# The libraries behind src/core/rt4 and src/core/rt3, from the sources with
# the ARTS changes listed in README.ARTS.
if(ENABLE_RT4)
  set(rt4_legacy_sources radtran4.f radutil4.f radintg4.f radscat4.f)
  add_library(rt4 ${polradtran_library_type} ${rt4_legacy_sources} rt4_c_interface.f90)
  set_target_properties(rt4 PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(rt4 PUBLIC polradtran_radmat)
  set_source_files_properties(${rt4_legacy_sources} PROPERTIES COMPILE_OPTIONS "${polradtran_legacy_options}")
endif()
if(ENABLE_RT3)
  set(rt3_legacy_sources radtran3.f radscat3.f radintg3.f radutil3.f)
  add_library(rt3 ${polradtran_library_type} ${rt3_legacy_sources} rt3_c_interface.f90)
  set_target_properties(rt3 PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(rt3 PUBLIC polradtran_radmat)
  set_source_files_properties(${rt3_legacy_sources} PROPERTIES COMPILE_OPTIONS "${polradtran_legacy_options}")
endif()

# Evans' original programs, from the tar's sources (the *.orig files copied
# to *.f so that CMake compiles them as Fortran, and the unchanged rt3.f,
# rt4.f and scatcnv.f).  They run his benchmark scripts runmietest, runtesta
# (rt3-evans), runtestr (scatcnv-evans and rt4-evans) and runtestc
# (rt4-evans); see src/tests/polradtran/evans-scripts.h.  They are not part of
# the libraries.
function(polradtran_original_program name)
  set(sources)
  foreach(source ${ARGN})
    if(source MATCHES "[.]orig$")
      string(REGEX REPLACE "[.]orig$" "" copy ${source})
      configure_file(${source} ${CMAKE_CURRENT_BINARY_DIR}/original/${copy} COPYONLY)
      list(APPEND sources ${CMAKE_CURRENT_BINARY_DIR}/original/${copy})
    else()
      list(APPEND sources ${CMAKE_CURRENT_SOURCE_DIR}/${source})
    endif()
  endforeach()
  add_executable(${name} ${sources})
  set_source_files_properties(${sources} PROPERTIES COMPILE_OPTIONS "${polradtran_legacy_options}")
endfunction()
if(ENABLE_RT3)
  polradtran_original_program(rt3-evans
    rt3.f radtran3.f.orig radscat3.f.orig radintg3.f.orig radutil3.f.orig radmat.f.orig)
endif()
if(ENABLE_RT4)
  polradtran_original_program(rt4-evans
    rt4.f radtran4.f.orig radscat4.f.orig radintg4.f.orig radutil4.f.orig radmat.f.orig)
  polradtran_original_program(scatcnv-evans scatcnv.f)
endif()
