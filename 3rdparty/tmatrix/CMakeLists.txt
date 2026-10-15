# Original ARTS2 sources; wrapper corrections are listed in
# doc/arts/dev.tmatrix.rst.
if(NOT ENABLE_TMATRIX)
  message(STATUS "T-matrix is disabled")
  return()
else()
  message(STATUS "T-matrix is enabled")
endif()
if(ENABLE_ARTS_LGPL)
  message(FATAL_ERROR "T-matrix is unavailable in ENABLE_ARTS_LGPL builds")
endif()
# Mishchenko's licence permits use in not-for-profit research and has no
# SPDX identifier.
arts_add_bundled_license(tmatrix SPDX LicenseRef-Mishchenko-T-matrix FILES License.txt)
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/tmatrix_fortran.h
  MACRO_NAMESPACE "TMATRIX_FC_")
# Keep GNU's Darwin unwind settings confined to Fortran code. Applying
# -no_compact_unwind to the mixed C++ link breaks Clang exception unwinding.
set(tmatrix_library_type STATIC)
if(APPLE AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(tmatrix_library_type SHARED)
endif()
if(ENABLE_TMATRIX_QUAD)
  add_library(tmatrix ${tmatrix_library_type} tmq.arts.f amplq.arts.f lpq.f avg_tmatrixq.f)
else()
  add_library(tmatrix ${tmatrix_library_type} tmd.arts.f ampld.arts.f lpd.f avg_tmatrix.f)
endif()
set_target_properties(tmatrix PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(tmatrix INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(tmatrix PRIVATE -std=legacy -fallow-argument-mismatch)
  if(ENABLE_TMATRIX_QUAD)
    # GNU equivalents of the original compiler-specific quad intrinsics.
    target_compile_options(tmatrix PRIVATE -cpp -ffixed-line-length-none
      -Dqabs=abs -DQARCOS=ACOS -DQACOS=ACOS -DQABS=ABS -DQATAN=ATAN
      -DQCOS=COS -DQCOSH=COSH -DQSIN=SIN -DQSINH=SINH
      -DQSQRT=SQRT -DQMAX1=MAX -DQIMAG=AIMAG -DQCONJG=CONJG
      "-DQFLOAT(x)=REAL(x,16)" "-DQCMPLX(x,y)=CMPLX(x,y,16)")
  endif()
endif()
# Keep the independent original programs available for reference comparisons.
foreach(example IN ITEMS tmd ampld)
  add_executable(tmatrix_${example} EXCLUDE_FROM_ALL ${example}.lp.f lpd.f)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(tmatrix_${example} PRIVATE -std=legacy -fallow-argument-mismatch)
  endif()
endforeach()
