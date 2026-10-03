# Evans' RT3 (polradtran) with the ARTS3 changes listed in README.  rt3.f
# (the original main program), runmietest and runtesta (Evans' benchmark
# scripts, which hold the expected outputs used by src/core/rt3/test) and
# the *.orig files are kept for provenance only and are not built.  The
# matrix routines (radmat.f) are the shared ../radmat.f (polradtran_radmat).
if(NOT ENABLE_RT3)
  return()
endif()

set(rt3_legacy_sources radtran3.f radscat3.f radintg3.f radutil3.f)
add_library(rt3 ${polradtran_library_type} ${rt3_legacy_sources} rt3_c_interface.f90)
set_target_properties(rt3 PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(rt3 PUBLIC polradtran_radmat)
set_source_files_properties(${rt3_legacy_sources} PROPERTIES COMPILE_OPTIONS "${polradtran_legacy_options}")
