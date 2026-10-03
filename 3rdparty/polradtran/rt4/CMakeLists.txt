# Evans' RT4 (polradtran) as modified for ARTS 2.6; the ARTS3 changes are
# listed in ../README.  rt4.f (the original main program) and the *.orig
# files are kept for provenance only and are not built.  The matrix
# routines (radmat.f) are the shared ../radmat.f (polradtran_radmat).
if(NOT ENABLE_RT4)
  return()
endif()

set(rt4_legacy_sources radtran4.f radutil4.f radintg4.f radscat4.f)
add_library(rt4 ${polradtran_library_type} ${rt4_legacy_sources} rt4_c_interface.f90)
set_target_properties(rt4 PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(rt4 PUBLIC polradtran_radmat)
set_source_files_properties(${rt4_legacy_sources} PROPERTIES COMPILE_OPTIONS "${polradtran_legacy_options}")
