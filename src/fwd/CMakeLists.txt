add_library(fwd STATIC
  fwd_abs.cpp
  fwd_cia.cpp
  fwd_lbl.cpp
  fwd_lbl_mtckd_voigt.cpp
  fwd_predef.cpp
  fwd_profile.cpp
)
target_link_libraries(fwd PUBLIC matpack Faddeeva predef)
target_include_directories(fwd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
