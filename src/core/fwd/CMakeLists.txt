add_library(fwd STATIC
  fwd_cia.cpp
  fwd_hxsec.cpp
  fwd_predef.cpp
  fwd_propmat.cpp
  fwd_spectral_radiance.cpp
)
target_link_libraries(fwd PUBLIC predef lbl)
target_include_directories(fwd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
