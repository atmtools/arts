add_library(fwd STATIC
  fwd_cia.cpp
  fwd_hxsec.cpp
  fwd_path.cpp
  fwd_predef.cpp
  fwd_propmat.cpp
  fwd_spectral_radiance.cpp
)
target_link_libraries(fwd PUBLIC path absorption)
target_include_directories(fwd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
