
add_library(montecarlo STATIC mc_antenna.cc)

target_link_libraries(montecarlo PUBLIC
  arts_enum_options
  matpack
  rng
  rtepack
)

target_include_directories(montecarlo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
