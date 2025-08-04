
# #####
add_subdirectory(absorption)
add_subdirectory(atm)
add_subdirectory(physics)
add_subdirectory(binio)
add_subdirectory(matpack)
add_subdirectory(rtepack)
add_subdirectory(rng)
add_subdirectory(predefined)
add_subdirectory(operators)
add_subdirectory(jacobian)
add_subdirectory(fwd)
add_subdirectory(lbl)
add_subdirectory(lookup)
add_subdirectory(quantum)
add_subdirectory(scattering)
add_subdirectory(spec)
add_subdirectory(artstime)
add_subdirectory(geodesy)
add_subdirectory(surface)
add_subdirectory(subsurface)
add_subdirectory(path)
add_subdirectory(util)
add_subdirectory(options)
add_subdirectory(sensor)
add_subdirectory(disort-cpp)
add_subdirectory(legendre)
add_subdirectory(strings)
add_subdirectory(xml)
add_subdirectory(coretypes)
add_subdirectory(tests)
add_subdirectory(sun)
add_subdirectory(functional)
add_subdirectory(igrf)

# #####
add_library(artscore STATIC
  gas_scattering.cc
  spectral_radiance_transform_operator.cc
)

target_link_libraries(artscore PUBLIC
  arts_enum_options
  absorption
  artstime
  atm
  binio
  coretypes
  fwd
  lbl
  lookup
  matpack
  operators
  path
  physics
  predef
  quantum
  rtepack
  sensor
  species
  species_tags
  surface
  subsurface
  util
  disort-cpp
  legendre
  jacobian
  xml_io
  rng
  sun
)

target_include_directories(artscore PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

if(ENABLE_PCH)
  target_precompile_headers(artscore PRIVATE <matpack.h>)
endif()

