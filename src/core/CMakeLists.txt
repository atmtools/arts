
# #####
add_subdirectory(absorption)
add_subdirectory(artstime)
add_subdirectory(atm)

add_subdirectory(binio)

add_subdirectory(coretypes)

add_subdirectory(disort-cpp)

add_subdirectory(functional)
add_subdirectory(fwd)

add_subdirectory(geodesy)

add_subdirectory(igrf)

add_subdirectory(jacobian)

add_subdirectory(lbl)
add_subdirectory(licenses)
add_subdirectory(legendre)
add_subdirectory(lookup)

add_subdirectory(matpack)
add_subdirectory(montecarlo)

add_subdirectory(operators)
add_subdirectory(options)

add_subdirectory(path)
add_subdirectory(physics)
add_subdirectory(predefined)

add_subdirectory(quantum)

add_subdirectory(rng)
add_subdirectory(rtepack)

add_subdirectory(scattering)
add_subdirectory(sensor)
add_subdirectory(spec)
add_subdirectory(strings)
add_subdirectory(sun)
add_subdirectory(surface)
add_subdirectory(subsurface)

add_subdirectory(polradtran)
add_subdirectory(tmatrix)

add_subdirectory(util)

add_subdirectory(tests)

add_subdirectory(xml)

# #####
add_library(artscore INTERFACE)

target_link_libraries(artscore INTERFACE
  arts_licenses
  absorption
  arts_enum_options
  artstime
  atm
  binio
  coretypes
  disort-cpp
  fwd
  jacobian
  lbl
  legendre
  lookup
  matpack
  montecarlo
  operators
  path
  physics
  predef
  quantum
  rng
  rtepack
  sensor
  species
  species_tags
  subsurface
  sun
  surface
  arts_rt3
  arts_rt4
  arts_tmatrix
  util
  xml_io
)

target_include_directories(artscore INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

if(ENABLE_PCH)
  target_precompile_headers(artscore PRIVATE <matpack.h>)
endif()

