# #### Give Index, Numeric, String, and Array<> types
set(CORETYPE_EXTRA)

add_library(coretypes STATIC
  ${CORETYPE_EXTRA}
  file.cc
  parameters.cc
)
target_include_directories(coretypes PUBLIC "${ARTS_SOURCE_DIR}/3rdparty")
target_include_directories(coretypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(coretypes PUBLIC "${ARTS_BINARY_DIR}/src")
target_link_libraries(coretypes PUBLIC util strings)

# #####
add_subdirectory(absorption)
add_subdirectory(atm)
add_subdirectory(physics)
add_subdirectory(binio)
add_subdirectory(matpack)
add_subdirectory(rtepack)
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
add_subdirectory(surface)
add_subdirectory(subsurface)
add_subdirectory(path)
add_subdirectory(util)
add_subdirectory(options)
add_subdirectory(sensor)
add_subdirectory(disort-cpp)
add_subdirectory(legendre)
add_subdirectory(strings)

# #####
add_library(artscore STATIC
  covariance_matrix.cc
  gas_scattering.cc
  geodetic.cc
  igrf13.cc
  interpolation.cc
  minimize.cc
  optproperties.cc
  sun.cc
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
)

target_include_directories(artscore PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

if(ENABLE_PCH)
  target_precompile_headers(artscore PRIVATE <matpack.h>)
endif()

