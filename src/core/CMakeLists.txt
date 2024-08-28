# #### Give Index, Numeric, String, and Array<> types
set(CORETYPE_EXTRA)

add_library(coretypes STATIC
  ${CORETYPE_EXTRA}
  file.cc
  mystring.cpp
  parameters.cc
)
target_include_directories(coretypes PUBLIC "${ARTS_SOURCE_DIR}/3rdparty")
target_include_directories(coretypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(coretypes PUBLIC "${ARTS_BINARY_DIR}/src")
target_link_libraries(coretypes PUBLIC util)

# #####
add_subdirectory(absorption)
add_subdirectory(atm)
add_subdirectory(physics)
add_subdirectory(binio)
add_subdirectory(matpack)
add_subdirectory(rtepack)
add_subdirectory(predefined)
add_subdirectory(operators)
add_subdirectory(fwd)
add_subdirectory(lbl)
add_subdirectory(quantum)
add_subdirectory(spec)
add_subdirectory(artstime)
add_subdirectory(surface)
add_subdirectory(path)
add_subdirectory(util)
add_subdirectory(options)
add_subdirectory(sensor)
add_subdirectory(disort-cpp)
add_subdirectory(legendre)

# #####
add_library(artscore STATIC
  absorptionlines.cc
  check_input.cc
  cloudbox.cc
  covariance_matrix.cc
  fastem.cc
  gas_abs_lookup.cc
  gas_scattering.cc
  geodetic.cc
  igrf13.cc
  interpolation.cc
  jacobian.cc
  jpl_species.cc
  lineshapemodel.cc
  mc_antenna.cc
  mc_interp.cc
  microphysics.cc
  minimize.cc
  nlte.cc
  optproperties.cc
  radiation_field.cc
  raw.cc
  refraction.cc
  retrieval_target.cc
  sensor.cc
  special_interp.cc
  sun.cc
  telsem.cc
  tessem.cc
  tmatrix.cc
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
  util
  disort-cpp
  legendre
)

target_include_directories(artscore PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
