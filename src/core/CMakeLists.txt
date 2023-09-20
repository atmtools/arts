##### Give Index, Numeric, String, and Array<> types
set (CORETYPE_EXTRA)
if (NOT HAVE_GETOPT_H)
  list (APPEND CORETYPE_EXTRA
    arts_getopt.c arts_getopt1.c)
endif (NOT HAVE_GETOPT_H)

add_library(coretypes STATIC
  ${CORETYPE_EXTRA}
  array.cpp
  file.cc
  mystring.cpp
  parameters.cc
)
target_include_directories(coretypes PUBLIC "${ARTS_SOURCE_DIR}/3rdparty")
target_include_directories(coretypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(coretypes PUBLIC "${ARTS_BINARY_DIR}/src")

######

add_subdirectory(binio)
add_subdirectory(matpack)
add_subdirectory(rtepack)
add_subdirectory(predefined)
add_subdirectory(operators)
add_subdirectory(fwd)
add_subdirectory(quantum)
add_subdirectory(spec)
add_subdirectory(griddedfield)
add_subdirectory(artstime)

######

add_library(artstypes STATIC
  absorptionlines.cc
  atm.cc
  cia.cc
  check_input.cc
  covariance_matrix.cc
  gas_abs_lookup.cc
  hitran_species.cc
  interpolation.cc
  jacobian.cc
  jpl_species.cc
  linemixing.cc
  linemixing_hitran.cc
  linescaling.cc
  lineshape.cc
  lineshapemodel.cc
  mc_antenna.cc
  minimize.cc
  nlte.cc
  physics_funcs.cc
  optproperties.cc
  ppath_struct.cc
  predefined_absorption_models.cc
  sun.cc
  surf.cc
  telsem.cc
  tessem.cc
  xsec_fit.cc
  zeemandata.cc
)

target_link_libraries(artstypes PUBLIC
  artstime
  binio
  coretypes
  griddedfield
  fwd
  matpack
  operators
  predef
  quantum
  rtepack
  species
  species_tags
)

######

add_library(artsfun STATIC
  absorption.cc
  atm_path.cc
  cloudbox.cc
  fastem.cc
  gas_scattering.cc
  geodetic.cc
  igrf13.cc
  legendre.cc
  mc_interp.cc
  microphysics.cc
  radiation_field.cc
  raw.cc
  refraction.cc
  sensor.cc
  special_interp.cc
  tmatrix.cc
  wigner_functions.cc
  zeeman.cc
)

target_link_libraries(artsfun PUBLIC
  artstypes
  wigner
)
