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

add_subdirectory(atm)
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

######

add_library(artscore STATIC
  absorption.cc
  absorptionlines.cc
  cia.cc
  check_input.cc
  cloudbox.cc
  covariance_matrix.cc
  fastem.cc
  gas_abs_lookup.cc
  gas_scattering.cc
  geodetic.cc
  hitran_species.cc
  igrf13.cc
  interpolation.cc
  jacobian.cc
  jpl_species.cc
  legendre.cc
  linemixing.cc
  linemixing_hitran.cc
  linescaling.cc
  lineshape.cc
  lineshapemodel.cc
  mc_antenna.cc
  mc_interp.cc
  microphysics.cc
  minimize.cc
  nlte.cc
  physics_funcs.cc
  optproperties.cc
  ppath_struct.cc
  predefined_absorption_models.cc
  radiation_field.cc
  raw.cc
  refraction.cc
  sensor.cc
  special_interp.cc
  sun.cc
  telsem.cc
  tessem.cc
  tmatrix.cc
  wigner_functions.cc
  xsec_fit.cc
  zeeman.cc
  zeemandata.cc
)

target_link_libraries(artscore PUBLIC
  artstime
  atm
  binio
  coretypes
  fwd
  lbl
  matpack
  operators
  path
  predef
  quantum
  rtepack
  species
  species_tags
  surface
  wigner
)

target_include_directories(artscore PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
