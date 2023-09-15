add_library(compoundtypes STATIC 
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
  linescaling.cc
  lineshapemodel.cc
  mc_antenna.cc
  nlte.cc
  optproperties.cc
  ppath_struct.cc
  sun.cc
  surf.cc
  telsem.cc
  tessem.cc
  xsec_fit.cc
  zeemandata.cc
)

target_link_libraries(compoundtypes PUBLIC
  artstime
  binio
  coretypes
  griddedfield
  matpack
  quantum
  rtepack
  species_tags
)

target_include_directories(compoundtypes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
