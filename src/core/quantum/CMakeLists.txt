add_library(quantum STATIC
            quantum.cc
            quantum_term_symbol.cc
            hitran_species.cc
            jpl_species.cc
            auto_hitran_species_map.h
            auto_jpl_species_map.h
            hitran_species_info.cc
            jpl_species_info.cc
            )
target_link_libraries(quantum PUBLIC species arts_enum_options matpack binio)
target_include_directories(quantum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantum PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_executable(make_molparam_map make_molparam_map.cc hitran_species_info.cc)
target_link_libraries(make_molparam_map PUBLIC species)
target_include_directories(make_molparam_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

arts_generate(make_molparam_map
  ARGS ${ARTS_CAT_DATA_DIR}
  OUTPUT auto_hitran_species_map.h
  COMMENT "Generating automatic HITRAN species map"
)

add_executable(make_jpl_data make_jpl_data.cc jpl_species_info.cc)
target_link_libraries(make_jpl_data PUBLIC species)
target_include_directories(make_jpl_data PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

arts_generate(make_jpl_data
  ARGS ${ARTS_CAT_DATA_DIR}
  OUTPUT auto_jpl_species_map.h
  COMMENT "Generating automatic JPL species map"
)
