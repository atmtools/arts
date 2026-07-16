add_library(species STATIC
  isotopologues.cc isotopologues_class.cc auto_isotopologues.h auto_isotopologues.cc auto_species_enum_info.cc
)
target_link_libraries(species PUBLIC arts_enum_options util xml_io species_info species_enum_info)
target_include_directories(species PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(species_tags STATIC species_tags.cc)
target_link_libraries(species_tags PUBLIC species arts_enum_options xml_io)
target_include_directories(species_tags PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(make_species make_species.cc)
target_link_libraries(make_species PUBLIC species_info arts_enum_options xml_io)

add_library(species_info STATIC species_info.cc)
target_link_libraries(species_info PUBLIC xml_io arts_enum_options)
target_include_directories(species_info PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(species_enum_info STATIC species_enum_info.cc)
target_link_libraries(species_enum_info PUBLIC util xml_io)
target_include_directories(species_enum_info PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(make_species_enum make_species_enum.cc)
target_link_libraries(make_species_enum PUBLIC species_enum_info)
target_include_directories(make_species_enum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_command(
  OUTPUT auto_species_enum_info.cc auto_species_enum_info.h auto_isotopologues.cc auto_partition_functions.h
  COMMAND make_species_enum ${ARTS_CAT_DATA_DIR}
  DEPENDS make_species_enum
)

add_library(species_enum_info_gen STATIC auto_species_enum_info.cc)
target_link_libraries(species_enum_info_gen PUBLIC species_enum_info)
target_include_directories(species_enum_info_gen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(
  OUTPUT auto_isotopologues.h
  COMMAND make_species ${ARTS_CAT_DATA_DIR}
  DEPENDS make_species
)
