add_library(species STATIC
  isotopologues.cc
)
target_link_libraries(species PUBLIC util arts_enum_options xml_io)
target_include_directories(species PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(species_tags STATIC species_tags.cc)
target_link_libraries(species_tags PUBLIC species xml_io)
target_include_directories(species_tags PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
