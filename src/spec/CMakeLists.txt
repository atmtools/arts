add_library(species STATIC
            isotopologues.cc
            )
target_link_libraries(species PUBLIC auto_partfun)
target_include_directories(species PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(species_tags STATIC species_tags.cc)
target_link_libraries(species_tags PUBLIC partfun)
