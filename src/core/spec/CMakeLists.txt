add_library(species STATIC
            isotopologues.cc
            species.cc
            )
target_link_libraries(species PUBLIC matpack)
target_include_directories(species PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(species_tags STATIC species_tags.cc)
target_link_libraries(species_tags PUBLIC species)
target_include_directories(species_tags PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
