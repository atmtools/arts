add_library(lookup STATIC lookup_map.cpp)

target_link_libraries(lookup PUBLIC matpack arts_options lbl atm)
target_include_directories(lookup PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
