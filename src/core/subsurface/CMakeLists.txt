add_library(subsurface STATIC subsurf_field.cpp xml_subsurf.cpp)

target_link_libraries(subsurface PUBLIC arts_enum_options matpack operators)
target_include_directories(subsurface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
