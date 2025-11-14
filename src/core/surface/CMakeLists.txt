add_library(surface STATIC surf_field.cpp xml_surf.cpp)

target_link_libraries(surface PUBLIC arts_enum_options matpack operators)
target_include_directories(surface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
