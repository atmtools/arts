add_library(surface STATIC surf_field.cpp telsem.cpp tessem.cpp xml_surf.cpp)

target_link_libraries(surface PUBLIC arts_enum_options matpack operators functional_numeric_ternary)
target_include_directories(surface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
