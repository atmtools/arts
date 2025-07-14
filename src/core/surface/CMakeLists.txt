add_library(surface STATIC surface_field.cpp xml_surface.cpp)

target_link_libraries(surface PUBLIC arts_enum_options matpack operators)
target_include_directories(surface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
