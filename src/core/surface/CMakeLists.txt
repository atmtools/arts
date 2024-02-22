add_library(surface STATIC surf.cpp)

target_link_libraries(surface PUBLIC arts_enum_options matpack)
target_include_directories(surface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
