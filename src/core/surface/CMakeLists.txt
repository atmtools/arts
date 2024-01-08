add_library(surface STATIC surf.cpp)

target_link_libraries(surface PUBLIC matpack griddedfield)
target_include_directories(surface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
