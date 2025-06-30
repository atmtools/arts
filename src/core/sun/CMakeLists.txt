add_library(sun STATIC sun.cc)
target_link_libraries(sun PUBLIC rtepack surface path)
target_include_directories(sun PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})