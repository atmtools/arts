add_library(path STATIC path_point.cpp atm_path.cpp)

target_link_libraries(path PUBLIC matpack atm surface)
target_include_directories(path PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
