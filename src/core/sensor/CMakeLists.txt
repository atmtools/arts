add_library(sensor STATIC obsel.cpp)
target_include_directories(sensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(sensor PUBLIC matpack rtepack)
