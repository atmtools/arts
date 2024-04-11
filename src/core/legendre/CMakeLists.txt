add_library(legendre STATIC legendre.cpp)

target_link_libraries(legendre PUBLIC matpack)
target_include_directories(legendre PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
