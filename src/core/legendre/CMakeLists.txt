add_library(legendre STATIC legendre.cpp fastgl.cpp)

target_link_libraries(legendre PUBLIC matpack Boost::headers)
target_include_directories(legendre PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
