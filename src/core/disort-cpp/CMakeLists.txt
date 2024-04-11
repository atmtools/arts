add_library(disort-cpp STATIC disort.cpp)

target_link_libraries(disort-cpp PUBLIC matpack legendre artsworkspace_gpl_selector)
target_include_directories(disort-cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(disort-cpp-test disort-test.cpp)
target_link_libraries(disort-cpp-test disort-cpp)
