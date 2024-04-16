add_library(disort-cpp STATIC disort.cpp)

target_link_libraries(disort-cpp PUBLIC matpack legendre artsworkspace_gpl_selector)
target_include_directories(disort-cpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add tests, split them because they contain too much data otherwise
add_executable(disort-cpp-test-1 disort-test-1.cpp)
add_executable(disort-cpp-test-2 disort-test-2.cpp)
add_executable(disort-cpp-test-3 disort-test-3.cpp)
add_executable(disort-cpp-test-4 disort-test-4.cpp)

target_link_libraries(disort-cpp-test-1 disort-cpp)
target_link_libraries(disort-cpp-test-2 disort-cpp)
target_link_libraries(disort-cpp-test-3 disort-cpp)
target_link_libraries(disort-cpp-test-4 disort-cpp)

add_dependencies(check-deps disort-cpp-test-1)
add_dependencies(check-deps disort-cpp-test-2)
add_dependencies(check-deps disort-cpp-test-3)
add_dependencies(check-deps disort-cpp-test-4)