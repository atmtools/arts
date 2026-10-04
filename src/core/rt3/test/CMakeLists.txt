add_executable(rt3-test rt3-test.cpp)
target_link_libraries(rt3-test arts_rt3 physics)
add_test(NAME "cpp.fast.rt3-test" COMMAND rt3-test)
add_dependencies(check-deps rt3-test)

add_executable(rt3-arts-test rt3-arts-test.cpp)
target_link_libraries(rt3-arts-test arts_rt3 physics legendre)
add_test(NAME "cpp.fast.rt3-arts-test" COMMAND rt3-arts-test)
add_dependencies(check-deps rt3-arts-test)
