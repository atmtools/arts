add_executable(rt3-test rt3-test.cpp)
target_link_libraries(rt3-test arts_rt3 physics)
add_test(NAME "cpp.fast.rt3-test" COMMAND rt3-test)
add_dependencies(check-deps rt3-test)
