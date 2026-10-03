add_executable(rt4-test rt4-test.cpp)
target_link_libraries(rt4-test arts_rt4 physics)
add_test(NAME "cpp.fast.rt4-test" COMMAND rt4-test)
add_dependencies(check-deps rt4-test)
