add_executable(rt4-test rt4-test.cpp)
target_link_libraries(rt4-test arts_rt4 physics)
add_test(NAME "cpp.fast.rt4-test" COMMAND rt4-test)
add_dependencies(check-deps rt4-test)

add_executable(rt4-arts-test rt4-arts-test.cpp)
target_link_libraries(rt4-arts-test arts_rt4 physics)
add_test(NAME "cpp.fast.rt4-arts-test" COMMAND rt4-arts-test)
add_dependencies(check-deps rt4-arts-test)

# The C++ RADTRANO against the Fortran one it ports
add_executable(rt4-radtrano-test rt4-radtrano-test.cpp)
target_link_libraries(rt4-radtrano-test arts_rt4 rt4 physics)
add_test(NAME "cpp.fast.rt4-radtrano-test" COMMAND rt4-radtrano-test)
add_dependencies(check-deps rt4-radtrano-test)
