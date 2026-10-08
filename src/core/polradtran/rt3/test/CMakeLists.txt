add_executable(rt3-test rt3-test.cpp)
target_link_libraries(rt3-test arts_rt3 physics)
add_test(NAME "cpp.fast.rt3-test" COMMAND rt3-test)
add_dependencies(check-deps rt3-test)

add_executable(rt3-arts-test rt3-arts-test.cpp)
target_link_libraries(rt3-arts-test arts_rt3 physics legendre)
add_test(NAME "cpp.fast.rt3-arts-test" COMMAND rt3-arts-test)
add_dependencies(check-deps rt3-arts-test)

# The C++ RADTRAN against the Fortran one it ports
add_executable(rt3-radtran-test rt3-radtran-test.cpp)
target_link_libraries(rt3-radtran-test arts_rt3 rt3)
add_test(NAME "cpp.fast.rt3-radtran-test" COMMAND rt3-radtran-test)
add_dependencies(check-deps rt3-radtran-test)
