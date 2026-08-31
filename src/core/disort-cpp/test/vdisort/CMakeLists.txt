add_executable(vdisort-cpp-test vdisort-test.cpp)
target_link_libraries(vdisort-cpp-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-cpp-test" COMMAND vdisort-cpp-test)
add_dependencies(check-deps vdisort-cpp-test)

add_executable(vdisort-conservative-test vdisort-test-conservative.cpp)
target_link_libraries(vdisort-conservative-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-conservative-test" COMMAND vdisort-conservative-test)
add_dependencies(check-deps vdisort-conservative-test)

add_executable(vdisort-fortran-test fortran-disort-test.cpp)
target_link_libraries(vdisort-fortran-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-fortran-test" COMMAND vdisort-fortran-test)
add_dependencies(check-deps vdisort-fortran-test)
