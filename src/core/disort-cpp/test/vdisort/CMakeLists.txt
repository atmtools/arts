add_executable(vdisort-cpp-test vdisort-test.cpp)
target_link_libraries(vdisort-cpp-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-cpp-test" COMMAND vdisort-cpp-test)
add_dependencies(check-deps vdisort-cpp-test)
