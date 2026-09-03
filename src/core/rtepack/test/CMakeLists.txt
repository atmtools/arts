add_executable(test_rtepack_perf test_rtepack_perf.cpp)

target_link_libraries(test_rtepack_perf PUBLIC rtepack rng)

add_executable(test_rtepack_magnus test_rtepack_magnus.cpp)
target_link_libraries(test_rtepack_magnus PUBLIC rtepack)
add_test(NAME "cpp.fast.test_rtepack_magnus" COMMAND test_rtepack_magnus)
add_dependencies(check-deps test_rtepack_magnus)

add_executable(test_rtepack_rtestep test_rtepack_rtestep.cpp)
target_link_libraries(test_rtepack_rtestep PUBLIC rtepack)
add_test(NAME "cpp.fast.test_rtepack_rtestep" COMMAND test_rtepack_rtestep)
add_dependencies(check-deps test_rtepack_rtestep)
