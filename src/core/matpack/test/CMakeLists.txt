add_executable(test_matpack_perf test_matpack_perf.cpp)
target_link_libraries(test_matpack_perf PUBLIC matpack time_report rng)

add_executable(test_complex_linalg test_complex_linalg.cpp)
target_link_libraries(test_complex_linalg PRIVATE matpack)
add_test(NAME "cpp.fast.test_complex_linalg" COMMAND test_complex_linalg)
add_dependencies(check-deps test_complex_linalg)
