add_executable(test_rtepack_perf test_rtepack_perf.cpp)

target_link_libraries(test_rtepack_perf PUBLIC rtepack rng)
