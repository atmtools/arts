add_subdirectory(pythonic-disort)
add_subdirectory(reference-disort)
add_subdirectory(vdisort)

add_executable(disort-test-clearsky-multilayer disort-test-clearsky-multilayer.cpp)
target_link_libraries(disort-test-clearsky-multilayer disort-cpp artstime rng)
add_test(NAME "cpp.fast.disort-test-clearsky-multilayer" COMMAND disort-test-clearsky-multilayer)
add_dependencies(check-deps disort-test-clearsky-multilayer)
