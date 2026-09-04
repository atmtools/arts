add_subdirectory(pythonic-disort)
add_subdirectory(reference-disort)
add_subdirectory(vdisort)

add_executable(disort-common-test common-test.cpp)
target_link_libraries(disort-common-test disort-cpp-common)
add_test(NAME "cpp.fast.disort-common-test" COMMAND disort-common-test)
add_dependencies(check-deps disort-common-test)

add_executable(disort-test-delta-m-plus-settings disort-test-delta-m-plus-settings.cpp)
target_link_libraries(disort-test-delta-m-plus-settings disort-cpp)
add_test(NAME "cpp.fast.disort-test-delta-m-plus-settings" COMMAND disort-test-delta-m-plus-settings)
add_dependencies(check-deps disort-test-delta-m-plus-settings)

add_executable(disort-test-clearsky-multilayer disort-test-clearsky-multilayer.cpp)
target_link_libraries(disort-test-clearsky-multilayer disort-cpp artstime rng)
add_test(NAME "cpp.fast.disort-test-clearsky-multilayer" COMMAND disort-test-clearsky-multilayer)
add_dependencies(check-deps disort-test-clearsky-multilayer)

add_executable(disort-test-conservative disort-test-conservative.cpp)
target_link_libraries(disort-test-conservative disort-cpp artstime)
add_test(NAME "cpp.fast.disort-test-conservative" COMMAND disort-test-conservative)
add_dependencies(check-deps disort-test-conservative)
