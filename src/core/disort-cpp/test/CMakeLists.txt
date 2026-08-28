# Keep the large reference cases in separate translation units.
foreach (case IN ITEMS 1 2 3 4 5 8 9 11)
  add_executable(disort-cpp-test-${case} disort-test-${case}.cpp)
  target_link_libraries(disort-cpp-test-${case} disort-cpp artstime rng)
  add_test(NAME "cpp.fast.disort-cpp-test-${case}" COMMAND disort-cpp-test-${case})
  add_dependencies(check-deps disort-cpp-test-${case})
endforeach ()

add_executable(disort-cpp-disotest disotest.cpp)
target_link_libraries(disort-cpp-disotest disort-cpp artstime rng)
add_test(NAME "cpp.fast.disort-cpp-disotest" COMMAND disort-cpp-disotest)
add_dependencies(check-deps disort-cpp-disotest)

add_executable(vdisort-cpp-test vdisort-test.cpp)
target_link_libraries(vdisort-cpp-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-cpp-test" COMMAND vdisort-cpp-test)
add_dependencies(check-deps vdisort-cpp-test)

add_executable(disort-test-clearsky-multilayer disort-test-clearsky-multilayer.cpp)
target_link_libraries(disort-test-clearsky-multilayer disort-cpp artstime rng)
add_test(NAME "cpp.fast.disort-test-clearsky-multilayer" COMMAND disort-test-clearsky-multilayer)
add_dependencies(check-deps disort-test-clearsky-multilayer)

if (NOT ENABLE_ARTS_LGPL AND NOT CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  add_executable(test-old-impl test-old-impl.cpp)
  target_link_libraries(test-old-impl disort-cpp artstime cdisort rng)
endif ()
