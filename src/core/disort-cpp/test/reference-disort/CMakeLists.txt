# Cases derived from cDISORT and the original Fortran DISORT suite.
add_executable(disort-cpp-disotest disotest.cpp)
target_link_libraries(disort-cpp-disotest disort-cpp artstime rng)
add_test(NAME "cpp.fast.disort-cpp-disotest" COMMAND disort-cpp-disotest)
add_dependencies(check-deps disort-cpp-disotest)

add_executable(fortran-disort-test fortran-disort-test.cpp)
target_link_libraries(fortran-disort-test disort-cpp artstime)
add_test(NAME "cpp.fast.fortran-disort-test" COMMAND fortran-disort-test)
add_dependencies(check-deps fortran-disort-test)

if (NOT ENABLE_ARTS_LGPL AND NOT CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  add_executable(test-old-impl test-old-impl.cpp)
  target_link_libraries(test-old-impl disort-cpp artstime cdisort rng)
endif ()
