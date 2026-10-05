add_executable(vdisort-cpp-test vdisort-test.cpp)
target_link_libraries(vdisort-cpp-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-cpp-test" COMMAND vdisort-cpp-test)
add_dependencies(check-deps vdisort-cpp-test)

add_executable(vdisort-conservative-test vdisort-test-conservative.cpp)
target_link_libraries(vdisort-conservative-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-conservative-test" COMMAND vdisort-conservative-test)
add_dependencies(check-deps vdisort-conservative-test)

add_executable(vdisort-arts-test vdisort-arts-test.cpp)
target_link_libraries(vdisort-arts-test vdisort-cpp physics legendre)
add_test(NAME "cpp.fast.vdisort-arts-test" COMMAND vdisort-arts-test)
add_dependencies(check-deps vdisort-arts-test)

add_executable(vdisort-fortran-test fortran-disort-test.cpp)
target_link_libraries(vdisort-fortran-test vdisort-cpp disort-cpp artstime)
add_test(NAME "cpp.fast.vdisort-fortran-test" COMMAND vdisort-fortran-test)
add_dependencies(check-deps vdisort-fortran-test)

# VDISORT against Evans' RT4 doubling-adding solver (src/core/rt4), an external reference
if (ENABLE_RT4)
  add_executable(vdisort-rt4-test vdisort-rt4-comparison.cpp)
  target_link_libraries(vdisort-rt4-test vdisort-cpp disort-cpp arts_rt4 physics artstime)
  add_test(NAME "cpp.fast.vdisort-rt4-test" COMMAND vdisort-rt4-test)
  add_dependencies(check-deps vdisort-rt4-test)
endif ()

# VDISORT against Evans' RT3 doubling-adding solver (src/core/rt3), an external reference
if (ENABLE_RT3)
  add_executable(vdisort-rt3-test vdisort-rt3-comparison.cpp)
  target_link_libraries(vdisort-rt3-test vdisort-cpp disort-cpp arts_rt3 physics legendre artstime)
  # E: Evans' benchmark scripts, read by src/tests/polradtran/evans-scripts.h
  target_include_directories(vdisort-rt3-test PRIVATE ${ARTS_SOURCE_DIR}/src/tests/polradtran)
  target_compile_definitions(vdisort-rt3-test PRIVATE POLRADTRAN_DIR="${ARTS_SOURCE_DIR}/3rdparty/polradtran")
  add_test(NAME "cpp.fast.vdisort-rt3-test" COMMAND vdisort-rt3-test)
  add_dependencies(check-deps vdisort-rt3-test)
endif ()

# RT4, RT3 and VDISORT on one ARTS atmosphere through their ARTS-native path builders
if (ENABLE_RT3 AND ENABLE_RT4)
  add_executable(vdisort-arts-comparison vdisort-arts-comparison.cpp)
  target_link_libraries(vdisort-arts-comparison vdisort-cpp arts_rt3 arts_rt4 physics legendre path)
  add_test(NAME "cpp.fast.vdisort-arts-comparison" COMMAND vdisort-arts-comparison)
  add_dependencies(check-deps vdisort-arts-comparison)
endif ()
