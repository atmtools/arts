# The T-matrix of a sphere against the independent Mie solver of the scattering module
add_executable(tmatrix-sphere-test tmatrix-sphere-test.cpp)
target_link_libraries(tmatrix-sphere-test arts_tmatrix physics)
target_include_directories(tmatrix-sphere-test PRIVATE ${PROJECT_SOURCE_DIR}/src/core/scattering)
add_test(NAME "cpp.fast.tmatrix-sphere-test" COMMAND tmatrix-sphere-test)
add_dependencies(check-deps tmatrix-sphere-test)

# The numerical maximum dimension of the shapes against the closed forms, and the Chebyshev radius ratio
add_executable(tmatrix-shape-test tmatrix-shape-test.cpp)
target_link_libraries(tmatrix-shape-test arts_tmatrix legendre)
add_test(NAME "cpp.fast.tmatrix-shape-test" COMMAND tmatrix-shape-test)
add_dependencies(check-deps tmatrix-shape-test)

# The amplitude matrix towards the particle axis, the laboratory poles and beta = 90 degrees
add_executable(tmatrix-singular-test tmatrix-singular-test.cpp)
target_link_libraries(tmatrix-singular-test arts_tmatrix)
add_test(NAME "cpp.fast.tmatrix-singular-test" COMMAND tmatrix-singular-test)
add_dependencies(check-deps tmatrix-singular-test)

# Azimuthally random orientation against Mie, random orientation, and its own exact quadratures
add_executable(tmatrix-aro-test tmatrix-aro-test.cpp)
target_link_libraries(tmatrix-aro-test arts_tmatrix scattering physics legendre)
target_include_directories(tmatrix-aro-test PRIVATE ${PROJECT_SOURCE_DIR}/src/core/scattering)
add_test(NAME "cpp.fast.tmatrix-aro-test" COMMAND tmatrix-aro-test)
add_dependencies(check-deps tmatrix-aro-test)

# The C++ port against recorded and printed outputs of the original Fortran
add_executable(tmatrix-reference-test tmatrix-reference-test.cpp)
target_link_libraries(tmatrix-reference-test arts_tmatrix physics)
add_test(NAME "cpp.fast.tmatrix-reference-test"
         COMMAND tmatrix-reference-test ${PROJECT_SOURCE_DIR}/3rdparty/tmatrix)
add_dependencies(check-deps tmatrix-reference-test)
