cmake_minimum_required (VERSION 2.6)
project (invlib)

include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")
include_directories("${PROJECT_SOURCE_DIR}/src/"
                    "${PROJECT_SOURCE_DIR}/test/"
                    "/home/simon/src/arts/src"
                    "/home/simon/src/arts/build/src"
                    "/home/simon/src/arts/3rdparty")

link_directories(${Boost_LIBRARIES}
                "/home/simon/src/arts/build/src")

# Executables
add_executable(identities identities.cpp)
target_link_libraries(identities matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
add_executable(solvers solvers.cpp)
target_link_libraries(solvers matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
add_executable(precision_matrix precision_matrix.cpp)
target_link_libraries(precision_matrix matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
add_executable(transformation transformation.cpp)
target_link_libraries(transformation matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
