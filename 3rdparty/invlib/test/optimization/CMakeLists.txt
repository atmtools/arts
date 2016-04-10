cmake_minimum_required (VERSION 2.6)
project (maplib)

include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")
include_directories("${PROJECT_SOURCE_DIR}/src/"
                    "${PROJECT_SOURCE_DIR}/test/"
                    "/home/simon/src/arts/src"
                    "/home/simon/src/arts/build/src"
                    "/home/simon/src/arts/3rdparty")

link_directories(${Boost_LIBRARIES}
                "/home/simon/src/arts/build/src")

link_directories(${Boost_LIBRARIES})

# Executables
add_executable(exact exact.cpp)
target_link_libraries(exact matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
add_executable(test_functions test_functions.cpp)
target_link_libraries(test_functions matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
