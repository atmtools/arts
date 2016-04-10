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
add_executable(sphere sphere.cpp)
target_link_libraries(sphere matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
add_executable(linear linear.cpp)
target_link_libraries(linear matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
