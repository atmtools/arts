cmake_minimum_required (VERSION 2.6)

add_subdirectory(optimization)
add_subdirectory(algebra)
add_subdirectory(forward_models)
add_subdirectory(profiling)
add_subdirectory(mpi)

include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")
include_directories("${PROJECT_SOURCE_DIR}/src/"
                    "${PROJECT_SOURCE_DIR}/test/"
                    "/home/simon/src/arts/src"
                    "/home/simon/src/arts/build/src"
                    "/home/simon/src/arts/3rdparty")

link_directories(${Boost_LIBRARIES}
                "/home/simon/src/arts/build/src")

add_executable(run_all run_all.cpp)
target_link_libraries(run_all matpack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
