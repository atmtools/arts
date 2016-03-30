cmake_minimum_required (VERSION 2.6)

add_subdirectory(optimization)
add_subdirectory(algebra)
add_subdirectory(forward_models)
add_subdirectory(profiling)

include_directories(SYSTEM "${Boost_INCLUDE_DIRS}")
include_directories("${PROJECT_SOURCE_DIR}/src/"
                    "${PROJECT_SOURCE_DIR}/test/")

link_directories(${Boost_LIBRARIES})

add_executable(run_all run_all.cpp)
