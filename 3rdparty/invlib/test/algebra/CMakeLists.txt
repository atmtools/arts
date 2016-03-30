cmake_minimum_required (VERSION 2.6)
project (invlib)

include_directories (SYSTEM "${Boost_INCLUDE_DIRS}")
include_directories ("${PROJECT_SOURCE_DIR}/src/"
                     "${PROJECT_SOURCE_DIR}/test/")

link_directories(${Boost_LIBRARIES})

# Executables
add_executable(identities identities.cpp)
add_executable(solvers solvers.cpp)
add_executable(precision_matrix precision_matrix.cpp)
add_executable(transformation transformation.cpp)
