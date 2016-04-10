cmake_minimum_required (VERSION 2.6)
project (invlib)

# Output build info.
message(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )

# C++11 support.
set (CMAKE_CXX_STANDARD 11)

# Check for BOOST library.
include(FindBLAS)
include(FindLAPACK)
find_package(Boost 1.60.0)

include_directories ("${PROJECT_SOURCE_DIR}/src/"
                     "${PROJECT_SOURCE_DIR}/test/"
                     "${Boost_INCLUDE_DIRS}")

# Build tests.
add_subdirectory(test)

