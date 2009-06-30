cmake_minimum_required (VERSION 2.6)

project (ARTS)
set (VERSION "1.12.44")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)

include (CheckIncludeFile)
include (CheckIncludeFileCXX)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

#include (MacroLibrary)

#include (ConvenienceLibs.cmake)

#include (ManualStuff.cmake)

include (CheckCXXCompilerFlag)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
  set (CMAKE_COMPILER_IS_INTELCXX 1)
endif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")

set (DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/arts)

if (NOT NUMERIC)
  set (NUMERIC double)
endif (NOT NUMERIC)

if (NOT INDEX)
  set (INDEX long)
endif (NOT INDEX)

add_definitions (-DHAVE_CONFIG_H)

if (NO_ASSERT)
  set (CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -DNDEBUG")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
endif (NO_ASSERT)

########### Check C header files ##########
check_include_file (stdlib.h    HAVE_STDLIB_H)
check_include_file (string.h    HAVE_STRING_H)
check_include_file (strings.h   HAVE_STRINGS_H)
check_include_file (sys/stat.h  HAVE_SYS_STAT_H)
check_include_file (sys/times.h HAVE_SYS_TIMES_H)
check_include_file (sys/types.h HAVE_SYS_TYPES_H)
check_include_file (unistd.h    HAVE_UNISTD_H)

########### Check C++ header files ##########
check_include_file_cxx (cstdlib HAVE_CSTDLIB)
check_include_file_cxx (cstring HAVE_CSTRING)
check_include_file_cxx (sstream HAVE_SSTREAM)
check_include_file_cxx (ctime   HAVE_CTIME)

########### Enable timer support, if we have all ##########
########### necessary header files               ##########
if (HAVE_SYS_TYPES_H AND HAVE_UNISTD_H AND HAVE_SYS_STAT_H AND HAVE_CTIME
    AND HAVE_SYS_TIMES_H)
  set (TIME_SUPPORT 1)
endif (HAVE_SYS_TYPES_H AND HAVE_UNISTD_H AND HAVE_SYS_STAT_H AND HAVE_CTIME
  AND HAVE_SYS_TIMES_H)


########### OpenMP support ##########
if (NOT NO_OPENMP)
  check_cxx_compiler_flag (-fopenmp CXX_OPENMP_SUPPORTED)
  if (CXX_OPENMP_SUPPORTED)
    if (CMAKE_COMPILER_IS_GNUCXX)
      set (OPENMP_FLAGS "-fopenmp")
    endif (CMAKE_COMPILER_IS_GNUCXX)

    if (CMAKE_COMPILER_IS_INTELCXX)
      set (OPENMP_FLAGS "-openmp")
    endif (CMAKE_COMPILER_IS_INTELCXX)

    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENMP_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_FLAGS}")
    set (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
    set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
  endif (CXX_OPENMP_SUPPORTED)

endif (NOT NO_OPENMP)

########### Special flags for certain compilers ##########
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-60")
endif (CMAKE_COMPILER_IS_GNUCXX)

########### Warning flags ##########
if (CMAKE_COMPILER_IS_INTELCXX)
  set (WARN_FLAGS "${WARN_FLAGS} -Wall -Wshadow -w1 -diag-disable 383,980,981,1418,1419")
endif (CMAKE_COMPILER_IS_INTELCXX)
if (CMAKE_COMPILER_IS_GNUCXX)
  set (WARN_FLAGS "${WARN_FLAGS} -W -Wall -Wconversion -Wshadow")
endif (CMAKE_COMPILER_IS_GNUCXX)

if (NOT NO_WERROR)
  set (WARN_FLAGS "${WARN_FLAGS} -Werror")
endif (NOT NO_WERROR)

if (NOT CXX_OPENMP_SUPPORTED)
  set (WARN_FLAGS "${WARN_FLAGS} -Wno-unknown-pragmas")
endif (NOT CXX_OPENMP_SUPPORTED)

string (TOUPPER ${CMAKE_BUILD_TYPE} MY_BUILD_TYPE)
set (COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${MY_BUILD_TYPE}}")

if (CMAKE_COMPILER_IS_INTELCXX)
  set (COMPILER "Intel C++")
else (CMAKE_COMPILER_IS_INTELCXX)
  if (CMAKE_COMPILER_IS_GNUCXX)
    set (COMPILER "GNU C++")
  else (CMAKE_COMPILER_IS_GNUCXX)
    set (COMPILER "unknown")
  endif (CMAKE_COMPILER_IS_GNUCXX)
endif (CMAKE_COMPILER_IS_INTELCXX)

set (COMPILER "${COMPILER} (${CMAKE_CXX_COMPILER})")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${WARN_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARN_FLAGS}")

########### Set Operating System name for config.h ##########
if (APPLE)
  set (OS_NAME "Mac OS")
elseif (UNIX)
  set (OS_NAME "Unix")
elseif (WIN)
  set (OS_NAME "Windows")
endif (APPLE)
  
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)

########### subdirectories ###############

add_subdirectory (src)
#add_subdirectory (aii)
#add_subdirectory (includes)
#add_subdirectory (doc)
#add_subdirectory (tests)
#add_subdirectory (3rdparty)

