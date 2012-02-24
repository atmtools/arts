if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message (FATAL_ERROR
"In-source builds are not allowed. Please create a directory and run cmake "
"from there, passing the path to this source directory as the last argument:"
"
"
"cmake path_to_arts"
"
"
"IMPORTANT: This process created the file 'CMakeCache.txt' and the directory "
"'CMakeFiles'.  Please delete them."
"
"
"If you get this error message although you are in a separate build directory, "
"please check that there is no file 'CMakeCache.txt' and no directory "
"'CMakeFiles' in the arts source directory (not the directory you're "
"currently in).  Please delete them."
)
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

cmake_minimum_required (VERSION 2.6.3)
set (CMAKE_LEGACY_CYGWIN_WIN32 0)

project (ARTS)

enable_testing ()

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)

include (CheckIncludeFile)
include (CheckIncludeFileCXX)

include (CheckCXXCompilerFlag)
include (CheckTypeSize)
include (CheckFunctionExists)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

include (ArtsTestcases)

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
check_include_file (getopt.h    HAVE_GETOPT_H)

########### Check C header files for libmicrohttpd ##########
check_include_file (fcntl.h     HAVE_FCNTL_H)
check_include_file (math.h      HAVE_MATH_H)
check_include_file (errno.h     HAVE_ERRNO_H)
check_include_file (limits.h    HAVE_LIMITS_H)
check_include_file (stdio.h     HAVE_STDIO_H)
check_include_file (locale.h    HAVE_LOCALE_H)
check_include_file (pthread.h   HAVE_PTHREAD_H)

########### Check additional C header files for libmicrohttpd ##########
check_include_file (sys/time.h    HAVE_SYS_TIME_H)
check_include_file (sys/msg.h     HAVE_SYS_MSG_H)
check_include_file (netdb.h       HAVE_NETDB_H)
check_include_file (netinet/in.h  HAVE_NETINET_IN_H)
check_include_file (netinet/tcp.h HAVE_NETINET_TCP_H)
check_include_file (sys/socket.h  HAVE_SYS_SOCKET_H)
check_include_file (sys/mman.h    HAVE_SYS_MMAN_H)
check_include_file (arpa/inet.h   HAVE_ARPA_INET_H)
check_include_file (sys/select.h  HAVE_SYS_SELECT_H)
check_include_file (poll.h        HAVE_POLL_H)

########### Check C++ header files ##########
check_include_file_cxx (cstdlib HAVE_CSTDLIB)
check_include_file_cxx (cstring HAVE_CSTRING)
check_include_file_cxx (sstream HAVE_SSTREAM)
check_include_file_cxx (ctime   HAVE_CTIME)

########### Check type sizes ##########
check_type_size("int"    SIZEOF_INT)
check_type_size("long"   SIZEOF_LONG)
check_type_size("size_t" SIZEOF_SIZE_T)
check_type_size("float"  SIZEOF_FLOAT)
check_type_size("double" SIZEOF_DOUBLE)

########### Check functions ##########
check_function_exists(unlink HAVE_REMOVE)

########### Enable timer support, if we have all ##########
########### necessary header files               ##########
if (HAVE_SYS_TYPES_H AND HAVE_UNISTD_H AND HAVE_SYS_STAT_H AND HAVE_CTIME
    AND HAVE_SYS_TIMES_H)
  set (TIME_SUPPORT 1)
endif (HAVE_SYS_TYPES_H AND HAVE_UNISTD_H AND HAVE_SYS_STAT_H AND HAVE_CTIME
  AND HAVE_SYS_TIMES_H)

########### Zlib support ##########
include (FindZLIB)
set (ENABLE_ZLIB ${ZLIB_FOUND})

########### NetCDF support ##########
include (FindNetCDF)

########### OpenMP support ##########
if (NOT NO_OPENMP AND NOT CMAKE_GENERATOR STREQUAL "Xcode")
  include (FindOpenMP)
endif (NOT NO_OPENMP AND NOT CMAKE_GENERATOR STREQUAL "Xcode")

########### Docserver support ##########

if (NOT NO_DOCSERVER)
  set (ENABLE_DOCSERVER true)
endif (NOT NO_DOCSERVER)

########### Disort support ##########
if (ENABLE_DISORT)
  string (TOUPPER ${INDEX} INDEX_UPPER)
  if (NOT SIZEOF_${INDEX_UPPER} EQUAL 4)
    message (STATUS "Forcing Index type to int for compatibility with the f2c library")
    set (INDEX int)
  endif (NOT SIZEOF_${INDEX_UPPER} EQUAL 4)

  check_include_file (f2c.h HAVE_F2C_H)
  if (NOT HAVE_F2C_H)
    message (FATAL_ERROR "cannot find f2c.h which is required for disort")
  endif (NOT HAVE_F2C_H)

  set (NO_WERROR 1)
endif (ENABLE_DISORT)

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
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (WARN_FLAGS "${WARN_FLAGS} -W -Wall -Wconversion -Wno-sign-conversion -Wshadow")
endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

if (WERROR)
  set (WARN_FLAGS "${WARN_FLAGS} -Werror")
endif (WERROR)

if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
  set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
else (OPENMP_FOUND)
  set (WARN_FLAGS "${WARN_FLAGS} -Wno-unknown-pragmas")
endif (OPENMP_FOUND)

# Workaround for f2c library in Debian/Ubuntu
if (ENABLE_DISORT AND CMAKE_COMPILER_IS_GNUCXX AND UNIX)
  set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -u MAIN__")
endif (ENABLE_DISORT AND CMAKE_COMPILER_IS_GNUCXX AND UNIX)

string (TOUPPER ${CMAKE_BUILD_TYPE} MY_BUILD_TYPE)
set (COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${MY_BUILD_TYPE}}")

if (CMAKE_COMPILER_IS_INTELCXX)
  set (COMPILER "Intel C++")
else (CMAKE_COMPILER_IS_INTELCXX)
  if (CMAKE_COMPILER_IS_GNUCXX)
    set (COMPILER "GNU C++")
  else (CMAKE_COMPILER_IS_GNUCXX)
    # OLE NOTE: cmake 2.8.1 and earlier identify Clang as GNUCXX!!!
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set (CMAKE_COMPILER_IS_CLANG 1)
      set (COMPILER "LLVM/Clang")
    else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set (COMPILER "unknown")
    endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  endif (CMAKE_COMPILER_IS_GNUCXX)
endif (CMAKE_COMPILER_IS_INTELCXX)

set (COMPILER "${COMPILER} (${CMAKE_CXX_COMPILER})")

if (CMAKE_GENERATOR STREQUAL "Xcode")
  set (COMPILER "Xcode")
endif (CMAKE_GENERATOR STREQUAL "Xcode")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${WARN_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARN_FLAGS}")

########### Set Operating System name for config.h ##########
if (APPLE)
  set (OS_NAME "Mac OS")
  set (OSX 1)
elseif ("${CMAKE_SYSTEM_NAME}" STREQUAL "CYGWIN")
  set (OS_NAME "Cygwin")
  set (CYGWIN 1)
elseif (UNIX)
  set (OS_NAME "Unix")
  set (LINUX 1)
elseif (WIN)
  set (OS_NAME "Windows")
  set (WINDOWS 1)
endif (APPLE)

set (ARTS_DEFAULT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/includes")

if (WITH_GUI)
  find_package (Qt4)
  if (QT4_FOUND)
    set (ENABLE_GUI true)
  endif (QT4_FOUND)
endif (WITH_GUI)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)

add_custom_command (
  OUTPUT auto_version.txt
  COMMAND ${CMAKE_COMMAND} -D "ARTS_SOURCE_DIR:STRING=\"${CMAKE_SOURCE_DIR}\"" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/update_auto_version_txt.cmake
  DEPENDS ChangeLog
)

add_custom_target (UpdateAutoVersion ALL DEPENDS auto_version.txt)

########### subdirectories ###############

add_subdirectory (src)
#add_subdirectory (aii)
add_subdirectory (includes)
add_subdirectory (doc)
add_subdirectory (tests)
#add_subdirectory (3rdparty)

if (NOT NUMERIC STREQUAL "double")
  message (WARNING "You chose ${NUMERIC} as the Numeric type. Using anything different than double is completely unsupported, unrecommended and will probably lead to compilation errors, crashes and wrong results. Sorry, you're on your own here. Godspeed!")
endif (NOT NUMERIC STREQUAL "double")

