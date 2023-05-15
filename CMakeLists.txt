if (ENABLE_CXX23)
  cmake_minimum_required (VERSION 3.20)
else ()
  cmake_minimum_required (VERSION 3.18)
endif ()

project (ARTS C CXX)

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" OR "${ARTS_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
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
endif ()

######## Check compiler versions ########
set (REQUIRED_GCC_VERSION 11)
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${REQUIRED_GCC_VERSION})
  message (FATAL_ERROR "GCC version must be at least ${REQUIRED_GCC_VERSION}!")
endif()

set (REQUIRED_CLANG_VERSION 13)
if (CMAKE_CXX_COMPILER_ID MATCHES Clang AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${REQUIRED_CLANG_VERSION})
  message (FATAL_ERROR "Clang version must be at least ${REQUIRED_CLANG_VERSION}!")
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(DEFAULT_BUILD_TYPE "Release")
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "Native")
endif()

# User options (auto-complete and cmake gui works better like this)
include(CMakeDependentOption)
option (ENABLE_ARTS_LGPL "Make the code LGPL compatible" OFF)
option (ENABLE_CCACHE "Turn on CCACHE" ON)
option (ENABLE_CXX23 "Turn on C++23 (Don't use, only for future updates)" OFF)
option (ENABLE_FORTRAN "Turn on Fortran code" OFF)
option (ENABLE_GUI "Turn on Debug GUI" OFF)
option (ENABLE_IPO "Turn on IPO (experimental)" OFF)
option (ENABLE_MPI "Turn on MPI" OFF)
option (ENABLE_NETCDF "Turn on NETCDF support" OFF)
option (ENABLE_PCH "Turn on precomnpiled headers" OFF)
option (IGNORE_CONDA_PATH "Turn on pathing guard for conda environment" OFF)
option (NO_ASSERT "Turn off all asserts" OFF)
option (NO_DOCSERVER "Turn off the document server" OFF)
cmake_dependent_option(
    NO_FASTEM "Turn off FASTEM even if Fortran is enabled" OFF "ENABLE_FORTRAN" ON)
cmake_dependent_option(
    NO_RT4 "Turn off RT4 even if Fortran is enabled" OFF "ENABLE_FORTRAN" ON)
cmake_dependent_option(
    NO_TMATRIX "Turn off TMATRIX even if Fortran is enabled" OFF "ENABLE_FORTRAN" ON)
option (NO_USER_ERRORS "Turn off all user error checking (reckless)" OFF)
option (RELEASE_O2 "Use O2 instead of O3 in release mode" OFF)

if (IGNORE_CONDA_PATH)
  # Avoid possible linker incompatibilities for
  # libstdc++ and libgomp on Linux if anaconda is installed
  find_program(CONDA_PROG NAMES conda)
  if (CONDA_PROG)
    get_filename_component(CONDA_PATH ${CONDA_PROG} DIRECTORY)
    get_filename_component(CONDA_PATH ${CONDA_PATH} DIRECTORY)
    if (NOT CONDA_PATH MATCHES "/usr.*")
      set(CMAKE_IGNORE_PATH ${CMAKE_IGNORE_PATH}
              ${CONDA_PATH}/lib
              ${CONDA_PATH}/include)
    endif()
  endif()
endif()

########### License settings ##########
if (ENABLE_ARTS_LGPL)
  add_definitions(-DARTS_LGPL=1)
else()
  add_definitions(-DARTS_LGPL=0)
endif()

if (ENABLE_FORTRAN)
 if (NOT NO_FASTEM AND ENABLE_ARTS_LGPL)
  message(FATAL_ERROR "FASTEM is not available as LGPL.")
 elseif(NOT NO_RT4 AND ENABLE_ARTS_LGPL)
  message(FATAL_ERROR "RT4 is not available as LGPL.")
 elseif(NOT NO_TMATRIX AND ENABLE_ARTS_LGPL)
  message(FATAL_ERROR "TMATRIX is not available as LGPL.")
 endif()
endif()

if(ENABLE_ARTS_LGPL)
  message(STATUS "Compiling ARTS under LGPL")
endif()

########### C++20/23 Support ##########

if (NOT CMAKE_CXX_STANDARD)
  if (ENABLE_CXX23)
    set(CMAKE_CXX_STANDARD 23)
  else ()
    set(CMAKE_CXX_STANDARD 20)
  endif ()
endif ()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing ()

list (APPEND CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

########### Find Python ##########
find_package (Python3 REQUIRED COMPONENTS Interpreter Development.Module)

include (CheckPythonModule)
check_python_modules()

set(ARTS_PYTHON_INTERPRETER "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/python" ${Python3_EXECUTABLE})

if (ENABLE_FORTRAN)
  enable_language (Fortran)

  get_filename_component (FORTRAN_COMPILER ${CMAKE_Fortran_COMPILER} NAME)

  if (FORTRAN_COMPILER MATCHES "gfortran.*")
    if(CMAKE_BUILD_TYPE STREQUAL "Native")
      set (FORTRAN_EXTRA_FLAGS "-fdefault-double-8 -fdefault-integer-8 -fdefault-real-8 -march=native -O3")
    else ()
      set (FORTRAN_EXTRA_FLAGS "-fdefault-double-8 -fdefault-integer-8 -fdefault-real-8")
    endif()
    if (ENABLE_TMATRIX_QUAD)
      message (FATAL_ERROR "Gfortran does not support quad-precision.")
    endif()
  elseif (FORTRAN_COMPILER MATCHES "ifort.*")
    set (FORTRAN_EXTRA_FLAGS "-i8 -r8")
  else ()
    message (FATAL_ERROR "Unsupported Fortran compiler. Must be gfortran or ifort.")
  endif (FORTRAN_COMPILER MATCHES "gfortran.*")

endif (ENABLE_FORTRAN)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)

include (CheckIncludeFile)
include (CheckIncludeFileCXX)

include (CheckTypeSize)
include (CheckFunctionExists)

include (ArtsTestcases)

set (DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/arts)

if (NOT NUMERIC)
  set (NUMERIC double)
endif (NOT NUMERIC)

if (NOT INDEX)
  set (INDEX std::int64_t)
endif (NOT INDEX)

add_definitions (-DHAVE_CONFIG_H)
add_definitions (-DEIGEN_DONT_PARALLELIZE)

if (NO_ASSERT)
  add_definitions (-DNDEBUG)
endif (NO_ASSERT)

if (NO_USER_ERRORS)
  message(STATUS "Experimental reckless build enabled (-DNO_USER_ERRORS=0 to disable)")
  add_definitions (-DNO_ARTS_USER_ERRORS)
endif (NO_USER_ERRORS)

########### ccache Support ##########

if (ENABLE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if (CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif (CCACHE_FOUND)
endif ()

########### IPO ##########

if (ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT IPO_SUPPORTED LANGUAGES C CXX OUTPUT IPO_ERROR)
  if (NOT IPO_SUPPORTED)
    message(FATAL_ERROR "IPO / LTO not supported: ${IPO_ERROR}")
  endif()
else()
  set(IPO_SUPPORTED OFF)
endif()

########### OEM Support ##########

if (NOT NO_OEM)
  set (OEM_SUPPORT true)
endif ()

########### Check MPI Support ############
if (ENABLE_MPI) # User must enable MPI using -DENABLE_MPI=ON
  find_package(MPI)
  if (NOT MPI_CXX_FOUND)
    message (FATAL_ERROR "No MPI with C++ support found.")
  endif ()
endif ()

########### Python ##########

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
find_package( Threads REQUIRED )
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

########### Zlib support ##########
include (FindZLIB)
set (ENABLE_ZLIB ${ZLIB_FOUND} REQUIRED)

########## BLAS/LAPACK support #########
find_package (LAPACK REQUIRED)

########## NetCDF support #########
if (ENABLE_NETCDF)
  include (FindNetCDF)
  if (NOT NETCDF_FOUND)
      message(FATAL_ERROR "NetCDF library not found")
  endif()
endif (ENABLE_NETCDF)

########## OpenMP support #########
if (NOT NO_OPENMP AND NOT CMAKE_GENERATOR STREQUAL "Xcode")
  include (FindOpenMP)
endif (NOT NO_OPENMP AND NOT CMAKE_GENERATOR STREQUAL "Xcode")

########### Docserver support ##########

if (NOT NO_DOCSERVER)
  find_package (MHD)
  if (MHD_FOUND)
    set (ENABLE_DOCSERVER true)
    message (STATUS "libmicrohttpd found. Enabling docserver.")
  else ()
    message (STATUS "No libmicrohttpd found. Disabling docserver.")
  endif ()
endif ()

########### Remove NDEBUG as set by cmake >2.8.10 ##########
if (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  string (REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELWITHDEBINFO
          "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  string (REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO
          "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")

if (CMAKE_BUILD_TYPE STREQUAL "Release" AND RELEASE_O2)
  string(REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif ()

########### Special flags for certain compilers ##########
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-1024")
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

########### Warning flags ##########
include (ArtsAddCompilerFlag)
if (CMAKE_GENERATOR STREQUAL "Xcode")
  set (CMAKE_XCODE_ATTRIBUTE_GCC_WARN_SHADOW "YES")
  set (CMAKE_XCODE_ATTRIBUTE_GCC_WARN_UNUSED_PARAMETER "YES")
  set (CMAKE_XCODE_ATTRIBUTE_GCC_WARN_UNUSED_VARIABLE "YES")
  set (CMAKE_XCODE_ATTRIBUTE_CLANG_WARN_SUSPICIOUS_IMPLICIT_CONVERSION "YES")
else (CMAKE_GENERATOR STREQUAL "Xcode")
  ARTS_ADD_COMPILER_FLAG (W)
  ARTS_ADD_COMPILER_FLAG (Wall)
  ARTS_ADD_COMPILER_FLAG (Wextra)
  ARTS_ADD_COMPILER_FLAG (Wshadow)
  ARTS_ADD_COMPILER_FLAG (Wconversion)
  ARTS_ADD_COMPILER_FLAG (Wdeprecated)
  ARTS_ADD_COMPILER_FLAG (Wdeprecated-declarations)
  ARTS_ADD_COMPILER_FLAG (Wno-sign-conversion)
  ARTS_ADD_COMPILER_FLAG (Wno-unknown-pragmas)
  ARTS_ADD_COMPILER_FLAG (Wno-return-type-c-linkage)
  ARTS_ADD_COMPILER_FLAG (Wno-strict-overflow)
  ARTS_ADD_COMPILER_FLAG (Wno-psabi)
endif (CMAKE_GENERATOR STREQUAL "Xcode")

if(CMAKE_BUILD_TYPE STREQUAL "Native")
    ARTS_ADD_COMPILER_FLAG (O3)
    ARTS_ADD_COMPILER_FLAG (march=native)
    ARTS_ADD_COMPILER_FLAG (fcx-fortran-rules)
    ARTS_ADD_COMPILER_FLAG (DNDEBUG)
endif ()

if (NO_USER_ERRORS)
  ARTS_ADD_COMPILER_FLAG (Wno-unused-parameter)
  ARTS_ADD_COMPILER_FLAG (Wno-unused-variable)
endif (NO_USER_ERRORS)

if (WERROR)
  ARTS_ADD_COMPILER_FLAG (Werror)
endif (WERROR)

if (OpenMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

set (COMPILER "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER})")

########### Xcode specific settings ##########
if (CMAKE_GENERATOR STREQUAL "Xcode")
  set (COMPILER "Xcode")
  set (CMAKE_OSX_SYSROOT "")
  set (LIBRARY_OUTPUT_PATH ${ARTS_BINARY_DIR}/build CACHE PATH
       "Single output directory for building all libraries.")
  set (EXECUTABLE_OUTPUT_PATH ${ARTS_BINARY_DIR}/build CACHE PATH
       "Single output directory for building all executables.")
  mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

  ######### Enable C++11 and link against the new libc++ ##########
  if (WITH_XCODE_LIBCPP)
    set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set (CMAKE_XCODE_ATTRIBUTE_OBJROOT ${ARTS_BINARY_DIR}/build/\${CONFIGURATION})
    set (CMAKE_XCODE_ATTRIBUTE_DEBUG_INFORMATION_FORMAT "dwarf-with-dsym")
    set (CMAKE_XCODE_ATTRIBUTE_GCC_VERSION "com.apple.compilers.llvm.clang.1_0")
  endif (WITH_XCODE_LIBCPP)
endif (CMAKE_GENERATOR STREQUAL "Xcode")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ARTS_C_FLAGS} ${ARTS_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ARTS_CXX_FLAGS} ${ARTS_FLAGS}")


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

set (ARTS_DEFAULT_INCLUDE_DIR "${ARTS_SOURCE_DIR}/controlfiles")

string (TOUPPER ${CMAKE_BUILD_TYPE} MY_BUILD_TYPE)
set (COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${MY_BUILD_TYPE}}")

add_custom_target (
  auto_version.txt
  COMMAND ${CMAKE_COMMAND} -D "ARTS_SOURCE_DIR:STRING=\"${ARTS_SOURCE_DIR}\"" -P ${ARTS_SOURCE_DIR}/cmake/scripts/update_auto_version_txt.cmake
)

add_custom_target (UpdateAutoVersion ALL DEPENDS auto_version.txt)


########### Tests ##########
if (TEST_JOBS)
  set (CTEST_JOBS "-j${TEST_JOBS}")
endif (TEST_JOBS)


########### Path to arts-xml-data ##########
if (ARTS_XML_DATA_PATH)
  get_filename_component (ARTS_XML_DATA_DIR "${ARTS_XML_DATA_PATH}" ABSOLUTE)
else (ARTS_XML_DATA_PATH)
  get_filename_component (ARTS_XML_DATA_DIR "${CMAKE_SOURCE_DIR}/../arts-xml-data" ABSOLUTE)
endif (ARTS_XML_DATA_PATH)

if (NOT IS_DIRECTORY "${ARTS_XML_DATA_DIR}")
  message(WARNING "Cannot determine location of arts-xml-data!
Depending tests will be disabled.
Please specify arts-xml-data location with
cmake -DARTS_XML_DATA_PATH=PATH_TO_ARTS_XML_DATA ..")
  unset (ARTS_XML_DATA_DIR)
else (NOT IS_DIRECTORY "${ARTS_XML_DATA_DIR}")
  message(STATUS "Found arts-xml-data in ${ARTS_XML_DATA_DIR}")
endif (NOT IS_DIRECTORY "${ARTS_XML_DATA_DIR}")

########### Path to arts-cat-data ##########
if (ARTS_CAT_DATA_PATH)
  get_filename_component (ARTS_CAT_DATA_DIR "${ARTS_CAT_DATA_PATH}" ABSOLUTE)
else (ARTS_CAT_DATA_PATH)
  get_filename_component (ARTS_CAT_DATA_DIR "${CMAKE_SOURCE_DIR}/../arts-cat-data" ABSOLUTE)
endif (ARTS_CAT_DATA_PATH)

if (NOT IS_DIRECTORY "${ARTS_CAT_DATA_DIR}")
  message(WARNING "Cannot determine location of arts-cat-data!
Depending tests will be disabled.
Please specify arts-cat-data location with
cmake -DARTS_CAT_DATA_PATH=PATH_TO_ARTS_CAT_DATA ..")
  unset (ARTS_CAT_DATA_DIR)
else (NOT IS_DIRECTORY "${ARTS_CAT_DATA_DIR}")
  message(STATUS "Found arts-cat-data in ${ARTS_CAT_DATA_DIR}")
endif (NOT IS_DIRECTORY "${ARTS_CAT_DATA_DIR}")

string(JOIN ":" ARTS_DATA_DIR ${ARTS_XML_DATA_DIR} ${ARTS_CAT_DATA_DIR})

########### Path to partition functions data ##########
if (ARTS_PARTFUN_DATA_PATH)
  get_filename_component (ARTS_PARTFUN_DATA_PATH "${ARTS_PARTFUN_DATA_PATH}" ABSOLUTE)
    message(STATUS "Using manual partition functions data (unset ARTS_PARTFUN_DATA_PATH to switch to built-in options)")
else (ARTS_PARTFUN_DATA_PATH)
  get_filename_component (ARTS_PARTFUN_DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/partfun/data" ABSOLUTE)
  message(STATUS "Using built-in partition functions data (use -DARTS_PARTFUN_DATA_PATH=/path/to/data/ to use other partition functions)")
endif (ARTS_PARTFUN_DATA_PATH)
message(STATUS "The partition functions data is at ${ARTS_PARTFUN_DATA_PATH}")

if (CMAKE_GENERATOR STREQUAL "Xcode")
  set (CTEST_MISC_OPTIONS -C \${CONFIGURATION})
endif (CMAKE_GENERATOR STREQUAL "Xcode")

add_custom_target(check-deps)

add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS}
    -E '\(python\\.|\\.pyarts\\.|\\.nocheck\\.|\\.slow\\.|\\.xmldata\\.|\\.planettoolbox\\.\)'
  DEPENDS check-deps)

add_custom_target(check-python
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS}
    -R '\(python\\.|\\.pyarts\\.|\\.cmdline\\.|.doc.\)'
    -E '\(\\.nocheck\\.|\\.slow\\.|\\.xmldata\\.|\\.planettoolbox\\.\)'
  DEPENDS check-deps python_tests)

add_custom_target(check-examples
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS}
    -R '\(\\pyarts.examples\\.\)'
    -E '\(\\.nocheck\\.\)'
  DEPENDS check-deps)

add_custom_target(check-xmldata
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS}
    -R '\(\\.xmldata\\.\)'
  DEPENDS check-deps)

add_dependencies(check-xmldata mkdir-arts-results)

add_custom_target(check-planettoolbox
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS}
    -R '\(\\.planettoolbox\\.\)'
  DEPENDS check-deps)

add_custom_target(check-all
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS}
    -E '\\.nocheck\\.'
  DEPENDS check-deps python_tests)

add_custom_target(check-all-controlfiles
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS}
    -E '^python\\.|\\.pyarts\\.|\\.nocheck\\.'
  DEPENDS check-deps)

add_custom_target(check-pyarts-only DEPENDS pyarts_tests)

add_dependencies(check-all mkdir-arts-results)


########### subdirectories ###############

add_subdirectory (3rdparty)

get_directory_property(ENABLE_RT4 DIRECTORY
                       3rdparty/polradtran/rt4
#                       3rdparty/rt4
                       DEFINITION ENABLE_RT4)

get_directory_property(ENABLE_FASTEM DIRECTORY
                       3rdparty/fastem
                       DEFINITION ENABLE_FASTEM)

get_directory_property(ENABLE_TMATRIX DIRECTORY
                       3rdparty/tmatrix
                       DEFINITION ENABLE_TMATRIX)

get_directory_property(ENABLE_WIGNER DIRECTORY
                       3rdparty/wigner
                       DEFINITION ENABLE_WIGNER)

if(FASTWIGNER)
  add_definitions( "-DFAST_WIGNER_PATH_3J=\"${ARTS_BINARY_DIR}/3rdparty/wigner/fastwigxj/fasttable.3j\"")
  add_definitions( "-DFAST_WIGNER_PATH_6J=\"${ARTS_BINARY_DIR}/3rdparty/wigner/fastwigxj/fasttable.6j\"")
endif()

if (OEM_SUPPORT)
  message (STATUS "OEM enabled")
else()
  message (STATUS "OEM disabled")
endif()

if (ENABLE_PCH)
  message (STATUS "Precompiled headers enabled")
  if (CCACHE_FOUND)
    message (STATUS "Make sure you set CCACHE_SLOPPINES properly for PCH:
   export CCACHE_SLOPPINESS=pch_defines,time_macros")
  endif()
else()
  message (STATUS "Precompiled headers disabled")
endif()

if (IPO_SUPPORTED)
  message(STATUS "IPO / LTO enabled")
else()
  message(STATUS "IPO / LTO disabled")
endif()

add_subdirectory (src)
add_subdirectory (doc)
add_subdirectory (controlfiles)
add_subdirectory (controlfiles-python)
add_subdirectory (python)
add_subdirectory (examples)
add_subdirectory (tests)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)

if (NOT "${NUMERIC}" STREQUAL "double")
  message (WARNING "You chose ${NUMERIC} as the Numeric type. Using anything different than double is completely unsupported, unrecommended and will probably lead to compilation errors, crashes and wrong results. Sorry, you're on your own here. Godspeed!")
endif()

message(STATUS "Building ARTS in ${CMAKE_BUILD_TYPE} mode")
if (CMAKE_BUILD_TYPE STREQUAL "Release" AND RELEASE_O2)
  message(STATUS "Optimization level set to O2 instead of O3 in Release mode")
endif ()
