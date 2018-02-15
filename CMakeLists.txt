
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

cmake_minimum_required (VERSION 3.1.0)
set (CMAKE_LEGACY_CYGWIN_WIN32 0)

if (NOT ENABLE_ALLPATHS)
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

project (ARTS C CXX)
enable_language(C)
enable_language(CXX)

enable_testing ()

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
      message (FATAL_ERROR "Gfortran does not support quad-precision. Only supported on Intel compiler.")
    endif()
  elseif (FORTRAN_COMPILER MATCHES "ifort.*")
    set (FORTRAN_EXTRA_FLAGS "-i8 -r8")
  else ()
    message (FATAL_ERROR "Unsupported Fortran compiler. Must be gfortran or ifort.")
  endif (FORTRAN_COMPILER MATCHES "gfortran.*")

  INCLUDE (CheckTypeSize)
  CHECK_TYPE_SIZE (long SIZEOF_LONG)
  if (NOT SIZEOF_LONG EQUAL 8)
    message (FATAL_ERROR "Fortran interface only works on 64-bit with long size of 8 bytes.\n"
      "Size of long on this system: ${SIZEOF_LONG} bytes")
  endif (NOT SIZEOF_LONG EQUAL 8)

endif (ENABLE_FORTRAN)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)

include (CheckIncludeFile)
include (CheckIncludeFileCXX)

include (CheckTypeSize)
include (CheckFunctionExists)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

include (ArtsTestcases)

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

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -fvisibility-inlines-hidden")

########### ccache Support ##########

if (ENABLE_CCACHE)
  find_program(CCACHE_FOUND ccache)
  if (CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif (CCACHE_FOUND)
endif ()

########### C++11 Support ##########

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Explicitly set C++ flags for Intel compiler because cmake
# only supports it from version 3.6.0+. We don't want to
# depend on this newer version of cmake because it was only
# released in June 2016.
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

########### Intel compiler specific settings ##########

# Turn off value-unsafe floating-point optimizations which lead
# to large descrepancies compared to other compilers
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fp-model precise")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model precise")
endif ()

if (ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise")
endif ()

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

########### C API ##########

if (ENABLE_C_API)
  set (C_API true)
  add_definitions(-fPIC)
endif ()

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

######## Enable timer support, if we have all necessary header files #######
if (HAVE_SYS_TYPES_H AND HAVE_UNISTD_H AND HAVE_SYS_STAT_H AND HAVE_CTIME
    AND HAVE_SYS_TIMES_H)
  set (TIME_SUPPORT 1)
endif (HAVE_SYS_TYPES_H AND HAVE_UNISTD_H AND HAVE_SYS_STAT_H AND HAVE_CTIME
  AND HAVE_SYS_TIMES_H)

########### Zlib support ##########
include (FindZLIB)
set (ENABLE_ZLIB ${ZLIB_FOUND})

########## BLAS support #########
if (${CMAKE_VERSION} VERSION_LESS "2.8.8")
  # Workaround: Older versions of cmake
  # rely on fortran being enabled to detect BLAS
  enable_language (Fortran)
endif()

########## BLAS/LAPACK support #########

find_package (LAPACK REQUIRED)

########## NetCDF support #########
if (ENABLE_NETCDF)
   include (FindNetCDF)
endif (ENABLE_NETCDF)

########## OpenMP support #########
if (NOT NO_OPENMP AND NOT CMAKE_GENERATOR STREQUAL "Xcode")
  include (FindOpenMP)
endif (NOT NO_OPENMP AND NOT CMAKE_GENERATOR STREQUAL "Xcode")

########### Docserver support ##########

if (NOT NO_DOCSERVER)
  set (ENABLE_DOCSERVER true)
endif (NOT NO_DOCSERVER)

########### Remove NDEBUG as set by cmake >2.8.10 ##########
if (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  string (REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELWITHDEBINFO
          "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  string (REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO
          "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")

########### Special flags for certain compilers ##########
if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-1024")
endif (CMAKE_COMPILER_IS_GNUCXX)

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
  ARTS_ADD_COMPILER_FLAG (Wshadow)
  ARTS_ADD_COMPILER_FLAG (Wconversion)
  ARTS_ADD_COMPILER_FLAG (Wno-sign-conversion)
  ARTS_ADD_COMPILER_FLAG (Wno-unknown-pragmas)
  ARTS_ADD_COMPILER_FLAG (Wno-return-type-c-linkage)
  ARTS_ADD_COMPILER_FLAG (Wno-strict-overflow)
endif (CMAKE_GENERATOR STREQUAL "Xcode")

if(CMAKE_BUILD_TYPE STREQUAL "Native")
  ARTS_ADD_COMPILER_FLAG (march=native)
  ARTS_ADD_COMPILER_FLAG (fcx-fortran-rules)
  ARTS_ADD_COMPILER_FLAG (O3)
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  ARTS_ADD_COMPILER_FLAG (w1)
  ARTS_ADD_COMPILER_FLAG (diag-disable=383,980,981,1418,1419,1682)
endif ()

if (WERROR)
  ARTS_ADD_COMPILER_FLAG (Werror)
endif (WERROR)

if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
  set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
endif (OPENMP_FOUND)

set (COMPILER "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER})")

########### Xcode specific settings ##########
if (CMAKE_GENERATOR STREQUAL "Xcode")
  set (COMPILER "Xcode")
  set (CMAKE_OSX_SYSROOT "")
  set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/build CACHE PATH
       "Single output directory for building all libraries.")
  set (EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/build CACHE PATH
       "Single output directory for building all executables.")
  mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

  ######### Enable C++11 and link against the new libc++ ##########
  if (WITH_XCODE_LIBCPP)
    set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set (CMAKE_XCODE_ATTRIBUTE_OBJROOT ${CMAKE_BINARY_DIR}/build/\${CONFIGURATION})
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

set (ARTS_DEFAULT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/controlfiles")

if (WITH_HITRAN2008)
  set (USE_HITRAN2008 true)
endif (WITH_HITRAN2008)

string (TOUPPER ${CMAKE_BUILD_TYPE} MY_BUILD_TYPE)
set (COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${MY_BUILD_TYPE}}")

add_custom_command (
  OUTPUT auto_version.txt
  COMMAND ${CMAKE_COMMAND} -D "ARTS_SOURCE_DIR:STRING=\"${CMAKE_SOURCE_DIR}\"" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/update_auto_version_txt.cmake
  DEPENDS ChangeLog
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
  get_filename_component (ARTS_XML_DATA_DIR "../arts-xml-data" ABSOLUTE)
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

if (CMAKE_GENERATOR STREQUAL "Xcode")
  set (CTEST_MISC_OPTIONS -C \${CONFIGURATION})
endif (CMAKE_GENERATOR STREQUAL "Xcode")

add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS} -E '\(.nocheck.|.slow.|.xmldata.|.classroom.|.planettoolbox.\)' DEPENDS arts)

add_custom_target(check-xmldata
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS} -R '\(.xmldata.|.classroom.\)' DEPENDS arts)

add_dependencies(check-xmldata mkdir-arts-results)

add_custom_target(check-classroom
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS} -R '\(.classroom.\)' DEPENDS arts)

add_dependencies(check-classroom mkdir-arts-results)

add_custom_target(check-planettoolbox
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS} -R '\(.planettoolbox.\)' DEPENDS arts)

add_custom_target(check-all
  COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_MISC_OPTIONS} --output-on-failure ${CTEST_JOBS} -E '\(.nocheck.\)' DEPENDS arts)

add_dependencies(check-all mkdir-arts-results)


########### subdirectories ###############

add_subdirectory (3rdparty)

get_directory_property(ENABLE_DISORT DIRECTORY
                       3rdparty/disort1.2
                       DEFINITION ENABLE_DISORT)

get_directory_property(ENABLE_RT4 DIRECTORY
                       3rdparty/polradtran/rt4
#                       3rdparty/rt4
                       DEFINITION ENABLE_RT4)

get_directory_property(ENABLE_FASTEM DIRECTORY
                       3rdparty/fastem
                       DEFINITION ENABLE_FASTEM)

get_directory_property(ENABLE_REFICE DIRECTORY
                       3rdparty/refice
                       DEFINITION ENABLE_REFICE)

get_directory_property(ENABLE_TMATRIX DIRECTORY
                       3rdparty/tmatrix
                       DEFINITION ENABLE_TMATRIX)

get_directory_property(ENABLE_RELMAT DIRECTORY
                       3rdparty/relmat
                       DEFINITION ENABLE_RELMAT)

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

if (C_API)
  message (STATUS "C API enabled")
else()
  message (STATUS "C API disabled (use -DENABLE_C_API=1 to enable)")
endif()

add_subdirectory (src)
add_subdirectory (doc)
add_subdirectory (controlfiles)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)

if (NOT "${NUMERIC}" STREQUAL "double")
  message (WARNING "You chose ${NUMERIC} as the Numeric type. Using anything different than double is completely unsupported, unrecommended and will probably lead to compilation errors, crashes and wrong results. Sorry, you're on your own here. Godspeed!")
endif()

