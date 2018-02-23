#define _GNU_SOURCE 1

/* Compile Flags */
#cmakedefine COMPILE_FLAGS "${COMPILE_FLAGS}"

/* Compiler */
#cmakedefine COMPILER "${COMPILER}"

/* Compiler */
#cmakedefine FORTRAN_COMPILER "${FORTRAN_COMPILER}"

/* Default directory for ARTS include files */
#cmakedefine ARTS_DEFAULT_INCLUDE_DIR "${ARTS_DEFAULT_INCLUDE_DIR}"

/* Define system constant */
#cmakedefine LINUX 1
#cmakedefine OSX 1
#cmakedefine WINDOWS 1

/* Threadprivate support */
#cmakedefine THREADPRIVATE_SUPPORTED

/* Define to compile with DISORT support */
#cmakedefine ENABLE_DISORT

/* Define to compile with RT4 support */
#cmakedefine ENABLE_RT4

/* Define to compile with FASTEM support */
#cmakedefine ENABLE_FASTEM

/* Define to compile with REFICE support */
#cmakedefine ENABLE_REFICE

/* Define to compile with relmat support */
#cmakedefine ENABLE_RELMAT

/* Define to compile with T-Matrix support */
#cmakedefine ENABLE_TMATRIX
#cmakedefine ENABLE_TMATRIX_QUAD

/* Define to compile with zlib support */
#cmakedefine ENABLE_ZLIB

/* Define to compile with documentation server support */
#cmakedefine ENABLE_DOCSERVER

/* Define to compile with NetCDF support */
#cmakedefine ENABLE_NETCDF ${NETCDF_FOUND}

/* Define to compile with FFTW support */
#cmakedefine ENABLE_FFTW ${FFTW_FOUND}

/* Define to compile with legacy HITRAN 2008 support */
#cmakedefine USE_HITRAN2008

/* define if the compiler supports ISO C++ standard library */
#cmakedefine HAVE_STD 

/* define if the compiler supports C++11 */
#cmakedefine CXX11_SUPPORT

/* define if OEM is enabled */
#cmakedefine OEM_SUPPORT

/* define if MPI was found */
#cmakedefine ENABLE_MPI

/* check existence of c header files */
#cmakedefine HAVE_STDLIB_H 1
#cmakedefine HAVE_STRINGS_H 1
#cmakedefine HAVE_STRING_H 1
#cmakedefine HAVE_SYS_STAT_H 1
#cmakedefine HAVE_SYS_TIMES_H 1
#cmakedefine HAVE_SYS_TYPES_H 1
#cmakedefine HAVE_UNISTD_H 1
#cmakedefine HAVE_GETOPT_H 1

/* check existence of c header files for libmicrohttpd */

#cmakedefine HAVE_FCNTL_H 1
#cmakedefine HAVE_MATH_H 1
#cmakedefine HAVE_ERRNO_H 1
#cmakedefine HAVE_LIMITS_H 1
#cmakedefine HAVE_STDIO_H 1
#cmakedefine HAVE_LOCALE_H 1
#cmakedefine HAVE_PTHREAD_H 1

/* check existence of additional c header files for libmicrohttpd */

#cmakedefine HAVE_SYS_TIME_H 1
#cmakedefine HAVE_SYS_MSG_H 1
#cmakedefine HAVE_NETDB_H 1
#cmakedefine HAVE_NETINET_IN_H 1
#cmakedefine HAVE_NETINET_TCP_H 1
#cmakedefine HAVE_SYS_SOCKET_H 1
#cmakedefine HAVE_SYS_MMAN_H 1
#cmakedefine HAVE_ARPA_INET_H 1
#cmakedefine HAVE_SYS_SELECT_H 1
#cmakedefine HAVE_POLL_H 1

/* check existence of c++ header files */
#cmakedefine HAVE_CSTDLIB 1
#cmakedefine HAVE_CSTRING 1
#cmakedefine HAVE_CTIME 1

/* check existence of functions */
#cmakedefine HAVE_REMOVE

/* availability of timer support */
#cmakedefine TIME_SUPPORT 1

/* Default Index type */
#cmakedefine INDEX ${INDEX}

/* Default Numeric type */
#cmakedefine NUMERIC ${NUMERIC}

/* Operating system name */
#cmakedefine OS_NAME "${OS_NAME}"

/* define if bool is a built-in type */
#define HAVE_BOOL 

/* define if the compiler supports const_cast<> */
#define HAVE_CONST_CAST 

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES 

/* define if the compiler supports static_cast<> */
#define HAVE_STATIC_CAST 

/* define if the compiler supports basic templates */
#define HAVE_TEMPLATES 

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
#undef NO_MINUS_C_MINUS_O

/* Operating system version */
#define OS_VERSION ""

/* Name of package */
#define PACKAGE "arts"

/* Define to the address where bug reports for this package should be sent. */
#undef PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#undef PACKAGE_NAME

/* Define to the full name and version of this package. */
#undef PACKAGE_STRING

/* Define to the one symbol short name of this package. */
#undef PACKAGE_TARNAME

/* Define to the version of this package. */
#undef PACKAGE_VERSION

/* The size of `int', as computed by sizeof. */
#cmakedefine SIZEOF_INT ${SIZEOF_INT}

/* The size of `long', as computed by sizeof. */
#cmakedefine SIZEOF_LONG ${SIZEOF_LONG}

/* The size of `size_t', as computed by sizeof. */
#cmakedefine SIZEOF_SIZE_T ${SIZEOF_SIZE_T}

/* The size of `double', as computed by sizeof. */
#cmakedefine SIZEOF_DOUBLE ${SIZEOF_DOUBLE}

/* The size of `float', as computed by sizeof. */
#cmakedefine SIZEOF_FLOAT ${SIZEOF_FLOAT}

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Whether double precision is in use */
#define USE_DOUBLE

/* Whether float precision is in use */
/* #undef USE_FLOAT */

/* Macro to ignore unused function parameters */
#define _U_ __attribute((unused))

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#undef inline
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

