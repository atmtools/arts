/* config.h.in.  Generated from configure.in by autoheader.  */

/* Compile Flags */
#cmakedefine COMPILE_FLAGS "${COMPILE_FLAGS}"

/* Compiler */
#cmakedefine COMPILER "${COMPILER}"

/* Define to compile with disort support */
#cmakedefine ENABLE_DISORT

/* Define to compile with NetCDF support */
#cmakedefine ENABLE_NETCDF

/* define if the compiler supports ISO C++ standard library */
#cmakedefine HAVE_STD 

/* check existence of c header files */
#cmakedefine HAVE_STDLIB_H 1
#cmakedefine HAVE_STRINGS_H 1
#cmakedefine HAVE_STRING_H 1
#cmakedefine HAVE_SYS_STAT_H 1
#cmakedefine HAVE_SYS_TIMES_H 1
#cmakedefine HAVE_SYS_TYPES_H 1
#cmakedefine HAVE_UNISTD_H 1

/* check existence of c++ header files */
#cmakedefine HAVE_CSTDLIB 1
#cmakedefine HAVE_CSTRING 1
#cmakedefine HAVE_SSTREAM 1
#cmakedefine HAVE_CTIME 1

/* availability of timer support */
#cmakedefine TIME_SUPPORT 1

/* Default Index type */
#cmakedefine INDEX ${INDEX}

/* Default Numeric type */
#cmakedefine NUMERIC ${NUMERIC}

/* Operating system name */
#cmakedefine OS_NAME "${OS_NAME}"

/* Define to compile with zlib support */
#define ENABLE_ZLIB 1

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

/* The size of `double', as computed by sizeof. */
#undef SIZEOF_DOUBLE

/* The size of `float', as computed by sizeof. */
#undef SIZEOF_FLOAT

/* The size of `int', as computed by sizeof. */
#undef SIZEOF_INT

/* The size of `long', as computed by sizeof. */
#undef SIZEOF_LONG

/* The size of `size_t', as computed by sizeof. */
#undef SIZEOF_SIZE_T

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Whether double precision is in use */
#define USE_DOUBLE

/* Whether float precision is in use */
/* #undef USE_FLOAT */

/* Version number of package */
#cmakedefine VERSION "${VERSION}"

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

