#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([SHTns],[3.5.2],[],[shtns],[https://bitbucket.org/nschaeff/shtns])
AC_LANG([C])
AC_CONFIG_SRCDIR([sht_init.c])
AC_CONFIG_HEADERS([sht_config.h])
target="libshtns.a"		# by default, build the library.
objs="sht_fly.o"
install="install-lib"	# by default, install the library.

# optional variables :
AC_ARG_VAR(PYTHON, [the python interpreter (defaults to 'python')])
# optional features with --enable-XXX
AC_ARG_ENABLE([verbose],
	AS_HELP_STRING([--enable-verbose=0-3], [define verbosity level for the library: 0=silent, 1=default, 2=debug, 3=full]),
	[],[enable_verbose=1])
AC_ARG_ENABLE([openmp],
	AS_HELP_STRING([--enable-openmp], [Enable multi-threading with OpenMP]), [], [enable_openmp=default])
AC_ARG_ENABLE([mkl],
	AS_HELP_STRING([--enable-mkl], [Try to link with intel MKL instead of FFTW (can be slightly faster, but NOT THREAD SAFE)]), [], [enable_mkl=default])
AC_ARG_ENABLE([knl],
	AS_HELP_STRING([--enable-knl], [Enable compilation for Xeon Phi (KNL)]), [], [enable_knl=no])
AC_ARG_ENABLE([cuda],
	AS_HELP_STRING([--enable-cuda@<:@=kepler|pascal|volta|ampere@:>@], [Enable compilation for Nvidia gpu using cuda, optionally specifying architecture]), [], [enable_cuda=no])
AC_ARG_ENABLE([magic-layout],
	AS_HELP_STRING([--enable-magic-layout], [Compile specific version for the MagIC code]), [], [enable_magic_layout=no])
AC_ARG_ENABLE([python],
	AS_HELP_STRING([--enable-python], [Build Python interface]), [], [enable_python=no])
AC_ARG_ENABLE([long-double],
	AS_HELP_STRING([--enable-long-double], [Use long double during initialization for (maybe) slightly better precision]))
AC_ARG_ENABLE([f77],
	AS_HELP_STRING([--disable-f77], [Do not include F77 wrapper to call SHTns library from Fortran]))
AC_ARG_ENABLE([simd],
	AS_HELP_STRING([--disable-simd], [Do not use vector extensions (SSE2, AVX, AVX2, AVX-512, or VSX)]))
AC_ARG_ENABLE([ishioka],
	AS_HELP_STRING([--disable-ishioka], [Do not use the new and faster recursion]), [], [enable_ishioka=yes])
AC_ARG_ENABLE([kernel-compiler],
	AS_HELP_STRING([--enable-kernel-compiler=gcc], [Use a different compiler for kernels (default=gcc)]), [], [enable_kernel_compiler=gcc])


dnl Sanitize $prefix. Autoconf does this by itself, but so late in the
dnl generated configure script that the expansion does not occur until
dnl after our eval magic below.
AS_IF([test "$prefix" = "NONE"],[prefix=$ac_default_prefix])

# Checks for programs.
if test "x$CFLAGS" = "x"; then
	CFLAGS="-O2"
fi
CFLAGS="$CFLAGS -I$prefix/include -L$prefix/lib"
FFLAGS="$CFLAGS"
AC_PROG_CC
AC_PROG_SED
if test "$SED" = :; then
	AC_MSG_ERROR([sed program required.])
fi

# define macrao AX_CHECK_COMPILE_FLAG, and add it to CFLAGS
AC_DEFUN([AX_CHECK_COMPILE_FLAG],[
  AC_MSG_CHECKING([whether $CC accepts $1])
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $1"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
    [_AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags
     AC_MSG_RESULT([yes])
     m4_default([$2], [CFLAGS="$CFLAGS $1"])],
    [_AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags
     AC_MSG_RESULT([no])
     m4_default([$3], :)])
])dnl AX_CHECK_COMPILE_FLAGS

# define macrao AX_CHECK_COMPILE_FLAG, and add it to CFLAGS
AC_DEFUN([AX_CHECK_CUDA_ARCH],[
  AC_MSG_CHECKING([whether nvcc supports gpu with compute capability $1])
  AS_IF([nvcc sht_gpu.cu -c -arch=sm_$1 --dryrun > /dev/null 2>&1], [
	AC_MSG_RESULT(yes)
	nvcc_gencode_flags="$nvcc_gencode_flags -gencode=arch=compute_$1,code=sm_$1"
  ],[
	AC_MSG_RESULT(no)
  ])
])dnl AX_CHECK_CUDA_ARCH

############################################
### we can use two compilers. CC and CC2 ###
############################################
CC2=$enable_kernel_compiler
# backup current stuff
backup_cflags=$CFLAGS
backup_cc=$CC
CC=$CC2
#### setup the kernel compiler
AS_IF([test "x$enable_knl" != "xyes"], [
	# add gcc compile options if supported.
	AX_CHECK_COMPILE_FLAG([-march=native])
	# tune for skylake, it gives significantly better performance on Haswell with gcc 9.2
	AX_CHECK_COMPILE_FLAG([-mtune=skylake])
	# this is an icc compile option, try if it is supported:
	AX_CHECK_COMPILE_FLAG([-qopt-zmm-usage=high])
  ],[
	# KNL native cross-compiling
	AX_CHECK_COMPILE_FLAG([-march=knl],[], [
	    CC=icc
	    AX_CHECK_COMPILE_FLAG([-xMIC-AVX512],[CC="$CC -xMIC-AVX512"], [AC_MSG_ERROR(["Xeon Phi (KNL) not supported by $CC. Use a recent gcc or icc compiler."])])
	])
])
CFLAGS="$CFLAGS -O3"
AX_CHECK_COMPILE_FLAG([-ffast-math])
AX_CHECK_COMPILE_FLAG([-std=gnu99])
AX_CHECK_COMPILE_FLAG([-fopenmp-simd])
AX_CHECK_COMPILE_FLAG([-qopenmp-simd])
# the following prevents gcc from inserting memset in tight loops for zeroing stuff. Must appear after -O3
AX_CHECK_COMPILE_FLAG([-fno-tree-loop-distribute-patterns])
# avoid insering VZEROUPPER instructions everywhere.
AX_CHECK_COMPILE_FLAG([-mno-vzeroupper])

## check for nasty bug in binutils, apply possible workaround. See https://bitbucket.org/nschaeff/shtns/issues/37/
bug_as=`echo 'asm("vmovaps 64(,%rax), %zmm0");' > x.c ; $CC -c x.c ; objdump -dw x.o | grep -c '1('`
rm -f x.c x.o
AS_IF([test "x$bug_as" == "x1"], [
    AX_CHECK_COMPILE_FLAG([-fno-tree-pre])
])

CC2=$CC
CFLAGS2=$CFLAGS
# restore the standard compiler
CC=$backup_cc
CFLAGS=$backup_cflags
################################
### kernel compiler is setup ###
################################

AS_IF([test "x$enable_knl" != "xyes"], [
	# add gcc compile options if supported.
	AX_CHECK_COMPILE_FLAG([-march=native])
	# tune for skylake, it gives significantly better performance on Haswell with gcc 9.2
	AX_CHECK_COMPILE_FLAG([-mtune=skylake])
	# this is an icc compile option, try if it is supported:
	AX_CHECK_COMPILE_FLAG([-qopt-zmm-usage=high])
  ],[
	target="libshtns_mic.a"
	# KNL native cross-compiling
	AX_CHECK_COMPILE_FLAG([-march=knl],[], [
	    CC=icc
	    AX_CHECK_COMPILE_FLAG([-xMIC-AVX512],[CC="$CC -xMIC-AVX512"], [AC_MSG_ERROR(["Xeon Phi (KNL) not supported by $CC. Use a recent gcc or icc compiler."])])
	])
	# With KNL, enable mkl and openmp by default.
	AS_IF([test "x$enable_mkl" = "xdefault"], [enable_mkl=yes])
	AS_IF([test "x$enable_openmp" = "xdefault"], [enable_openmp=yes])
])
AX_CHECK_COMPILE_FLAG([-ffast-math])
AX_CHECK_COMPILE_FLAG([-std=gnu99])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h stdio.h string.h math.h complex.h])

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,AC_MSG_ERROR([math library not found.]))

# With Python, enable openmp by default.
AS_IF([test "x$enable_python" != "xno"], [
	AS_IF([test "x$enable_openmp" = "xdefault"], [enable_openmp=yes])
])

# for MagIC code, also disable matrix transforms.
AS_IF([test "x$enable_magic_layout" != "xno"], [
	AC_DEFINE([SHTNS4MAGIC],[1],[I need the transforms compatible with the MagIC code, to speed it up!])
])

nvcc_flags=""
# Disable precomputed matrix support with openmp: to avoid bugs arising with icc + openmp.
AS_IF([test "x$enable_openmp" = "xyes"], [
	nvcc_flags="-Xcompiler -fopenmp"
	AC_OPENMP
	CFLAGS="$CFLAGS $OPENMP_CFLAGS"
	dnl AC_CHECK_HEADERS([omp.h])
	objs="$objs sht_omp.o"
	AS_IF([test "x$enable_knl" != "xno"], [
		target="libshtns_mic_omp.a"
	],[
		target="libshtns_omp.a"
	])
])

# Check for CUDA
AS_IF([test "x$enable_cuda" != "xno"], [
	CFLAGS="$CFLAGS -I$CUDA_PATH/include"
	CFLAGS2="$CFLAGS2 -I$CUDA_PATH/include"
	# assume 64 bit mode, add default directory to find CUDA
	LDFLAGS="$LDFLAGS -L$CUDA_PATH/lib64"
	echo $LDFLAGS
	AC_CHECK_LIB([cudart],[cudaMalloc],,AC_MSG_ERROR([cudart library not found. Try adding LDFLAGS="-L/usr/local/cuda/lib64"]))
	AC_CHECK_LIB([cufft],[cufftPlanMany],,AC_MSG_ERROR([cufft library not found.]))
	objs="$objs sht_gpu.o"
	target="libshtns_cuda.a"

	nvcc_gencode_flags=""
	# find cuda arch:
	AS_IF([test "x$enable_cuda" == "xyes"],[
	    AX_CHECK_CUDA_ARCH([30])
	    AX_CHECK_CUDA_ARCH([60])
	    AX_CHECK_CUDA_ARCH([70])
	    AX_CHECK_CUDA_ARCH([80])
	],[
	    AS_CASE([$enable_cuda], "kepler", [cuda_cc=30], "pascal", [cuda_cc=60], "volta", [cuda_cc=70], 'ampere', [cuda_cc=80], [cuda_cc=$enable_cuda])
	    AX_CHECK_CUDA_ARCH([$cuda_cc])
	])
	AS_IF([test "x$nvcc_gencode_flags" == "x"],
	    AC_MSG_ERROR([nvcc does not support any compatible gpu. Is nvcc available ?])
	])
	nvcc_flags="$nvcc_flags $nvcc_gencode_flags"
	LIBS="$LIBS -lstdc++"	# cuda uses the c++ standard library. This is needed to link with gcc
])
libname=$target		# set the libname befor checking for python.

# Check for FFTW (MKL preferred, FFTW optional)
AS_IF([test "x$enable_mkl" == "xyes"], [
  mkl_found=no
  AX_CHECK_COMPILE_FLAG([-mkl], [
	LIBS_bak="$LIBS"
	LIBS="-mkl $LIBS"	# easy way, for intel compiler !
	AC_CHECK_FUNC(fftw_plan_many_dft, [
	    AC_MSG_NOTICE(["FFTW interface found in MKL, link with -mkl"])
	    mkl_found=yes
	],[
	    LIBS="$LIBS_bak"
	    AC_MSG_NOTICE(["FFTW interface not found in MKL libraries (-mkl)."])
	])
  ])
  AS_IF([test "x$mkl_found" == "xno"], [
	# assume 64 bit mode, add default directory to find the MKL
	LDFLAGS="$LDFLAGS -L$MKLROOT/lib/intel64"
	AS_IF([test "x$enable_openmp" = "xyes"], [
		# multi-threaded MKL libs
		#mkl_libs="-Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread"
		mkl_libs="-Wl,--start-group -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group"
	],[
		# sequential MKL libs
		mkl_libs="-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group"
	])
	AC_CHECK_LIB(mkl_intel_lp64, fftw_plan_many_dft, [
	    LIBS="$mkl_libs $LIBS"
	    AC_MSG_NOTICE(["FFTW interface found in MKL, link with $mkl_libs"])
	    mkl_found=yes
	],[
	    AC_MSG_NOTICE(["FFTW interface not found in MKL libraries ($mkl_libs). Try to add 'LDFLAGS=-L/path/to/mkl' to ./configure command line."])
	],["$mkl_libs"])
  ])
  AS_IF([test "x$mkl_found" == "xno"], [
    AS_IF([test "x$enable_mkl" == "xyes"], [AC_MSG_ERROR(["FFTW interface not found in MKL."])], [enable_mkl=no])
  ],[
    AC_DEFINE(HAVE_LIBFFTW3_OMP, 1, [MKL has the multi-thread fftw interface.])
  ])
])

# Check for FFTW, if MKL was not required and not found.
AS_IF([test "x$enable_mkl" != "xyes"], [
	# fftw3
	AC_CHECK_LIB([fftw3],[fftw_plan_many_dft],,AC_MSG_ERROR([FFTW3 library required. Or try using the intel MKL library with --enable-mkl]))
	AS_IF([test "x$enable_openmp" = "xyes"], [
		AC_CHECK_LIB([fftw3_omp], [fftw_init_threads],,AC_MSG_ERROR([FFTW3 does not support OpenMP. Did you compile it with --enable-openmp ?.]))
	])
])


# Checks related to Python and NumPy paths:
AS_IF([test "x$enable_python" != "xno"], [
	AS_IF([test "x$PYTHON" = "x"], [
		AC_MSG_CHECKING(for python with numpy package)
		py_test="python python2 python3"
	],[
		AC_MSG_CHECKING(if $PYTHON has numpy package)
		py_test="$PYTHON"
	])
	for py in $py_test
	do :
		numpy_inc=`$py -c "from numpy import get_include; print(get_include())" 2>/dev/null`
		AS_IF([test "x$numpy_inc" != "x"],[break])
	done
	AS_IF([test "x$numpy_inc" = "x"], [
		AC_MSG_RESULT(no)
		AC_MSG_ERROR([NumPy package is required for the python interface.])
	],[
		AC_MSG_RESULT($py)
		PYTHON=$py
	])
	python_inc=`$PYTHON -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())"  2>/dev/null`
	echo "  python include path = $python_inc"
	echo "  numpy include path = $numpy_inc"
	AC_SUBST([numpy_inc])
	AC_SUBST([python_inc])
	CFLAGS="$CFLAGS -fpic"	# required compile flag for python extensions.
	CFLAGS2="$CFLAGS2 -fpic"	# required compile flag for python extensions.
	target="_shtns.so"	# build python extension instead of C library.
	install="install-py"	# install python extension instead of C library.
])

# Checks related to long double
AS_IF([test "x$enable_long_double" == "xyes"], [
	AC_TYPE_LONG_DOUBLE_WIDER])

# Disable Fortran interface ?
AS_IF([test "x$enable_f77" != "xno"], [
	AC_DEFINE([SHT_F77_API],[1],[Compile the Fortran API])
])

# Disable SIMD ?
AS_IF([test "x$enable_simd" != "xno"], [
	AC_DEFINE([_GCC_VEC_],[1],[I compile with GCC 4 or ICC 14 or later, and I would like fast vectorized code (if SSE2, AVX or MIC is supported) !])
])


AS_IF([test "x$enable_ishioka" == "xyes"], [
       AC_DEFINE([SHTNS_ISHIOKA],[1],[Enable the new recurrence proposed by Ishioka (2018) see https://doi.org/10.2151/jmsj.2018-019 (faster)])
])

# Verbosity setting
AS_IF([test "x$enable_verbose" == "xno"], [enable_verbose=0],
	[test "x$enable_verbose" == "xyes"], [enable_verbose=1])
AC_DEFINE_UNQUOTED([SHT_VERBOSE],$enable_verbose,[0:no output, 1:output info to stdout, 2:more output (debug info), 3:also print fftw plans.])
echo "  verbose level = $enable_verbose"

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([fftw_cost])

### for cycle.h ###
   AC_C_INLINE
   AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])

   AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])

   AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])

   dnl Cray UNICOS _rtc() (real-time clock) intrinsic
   AC_MSG_CHECKING([for _rtc intrinsic])
   rtc_ok=yes
   AC_LINK_IFELSE([AC_LANG_PROGRAM([[#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif]], [[_rtc()]])],[AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])],[rtc_ok=no])
   AC_MSG_RESULT($rtc_ok)
### end cycle.h ###

#### for the Fortran example ####
AC_PROG_FC
#### end Fortran #####

echo "prefix=$prefix"
AC_SUBST([target])
AC_SUBST([install])
AC_SUBST([objs])
AC_SUBST([libname])
AC_SUBST([CC2])
AC_SUBST([CFLAGS2])
AC_SUBST([FFLAGS])
AC_SUBST([nvcc_flags])
AC_CONFIG_FILES([Makefile setup.py])
AC_OUTPUT

## if binutils bug: warn the user.
AS_IF([test "x$bug_as" == "x1"], [
    AC_MSG_NOTICE(["WARNING! your version of binutils is buggy. Please update the binutils package to v2.32 or more recent.
	Meanwhile, I try to apply a workaround (no guarantee, possible performance loss).
	See https://bitbucket.org/nschaeff/shtns/issues/37/"])
])

