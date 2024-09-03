add_subdirectory(SHT)

include(CheckSymbolExists)
check_include_file(complex.h HAVE_COMPLEX_H)
check_include_file(c_asm.h HAVE_C_ASM_H)
check_include_file(intrinsics.h HAVE_INTRINSICS_H)
check_include_file(inttypes.h HAVE_INTTYPES_H)
check_symbol_exists(gethrtime sys/time.h HAVE_GETHRTIME)
check_symbol_exists(mach_absolute_time mach/mach_time.h HAVE_MACH_ABSOLUTE_TIME)
check_include_file(mach/mach_time.h HAVE_MACH_MACH_TIME_H)
check_symbol_exists(hrtime_t time.h HAVE_HRTIME_T)
check_symbol_exists(clock_gettime time.h HAVE_CLOCK_GETTIME)
find_library(HAVE_LIBM NAMES m)
if (FFTW_LIBRARY)
  set(HAVE_LIBFFTW3 1)
endif()
set(PACKAGE_NAME "SHTns")
set(PACKAGE_VERSION "3.5.2")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_URL "https://bitbucket.org/nschaeff/shtns")
set(SHTNS_ISHIOKA 1)
set(STDC_HEADERS 1)
set(SHT_VERBOSE 0)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/sht_config_cmake.h.in
                ${CMAKE_CURRENT_BINARY_DIR}/sht_config.h)

add_library(
  shtns
  sht_init.c
  sht_kernels_a.c
  sht_kernels_s.c
  sht_odd_nlat.c
  sht_fly.c
  sht_omp.c
  )

Check_C_Compiler_Flag(-fopenmp-simd CCFLAG_OPENMP_SIMD)
if(CCFLAG_OPENMP_SIMD)
  target_compile_options(shtns PRIVATE -fopenmp-simd)
endif()

Check_C_Compiler_Flag(-mno-vzeroupper CCFLAG_NO_VZEROUPPER)
if(CCFLAG_NO_VZEROUPPER)
  target_compile_options(shtns PRIVATE -mno-vzeroupper)
endif()

target_compile_options(
  shtns
  PRIVATE
  -ffast-math
  -Wno-unused -Wno-shadow -Wno-unused-parameter
  -Wno-conversion -Wno-sign-compare -Wno-switch -Wno-misleading-indentation
  -Wno-uninitialized -Wno-implicit-fallthrough -Wno-unused-command-line-argument
)
target_compile_definitions(shtns PRIVATE _GNU_SOURCE SHTNS_GIT=\"ARTS\")
target_include_directories(shtns PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_dependencies(shtns shtns_generated_files)

