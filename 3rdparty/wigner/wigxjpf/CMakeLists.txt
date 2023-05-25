include(CheckCSourceCompiles)

check_c_source_compiles("
    #include <math.h>
    int main() {
        long double a;
        long double b = 1.14, c = 2.00159;
        long double d;

        a = b + c;
        d = ldexpl(c, 5);
        return 0;
    }
    "
    WIGXJPF_IMPL_LONG_DOUBLE
)

check_c_source_compiles("
    #include \"quadmath.h\"
    int main() {
        __float128 a;
        __float128 b = 1.14, c = 2.00159;
        char s[64];

        a = b + c;

        quadmath_snprintf(s, sizeof(s), \"%Qf\", a);
        return 0;
    }
    "
    WIGXJPF_IMPL_FLOAT128
)

check_c_source_compiles("
    __thread int global = 0;
    int main() {
        global = 2;
        return 0;
    }
    "
    WIGXJPF_HAVE_THREAD
)

check_c_source_compiles("
    int main() {
        __int128 a;
        __uint128_t b;
    }
    "
    MULTI_WORD_INT_SIZEOF_ITEM
)


configure_file (${CMAKE_CURRENT_SOURCE_DIR}/wigxjpf_auto_config.h.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h)

add_library(wigxjpf STATIC
    src/prime_factor.c
    src/calc.c
    src/trivial_zero.c
    src/c_wrap.c
    src/fortran_wrap.c
    src/wigxjpf_error.c
)

target_include_directories(wigxjpf
    PUBLIC
        inc/
    PRIVATE
        cfg/
        src/
        ${CMAKE_CURRENT_BINARY_DIR}
)

target_compile_definitions(wigxjpf PUBLIC -DERRNO_ERROR_HANDLING)

set_target_properties(wigxjpf PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS}")

target_compile_options(wigxjpf PRIVATE -Wno-shadow)

if (IPO_SUPPORTED)
  set_property(TARGET wigxjpf PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
