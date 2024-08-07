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
    HAS_LONG_DOUBLE
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
    HAS_FLOAT128
)

check_c_source_compiles("
    __thread int global = 0;
    int main() {
        global = 2;
        return 0;
    }
    "
    HAS_THREAD
)

check_c_source_compiles("
    int main() {
        __int128 a;
        __uint128_t b;
    }
    "
    HAS_UINT128
)

check_c_source_compiles("
    struct trail_struct {
        int a;
        int trail[];
    };
    int main() {
        struct trail_struct ts;
    }
    "
    TEST_TRAIL_ARRAY_EMPTY
)

check_c_source_compiles("
    struct trail_struct {
        int a;
        int trail[0];
    };
    int main() {
        struct trail_struct ts;
    }
    "
    TEST_TRAIL_ARRAY_ZERO
)

check_c_source_compiles("
    struct trail_struct {
        int a;
        int trail[1];
    };
    int main() {
        struct trail_struct ts;
    }
    "
    TEST_TRAIL_ARRAY_ONE
)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated definitions file\n\n")

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated HAS_LONG_DOUBLE rule:\n")
if (HAS_LONG_DOUBLE)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
        "#define WIGXJPF_IMPL_LONG_DOUBLE 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated HAS_FLOAT128 rule:\n")
if (HAS_FLOAT128)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
        "#define WIGXJPF_IMPL_FLOAT128 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated HAVE_THREAD rule:\n")
if (HAS_THREAD)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
        "#define WIGXJPF_HAVE_THREAD 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated HAS_UINT128 rule:\n")
if (HAS_UINT128)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
        "#define MULTI_WORD_INT_SIZEOF_ITEM 8\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated TRAIL_ARRAY_EMPTY rule:\n")
if (TEST_TRAIL_ARRAY_EMPTY)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
        "#define TRAIL_ARRAY_EMPTY 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated TRAIL_ARRAY_ZERO rule:\n")
if (TEST_TRAIL_ARRAY_ZERO)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
        "#define TRAIL_ARRAY_ZERO 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
  "// auto-generated TRAIL_ARRAY_ONE rule:\n")
if (TEST_TRAIL_ARRAY_ONE)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/wigxjpf_auto_config.h
        "#define TRAIL_ARRAY_ONE 1\n")
endif()

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
        cfg/
        ${CMAKE_CURRENT_BINARY_DIR}
    PRIVATE
        src/
)

target_compile_definitions(wigxjpf PUBLIC -DERRNO_ERROR_HANDLING)

set_target_properties(wigxjpf PROPERTIES COMPILE_FLAGS
  "-O3 -fPIC ${CMAKE_C_FLAGS}")

if (NOT CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  target_compile_options(wigxjpf PRIVATE -Wno-shadow)
endif()
