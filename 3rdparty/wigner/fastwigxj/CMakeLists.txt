if(FASTWIGNER LESS 0)
  message (SEND_ERROR "ERROR:  FASTWIGNER must be positive if active")
  return()
endif()

check_c_source_compiles("
    #include <pthread.h>
    int main() {
      pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
      pthread_mutex_lock(&mutex);
      return 0;
    }
    "
    HAS_THREAD_MUTEX
)

check_c_source_compiles("
    #include <x86intrin.h>
    typedef long long v4di __attribute__ ((vector_size (32)));
    int main() {
    #if defined(__AVX2__)
        v4di a = { 1, 0., 0., 0. }, b = { 2, 0., 0., 0. };
        v4di c;
        c = _mm256_add_epi64(a,b);
    #endif
        return 0;
    }
    "
    HAS_AVX2
)

check_c_source_compiles([[
    int main() {
        asm __volatile__ ("    lfence  \\n\\t" : : : "memory");
        asm __volatile__ ("    sfence  \\n\\t" : : : "memory");
    }
    ]]
    HAS_LSFENCE
)

check_c_source_compiles("
    #include <x86intrin.h>
    typedef double v2df __attribute__ ((vector_size (16)));
    int main() {
        #if defined(__SSE4_1__)
            v2df x = { 1.2, 0.}, y = { 1.3, 0.}, c = { 0., 0.};
            v2df z;
            z = __builtin_ia32_blendvpd(y,x,c);
        #endif
    }
    "
    HAS_SSE4_1
)

check_c_source_compiles("
    #include <x86intrin.h>
    typedef long long int  v2di __attribute__ ((vector_size (16)));
    int main() {
        #if defined(__SSE4_2__)
            v2di a = { 1, 0. }, b = { 2, 0. };
            v2di c;
            c = (v2di) _mm_cmpgt_epi64((__m128i) a,(__m128i) b);
        #endif
    }
    "
    HAS_SSE4_2
)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated definitions file\n\n")

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_LONG_DOUBLE rule:\n")
if (HAS_LONG_DOUBLE)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_USE_LONG_DOUBLE 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_FLOAT128 rule:\n")
if (HAS_FLOAT128)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_USE_FLOAT128 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_THREAD rule:\n")
if (HAS_THREAD)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_HAVE_THREAD 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_THREAD_MUTEX rule:\n")
if (HAS_THREAD_MUTEX)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_HAVE_THREAD_MUTEX 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_LSFENCE rule:\n")
if (HAS_LSFENCE)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_HAVE_LSFENCE 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_SSE4_1 rule:\n")
if (HAS_SSE4_1)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_HAVE_SSE4_1 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_SSE4_2 rule:\n")
if (HAS_SSE4_2)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_HAVE_SSE4_2 1\n")
endif()

file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h "// auto-generated HAS_AVX2 rule:\n")
if (HAS_AVX2)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h
        "#define FASTWIGXJ_HAVE_AVX2 1\n")
endif()


add_library(fastwigxj STATIC
    src/fastwigxj.c
    src/canonicalise.c
    src/fortran_wrap.c
)

target_link_libraries(fastwigxj wigxjpf)

target_include_directories(fastwigxj
    PUBLIC 
        inc/
        ../wigxjpf/inc/
    PRIVATE
        cfg/
        src/
        ../wigxjpf/cfg/
        ../wigxjpf/src/
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/../wigxjpf/
)
    
set_target_properties(fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS} -march=native")

target_compile_options(fastwigxj PRIVATE -Wno-shadow)


add_library(hashing_fastwigxj STATIC
            src/unique_js.cc
            src/combine_js.cc
            src/gen_header.cc)

target_link_libraries(hashing_fastwigxj fastwigxj)

set_target_properties(hashing_fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS} -march=native")

add_executable(hash_fastwigxj src/hash_js.cc)

target_link_libraries(hash_fastwigxj hashing_fastwigxj fastwigxj)

set_target_properties(hash_fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS} -march=native")

add_custom_command(TARGET hash_fastwigxj POST_BUILD
  COMMAND "hash_fastwigxj" 
  ARGS "--max-E-3j=${FASTWIGNER}" "/dev/null" "./fasttable.3j"
  COMMAND "hash_fastwigxj" 
  ARGS "--max-E-6j=${FASTWIGNER}" "/dev/null" "./fasttable.6j"
  WORKING_DIRECTORY "${ARTS_BINARY_DIR}/3rdparty/wigner/fastwigxj/")
