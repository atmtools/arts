include(CheckCSourceCompiles)

check_c_source_compiles("
    #include <pthread.h>
    int main() {
      pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
      pthread_mutex_lock(&mutex);
      return 0;
    }
    "
    FASTWIGXJ_HAVE_THREAD_MUTEX
)

check_c_source_compiles("
    __thread int global = 0;
    int main() {
        global = 2;
        return 0;
    }
    "
    FASTWIGXJ_HAVE_THREAD
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
    FASTWIGXJ_HAVE_AVX2
)

check_c_source_compiles([[
    int main() {
        asm __volatile__ ("    lfence  \\n\\t" : : : "memory");
        asm __volatile__ ("    sfence  \\n\\t" : : : "memory");
    }
    ]]
    FASTWIGXJ_HAVE_LSFENCE
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
    FASTWIGXJ_HAVE_SSE4_1
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
    FASTWIGXJ_HAVE_SSE4_2
)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/fastwigxj_auto_config.h.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/fastwigxj_auto_config.h)

add_library(fastwigxj STATIC
    src/fastwigxj.c
    src/canonicalise.c
    src/fortran_wrap.c
)

target_include_directories(fastwigxj
    PUBLIC 
        inc/
        cfg/
        ${CMAKE_CURRENT_BINARY_DIR}
    PRIVATE
        src/
)

target_link_libraries(fastwigxj wigxjpf)
    
set_target_properties(fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS} -march=native")

if (NOT CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_options(fastwigxj PRIVATE -Wno-shadow)
endif()


add_library(hashing_fastwigxj STATIC
    src/fastwigxj.c
    src/canonicalise.c
    src/fortran_wrap.c
    src/unique_js.cc
    src/combine_js.cc
    src/gen_header.cc
)

target_include_directories(hashing_fastwigxj
    PUBLIC
        inc/
        cfg/
        src/
        ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(hashing_fastwigxj wigxjpf)

set_target_properties(hashing_fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS} -march=native")


add_executable(hash_fastwigxj src/hash_js.cc)

target_link_libraries(hash_fastwigxj hashing_fastwigxj)

set_target_properties(hash_fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS} -march=native")

if(FASTWIGNER LESS 1)
  message (SEND_ERROR "ERROR:  FASTWIGNER must be positive if active")
  return()
endif()

add_custom_command(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/fasttable.3j
    COMMAND
        hash_fastwigxj --max-E-3j=${FASTWIGNER} /dev/null ${CMAKE_CURRENT_BINARY_DIR}/fasttable.3j
    DEPENDS
        hash_fastwigxj
)

add_custom_command(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/fasttable.6j
    COMMAND
        hash_fastwigxj --max-E-6j=${FASTWIGNER} /dev/null ${CMAKE_CURRENT_BINARY_DIR}/fasttable.6j
    DEPENDS
        hash_fastwigxj
)

add_custom_target(fastwigxj_generate_tables
    DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/fasttable.3j
        ${CMAKE_CURRENT_BINARY_DIR}/fasttable.6j
)

add_dependencies(fastwigxj fastwigxj_generate_tables)

if (IPO_SUPPORTED)
  set_property(TARGET fastwigxj PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET hashing_fastwigxj PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
